//! Zone air predictor-corrector, surface conduction, and ideal-loads system.
//!
//! The zone temperature update is the EnergyPlus third-order scheme from
//! `ZoneTempPredictorCorrector.cc`:
//!
//! ```text
//! ZT = (TempIndCoef + Q_sys + AirPowerCap * (3*ZTM1 - 3/2*ZTM2 + 1/3*ZTM3))
//!      / (11/6 * AirPowerCap + TempDepCoef)
//! ```
//!
//! Ideal-loads capacity is the sensible load that holds the thermostat setpoint,
//! which is `tempDepLoad * Tset - tempIndLoad` in that same file. Convective
//! electric baseboard meets the remaining heating-setpoint load in heating
//! sequence order before or after that ideal-loads unit. Opaque conduction
//! uses the state-space conduction transfer functions from
//! `Construction::calculateTransferFunction`. Interior convection
//! is TARP. Exterior convection is the DOE-2 model. Opaque surfaces exchange
//! long-wave radiation through Hottel's Script F coefficients. Zone humidity
//! uses the same third-order corrector as EnergyPlus.

use crate::baseboard::convective_electric;
use crate::convection::{doe2_exterior_hc, tarp_interior_hc};
use crate::ctf::{conduction_transfer_functions, Ctf, CtfLayer};
use crate::epw::{WeatherFile, WeatherHour};
use crate::model::{
    infiltration_modifier, Boundary, Gain, InfilDensity, Model, Surface, SurfaceKind, ZoneEquipment,
};
use crate::psychrometrics::{
    psy_cp_air_fn_w, psy_hg_air_fn_w_tdb, psy_rho_air_fn_pb_tdb_w, psy_w_fn_tdb_rh_pb, PsychState,
};
use crate::radiosity::{class_from_tilt, net_longwave, tilt_from_normal, EnclosureSurface};
use crate::simple_glazing::{
    absorb_transmitted, simple_glazing_beam, window_face_temperatures, SolarReceiver,
};
use crate::solar::{daily_solar_coeffs, day_of_year, sun4, SunVector};
use crate::{KELVIN, SIGMA};

#[derive(Clone, Debug)]
pub struct SimOptions {
    pub warmup_days: u32,
}

impl Default for SimOptions {
    fn default() -> Self {
        Self { warmup_days: 7 }
    }
}

#[derive(Clone, Debug)]
pub struct ZoneStep {
    pub name: String,
    pub temperature: f64,
    pub humidity_ratio: f64,
    /// Sensible heating delivered to the zone, including baseboard heat.
    pub heating_w: f64,
    pub cooling_w: f64,
    /// Ideal-loads heating only. Baseboard heat is reported on [`BaseboardStep`].
    pub ideal_heating_w: f64,
    pub baseboards: Vec<BaseboardStep>,
}

#[derive(Clone, Debug)]
pub struct BaseboardStep {
    pub name: String,
    pub heating_w: f64,
    pub electric_w: f64,
}

#[derive(Clone, Debug)]
pub struct TimestepResult {
    pub year: i32,
    pub month: u32,
    pub day: u32,
    pub hour: u32,
    pub end_minute: u32,
    pub outdoor_db: f64,
    pub zones: Vec<ZoneStep>,
}

#[derive(Clone, Debug)]
pub struct Simulation {
    pub steps: Vec<TimestepResult>,
    pub heating_j: f64,
    pub cooling_j: f64,
}

struct ZoneRuntime {
    temperature: f64,
    humidity: f64,
    history: [f64; 3],
    humidity_history: [f64; 3],
}

struct SurfaceRuntime {
    ctf: Option<Ctf>,
    t_out: f64,
    t_in: f64,
    t_out_hist: Vec<f64>,
    t_in_hist: Vec<f64>,
    q_out_hist: Vec<f64>,
    q_in_hist: Vec<f64>,
    /// Zone timesteps remaining before the next CTF history update.
    steps_until_update: usize,
}

const SURFACE_DAMP: f64 = 5.0;

pub fn simulate(
    model: &Model,
    weather: &WeatherFile,
    options: &SimOptions,
) -> Result<Simulation, String> {
    if model.zones.is_empty() {
        return Err("model has no zones".to_string());
    }
    let location = model
        .location
        .clone()
        .unwrap_or_else(|| crate::model::Location {
            name: weather.city.clone(),
            latitude: weather.latitude,
            longitude: weather.longitude,
            time_zone: weather.time_zone,
            elevation: weather.elevation,
        });
    let (start, end) = resolve_period(model, weather)?;
    let mut zones: Vec<ZoneRuntime> = model
        .zones
        .iter()
        .map(|zone| {
            let t0 = zone
                .heating_schedule
                .map(|i| model.schedules[i].value(start.0, start.1, start.2, 0))
                .unwrap_or(20.0);
            ZoneRuntime {
                temperature: t0,
                humidity: 0.008,
                history: [t0, t0, t0],
                humidity_history: [0.008, 0.008, 0.008],
            }
        })
        .collect();
    let timestep_hours = 1.0 / model.timesteps_per_hour as f64;
    let mut surfaces = Vec::new();
    for surface in &model.surfaces {
        let t0 = zones[surface.zone].temperature;
        surfaces.push(build_surface(surface, timestep_hours, t0)?);
    }
    let mut psych = PsychState::default();

    for _ in 0..options.warmup_days {
        step_day(
            model,
            weather,
            &location,
            start,
            &mut zones,
            &mut surfaces,
            &mut psych,
            None,
        )?;
    }

    let mut steps = Vec::new();
    let mut heating_j = 0.0;
    let mut cooling_j = 0.0;
    for date in enumerate_dates(start, end) {
        let mut day_steps = Vec::new();
        step_day(
            model,
            weather,
            &location,
            date,
            &mut zones,
            &mut surfaces,
            &mut psych,
            Some(&mut day_steps),
        )?;
        for step in &day_steps {
            let dt = 3600.0 / model.timesteps_per_hour as f64;
            for zone in &step.zones {
                heating_j += zone.heating_w * dt;
                cooling_j += zone.cooling_w * dt;
            }
        }
        steps.extend(day_steps);
    }
    Ok(Simulation {
        steps,
        heating_j,
        cooling_j,
    })
}
