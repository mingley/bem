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
use crate::model::{Boundary, Gain, Model, Surface, SurfaceKind, ZoneEquipment};
use crate::psychrometrics::{
    psy_cp_air_fn_w, psy_hg_air_fn_w_tdb, psy_rho_air_fn_pb_tdb_w, psy_w_fn_tdb_rh_pb, PsychState,
};
use crate::radiosity::{class_from_tilt, net_longwave, tilt_from_normal, EnclosureSurface};
use crate::simple_glazing::{absorb_transmitted, SolarReceiver};
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

#[allow(clippy::too_many_arguments)]
fn step_day(
    model: &Model,
    weather: &WeatherFile,
    location: &crate::model::Location,
    date: (i32, u32, u32),
    zones: &mut [ZoneRuntime],
    surfaces: &mut [SurfaceRuntime],
    psych: &mut PsychState,
    mut record: Option<&mut Vec<TimestepResult>>,
) -> Result<(), String> {
    let (year, month, day) = date;
    let solar = daily_solar_coeffs(day_of_year(year, month, day), location.latitude);
    let n = model.timesteps_per_hour;
    let dt = 3600.0 / n as f64;
    for hour in 1..=24 {
        let wx = weather.hour(year, month, day, hour).ok_or_else(|| {
            format!("weather file has no hour {year}-{month:02}-{day:02} hour {hour}")
        })?;
        for ts in 1..=n {
            let end_minute = (hour - 1) * 60 + (ts as u32 * 60) / n as u32;
            let clock = (hour - 1) as f64 + (ts as f64 - 0.5) / n as f64;
            let sun = sun4(
                clock,
                &solar,
                location.latitude,
                location.longitude,
                location.time_zone,
            );
            let zone_steps =
                step_zone_timestep(model, wx, &sun, end_minute, dt, zones, surfaces, psych);
            if let Some(out) = record.as_mut() {
                out.push(TimestepResult {
                    year,
                    month,
                    day,
                    hour,
                    end_minute,
                    outdoor_db: wx.dry_bulb,
                    zones: zone_steps,
                });
            }
        }
    }
    Ok(())
}

fn step_zone_timestep(
    model: &Model,
    wx: &WeatherHour,
    sun: &SunVector,
    end_minute: u32,
    dt: f64,
    zones: &mut [ZoneRuntime],
    surfaces: &mut [SurfaceRuntime],
    psych: &mut PsychState,
) -> Vec<ZoneStep> {
    let pressure = if wx.pressure > 1000.0 {
        wx.pressure
    } else {
        101325.0
    };
    let rh = if wx.relative_humidity > 1.5 {
        wx.relative_humidity / 100.0
    } else {
        wx.relative_humidity
    };
    let w_out = psy_w_fn_tdb_rh_pb(wx.dry_bulb, rh.clamp(0.0, 1.0), pressure);
    let sky = sky_temp(wx.horizontal_ir, wx.dry_bulb);
    let (surface_solar, zone_solar) = interior_solar(model, sun, wx);
    let humidity0: Vec<f64> = zones.iter().map(|z| z.humidity).collect();
    let update_conduction: Vec<bool> = surfaces
        .iter()
        .map(|surface| surface.ctf.is_some() && surface.steps_until_update == 0)
        .collect();
    let mut interior: Vec<f64> = surfaces.iter().map(|surface| surface.t_in).collect();
    let mut h_interior = vec![0.1; surfaces.len()];
    let mut loads: Vec<Delivered> = model
        .zones
        .iter()
        .map(|zone| Delivered {
            q: 0.0,
            ideal_heating_w: 0.0,
            baseboards: baseboard_steps(zone),
        })
        .collect();

    for _ in 0..3 {
        for (zi, w) in humidity0.iter().enumerate() {
            zones[zi].humidity = *w;
        }
        let longwave_flux = enclosure_longwave(model, &interior);
        for (si, surface) in model.surfaces.iter().enumerate() {
            if !surface.active {
                continue;
            }
            let t_zone = zones[surface.zone].temperature;
            let (t_air, solar, use_sky) = boundary_driver(surface, zones, wx, sun);
            let t_so = surfaces[si].t_out;
            let t_si = surfaces[si].t_in;
            let h_ext = exterior_film(surface, t_so, t_air, wx);
            h_interior[si] = tarp_interior_hc(t_si, t_zone, surface.normal[2]);
            let q_rad = radiant_share(model, surface, end_minute, wx);
            let q_lw = longwave_flux[si] * surface.area;
            if update_conduction[si] {
                if let Some(t_int) = commit_conduction(
                    surface,
                    &mut surfaces[si],
                    t_air,
                    sky,
                    h_ext,
                    h_interior[si],
                    solar,
                    use_sky,
                    t_zone,
                    q_rad + q_lw + surface_solar[si],
                ) {
                    interior[si] = t_int;
                }
            }
        }
        for (zi, zone) in model.zones.iter().enumerate() {
            let update = zone_balance(
                model,
                zone,
                zi,
                &zones[zi],
                &model.surfaces,
                &interior,
                &h_interior,
                wx,
                w_out,
                pressure,
                end_minute,
                dt,
                psych,
                zone_solar[zi],
            );
            loads[zi] = Delivered {
                q: update.q,
                ideal_heating_w: update.ideal_heating_w,
                baseboards: update.baseboards,
            };
            zones[zi].temperature = update.temperature;
            zones[zi].humidity = update.humidity;
        }
    }

    for (runtime, update) in surfaces.iter_mut().zip(&update_conduction) {
        if *update {
            shift_conduction_history(runtime);
            let interval = runtime
                .ctf
                .as_ref()
                .map(|ctf| ctf.zone_steps_per_update)
                .unwrap_or(1);
            runtime.steps_until_update = interval;
        }
        if runtime.ctf.is_some() && runtime.steps_until_update > 0 {
            runtime.steps_until_update -= 1;
        }
    }

    for zone in zones.iter_mut() {
        zone.history[2] = zone.history[1];
        zone.history[1] = zone.history[0];
        zone.history[0] = zone.temperature;
        zone.humidity_history[2] = zone.humidity_history[1];
        zone.humidity_history[1] = zone.humidity_history[0];
        zone.humidity_history[0] = zone.humidity;
    }

    model
        .zones
        .iter()
        .enumerate()
        .map(|(zi, zone)| ZoneStep {
            name: zone.name.clone(),
            temperature: zones[zi].temperature,
            humidity_ratio: zones[zi].humidity,
            heating_w: loads[zi].q.max(0.0),
            cooling_w: (-loads[zi].q).max(0.0),
            ideal_heating_w: loads[zi].ideal_heating_w,
            baseboards: loads[zi].baseboards.clone(),
        })
        .collect()
}

fn zone_balance(
    model: &Model,
    zone: &crate::model::Zone,
    zi: usize,
    runtime: &ZoneRuntime,
    surfaces: &[Surface],
    interior: &[f64],
    h_interior: &[f64],
    wx: &WeatherHour,
    w_out: f64,
    pressure: f64,
    end_minute: u32,
    dt: f64,
    psych: &mut PsychState,
    zone_solar_air_w: f64,
) -> AirUpdate {
    let (year, month, day) = (wx.year, wx.month, wx.day);
    let minutes = end_minute;
    let mut conv = zone_solar_air_w;
    let mut radiant = 0.0;
    let mut latent_w = 0.0;
    for gain in &zone.gains {
        match gain {
            Gain::Power {
                schedule,
                convective_w,
                radiant_w,
                latent_w: lat,
            } => {
                let f = model.schedules[*schedule].value(year, month, day, minutes);
                conv += convective_w * f;
                radiant += radiant_w * f;
                latent_w += lat * f;
            }
            Gain::People {
                occupancy,
                activity,
                count,
                fraction_radiant,
                sensible_fraction,
            } => {
                let n = count * model.schedules[*occupancy].value(year, month, day, minutes);
                let act = activity
                    .map(|i| model.schedules[i].value(year, month, day, minutes))
                    .unwrap_or(120.0);
                let total = n * act;
                let sensible = total * *sensible_fraction;
                latent_w += total - sensible;
                radiant += sensible * *fraction_radiant;
                conv += sensible * (1.0 - *fraction_radiant);
            }
        }
    }
    let _ = radiant;

    let rho = psy_rho_air_fn_pb_tdb_w(pressure, runtime.temperature, runtime.humidity);
    let cp = psy_cp_air_fn_w(runtime.humidity);
    let mut sum_ha = 0.0;
    let mut sum_ha_tsurf = 0.0;
    let mut sum_mcp = 0.0;
    let mut sum_mcp_t = 0.0;
    for inf in &zone.infiltration_m3_s {
        let frac = inf
            .schedule
            .map(|i| model.schedules[i].value(year, month, day, minutes))
            .unwrap_or(1.0);
        let m = rho * inf.flow_m3_s * frac;
        sum_mcp += m * cp;
        sum_mcp_t += m * cp * wx.dry_bulb;
    }
    for (si, surface) in surfaces.iter().enumerate() {
        if surface.zone != zi || !surface.active {
            continue;
        }
        match &surface.kind {
            SurfaceKind::Opaque { .. } => {
                let h_in = h_interior[si];
                sum_ha += h_in * surface.area;
                sum_ha_tsurf += h_in * surface.area * interior[si];
            }
            SurfaceKind::Window {
                u_factor,
                inward_absorptance,
                ..
            } => {
                if matches!(surface.boundary, Boundary::Outdoors) {
                    sum_ha += u_factor * surface.area;
                    sum_ha_tsurf += u_factor * surface.area * wx.dry_bulb;
                    let _ = inward_absorptance;
                }
            }
        }
    }

    let air_cap = zone.volume * rho * cp / dt;
    let hist =
        3.0 * runtime.history[0] - 1.5 * runtime.history[1] + (1.0 / 3.0) * runtime.history[2];
    let temp_dep = sum_ha + sum_mcp;
    let temp_ind = conv + sum_ha_tsurf + sum_mcp_t;
    let dep_load = (11.0 / 6.0) * air_cap + temp_dep;
    let ind_load = air_cap * hist + temp_ind;
    let delivered = zone_system(model, zone, year, month, day, minutes, dep_load, ind_load);
    let q = delivered.q;
    let t_new = if dep_load.abs() < 1.0e-8 {
        runtime.temperature
    } else {
        (temp_ind + q + air_cap * hist) / dep_load
    };

    // Third-order humidity corrector from ZoneTempPredictorCorrector::correctHumRat.
    // Latent watts are converted with the vapor enthalpy PsyHgAirFnWTdb, matching EnergyPlus.
    let hg = psy_hg_air_fn_w_tdb(runtime.humidity, t_new).max(1.0e5);
    let m_lat = latent_w / hg;
    let mut m_inf = 0.0;
    for inf in &zone.infiltration_m3_s {
        let frac = inf
            .schedule
            .map(|i| model.schedules[i].value(year, month, day, minutes))
            .unwrap_or(1.0);
        m_inf += rho * inf.flow_m3_s * frac;
    }
    let rho_moist = psy_rho_air_fn_pb_tdb_w(pressure, t_new, runtime.humidity);
    let c_moist = (rho_moist * zone.volume / dt).max(1.0e-8);
    let a_flow = m_inf;
    let b_moist = m_lat + m_inf * w_out;
    let w_hist = 3.0 * runtime.humidity_history[0] - 1.5 * runtime.humidity_history[1]
        + (1.0 / 3.0) * runtime.humidity_history[2];
    let mut w_new = (b_moist + c_moist * w_hist) / ((11.0 / 6.0) * c_moist + a_flow);
    if w_new < 0.0 {
        w_new = 0.0;
    }
    let w_sat = psy_w_fn_tdb_rh_pb(t_new, 1.0, pressure);
    if w_new > w_sat {
        w_new = w_sat;
    }
    let _ = psych;
    AirUpdate {
        q,
        temperature: t_new,
        humidity: w_new,
        ideal_heating_w: delivered.ideal_heating_w,
        baseboards: delivered.baseboards,
    }
}

struct AirUpdate {
    q: f64,
    temperature: f64,
    humidity: f64,
    ideal_heating_w: f64,
    baseboards: Vec<BaseboardStep>,
}

struct Delivered {
    q: f64,
    ideal_heating_w: f64,
    baseboards: Vec<BaseboardStep>,
}

enum HeatKind {
    Ideal {
        capacity: f64,
    },
    Baseboard {
        index: usize,
        capacity: f64,
        efficiency: f64,
        available: bool,
    },
}

struct HeatSlot {
    sequence: i32,
    order: usize,
    kind: HeatKind,
}

fn baseboard_steps(zone: &crate::model::Zone) -> Vec<BaseboardStep> {
    zone.equipment
        .iter()
        .filter_map(|item| match item {
            ZoneEquipment::ElectricBaseboard { name, .. } => Some(BaseboardStep {
                name: name.clone(),
                heating_w: 0.0,
                electric_w: 0.0,
            }),
            ZoneEquipment::IdealLoads { .. } => None,
        })
        .collect()
}

fn zone_system(
    model: &Model,
    zone: &crate::model::Zone,
    year: i32,
    month: u32,
    day: u32,
    minutes: u32,
    dep_load: f64,
    ind_load: f64,
) -> Delivered {
    let baseboards = baseboard_steps(zone);
    let (Some(heat_i), Some(cool_i)) = (zone.heating_schedule, zone.cooling_schedule) else {
        return Delivered {
            q: 0.0,
            ideal_heating_w: 0.0,
            baseboards,
        };
    };
    let t_heat = model.schedules[heat_i].value(year, month, day, minutes);
    let t_cool = model.schedules[cool_i].value(year, month, day, minutes);
    let load_heat = dep_load * t_heat - ind_load;
    let load_cool = dep_load * t_cool - ind_load;
    if load_heat > 0.0 {
        heat_sequence(
            model, zone, year, month, day, minutes, load_heat, baseboards,
        )
    } else if load_cool < 0.0 {
        let q = cool_sequence(zone, load_cool);
        Delivered {
            q,
            ideal_heating_w: 0.0,
            baseboards,
        }
    } else {
        Delivered {
            q: 0.0,
            ideal_heating_w: 0.0,
            baseboards,
        }
    }
}

fn heat_sequence(
    model: &Model,
    zone: &crate::model::Zone,
    year: i32,
    month: u32,
    day: u32,
    minutes: u32,
    load_heat: f64,
    mut baseboards: Vec<BaseboardStep>,
) -> Delivered {
    let mut slots = heating_slots(model, zone, year, month, day, minutes);
    slots.sort_by(|left, right| {
        left.sequence
            .cmp(&right.sequence)
            .then(left.order.cmp(&right.order))
    });
    let mut remaining = load_heat;
    let mut ideal = 0.0;
    for slot in slots {
        match slot.kind {
            HeatKind::Ideal { capacity } => {
                let q = remaining.max(0.0).min(capacity);
                ideal += q;
                remaining -= q;
            }
            HeatKind::Baseboard {
                index,
                capacity,
                efficiency,
                available,
            } => {
                // SimElectricBaseboard passes RemainingOutputReqToHeatSP.
                let (q, elec) = convective_electric(remaining, capacity, efficiency, available);
                baseboards[index].heating_w = q;
                baseboards[index].electric_w = elec;
                remaining -= q;
            }
        }
    }
    let baseboard_heat: f64 = baseboards.iter().map(|step| step.heating_w).sum();
    Delivered {
        q: ideal + baseboard_heat,
        ideal_heating_w: ideal,
        baseboards,
    }
}

fn cool_sequence(zone: &crate::model::Zone, load_cool: f64) -> f64 {
    let mut slots: Vec<(i32, usize, f64)> = Vec::new();
    if zone.equipment.is_empty() {
        slots.push((1, 0, zone.cooling_capacity_w));
    } else {
        for (order, item) in zone.equipment.iter().enumerate() {
            if let ZoneEquipment::IdealLoads {
                cooling_capacity_w,
                cooling_sequence,
                ..
            } = item
            {
                if *cooling_sequence > 0 {
                    slots.push((*cooling_sequence, order, *cooling_capacity_w));
                }
            }
        }
    }
    slots.sort_by(|left, right| left.0.cmp(&right.0).then(left.1.cmp(&right.1)));
    let mut remaining = load_cool;
    let mut cool = 0.0;
    for (_, _, capacity) in slots {
        let q = if remaining < 0.0 {
            remaining.max(-capacity)
        } else {
            0.0
        };
        cool += q;
        remaining -= q;
    }
    cool
}

fn heating_slots(
    model: &Model,
    zone: &crate::model::Zone,
    year: i32,
    month: u32,
    day: u32,
    minutes: u32,
) -> Vec<HeatSlot> {
    if zone.equipment.is_empty() {
        return vec![HeatSlot {
            sequence: 1,
            order: 0,
            kind: HeatKind::Ideal {
                capacity: zone.heating_capacity_w,
            },
        }];
    }
    let mut slots = Vec::new();
    let mut baseboard_index = 0;
    for (order, item) in zone.equipment.iter().enumerate() {
        match item {
            ZoneEquipment::IdealLoads {
                heating_capacity_w,
                heating_sequence,
                ..
            } => {
                if *heating_sequence > 0 {
                    slots.push(HeatSlot {
                        sequence: *heating_sequence,
                        order,
                        kind: HeatKind::Ideal {
                            capacity: *heating_capacity_w,
                        },
                    });
                }
            }
            ZoneEquipment::ElectricBaseboard {
                capacity_w,
                efficiency,
                availability,
                heating_sequence,
                ..
            } => {
                if *heating_sequence > 0 {
                    let available = match availability {
                        None => true,
                        Some(index) => {
                            model.schedules[*index].value(year, month, day, minutes) > 0.0
                        }
                    };
                    slots.push(HeatSlot {
                        sequence: *heating_sequence,
                        order,
                        kind: HeatKind::Baseboard {
                            index: baseboard_index,
                            capacity: *capacity_w,
                            efficiency: *efficiency,
                            available,
                        },
                    });
                }
                baseboard_index += 1;
            }
        }
    }
    slots
}

fn radiant_share(model: &Model, surface: &Surface, end_minute: u32, wx: &WeatherHour) -> f64 {
    if !matches!(surface.kind, SurfaceKind::Opaque { .. }) || !surface.active {
        return 0.0;
    }
    let zone = &model.zones[surface.zone];
    let mut radiant = 0.0;
    for gain in &zone.gains {
        match gain {
            Gain::Power {
                schedule,
                radiant_w,
                ..
            } => {
                radiant += radiant_w
                    * model.schedules[*schedule].value(wx.year, wx.month, wx.day, end_minute);
            }
            Gain::People {
                occupancy,
                activity,
                count,
                fraction_radiant,
                sensible_fraction,
            } => {
                let n = count
                    * model.schedules[*occupancy].value(wx.year, wx.month, wx.day, end_minute);
                let act = activity
                    .map(|i| model.schedules[i].value(wx.year, wx.month, wx.day, end_minute))
                    .unwrap_or(120.0);
                radiant += n * act * sensible_fraction * fraction_radiant;
            }
        }
    }
    let weight: f64 = model
        .surfaces
        .iter()
        .filter(|s| {
            s.zone == surface.zone && s.active && matches!(s.kind, SurfaceKind::Opaque { .. })
        })
        .map(|s| s.area * s.thermal_absorptance_int.max(0.05))
        .sum();
    if weight <= 0.0 {
        0.0
    } else {
        radiant * surface.area * surface.thermal_absorptance_int.max(0.05) / weight
    }
}

fn boundary_driver(
    surface: &Surface,
    zones: &[ZoneRuntime],
    wx: &WeatherHour,
    sun: &SunVector,
) -> (f64, f64, bool) {
    match surface.boundary {
        Boundary::Adiabatic => (zones[surface.zone].temperature, 0.0, false),
        Boundary::Ground => (18.0, 0.0, false),
        Boundary::Zone(other) => (zones[other].temperature, 0.0, false),
        Boundary::Outdoors => {
            let incident =
                incident_solar(sun, surface.normal, wx.direct_normal, wx.diffuse_horizontal);
            (wx.dry_bulb, incident, true)
        }
    }
}

fn exterior_film(surface: &Surface, t_surf: f64, t_air: f64, wx: &WeatherHour) -> f64 {
    match surface.boundary {
        Boundary::Adiabatic => 0.0,
        Boundary::Outdoors => doe2_exterior_hc(
            t_surf,
            t_air,
            surface.normal[2],
            wx.wind_speed,
            wx.wind_direction,
            surface.azimuth_deg,
            surface.roughness,
        ),
        Boundary::Ground | Boundary::Zone(_) => {
            // No wind on a buried or interzone face. Natural convection uses the outward cosine.
            crate::convection::ashrae_tarp_natural(t_surf, t_air, surface.normal[2])
                .max(crate::convection::LOW_H_CONV)
        }
    }
}

fn interior_solar(model: &Model, sun: &SunVector, wx: &WeatherHour) -> (Vec<f64>, Vec<f64>) {
    let mut surface_watts = vec![0.0; model.surfaces.len()];
    let mut zone_watts = vec![0.0; model.zones.len()];
    for zone in 0..model.zones.len() {
        let mut beam_w = 0.0;
        let mut diffuse_w = 0.0;
        let mut inward_w = 0.0;
        for surface in model.surfaces.iter().filter(|s| s.zone == zone && s.active) {
            let SurfaceKind::Window {
                transmittance,
                inward_absorptance,
                ..
            } = surface.kind
            else {
                continue;
            };
            if !matches!(surface.boundary, Boundary::Outdoors) {
                continue;
            }
            let (beam, diffuse) =
                incident_parts(sun, surface.normal, wx.direct_normal, wx.diffuse_horizontal);
            beam_w += transmittance * beam * surface.area;
            diffuse_w += transmittance * diffuse * surface.area;
            inward_w += inward_absorptance * (beam + diffuse) * surface.area;
        }
        let mut receivers = Vec::new();
        let mut index = Vec::new();
        for (si, surface) in model.surfaces.iter().enumerate() {
            if surface.zone != zone || !surface.active {
                continue;
            }
            if !matches!(surface.kind, SurfaceKind::Opaque { .. }) {
                continue;
            }
            index.push(si);
            receivers.push(SolarReceiver {
                area: surface.area,
                inside_absorptance: surface.solar_absorptance_int,
                floor: surface.normal[2] < -0.5,
            });
        }
        let absorbed = absorb_transmitted(&receivers, beam_w, diffuse_w);
        let absorbed_sum: f64 = absorbed.iter().sum();
        for (local, si) in index.into_iter().enumerate() {
            surface_watts[si] += absorbed[local];
        }
        zone_watts[zone] += inward_w + (beam_w + diffuse_w - absorbed_sum).max(0.0);
    }
    (surface_watts, zone_watts)
}

fn incident_parts(sun: &SunVector, normal: [f64; 3], dni: f64, dhi: f64) -> (f64, f64) {
    const ALBEDO: f64 = 0.2;
    let cos_tilt = normal[2];
    let f_sky = ((1.0 + cos_tilt) * 0.5).clamp(0.0, 1.0);
    let f_ground = ((1.0 - cos_tilt) * 0.5).clamp(0.0, 1.0);
    let (beam, horizontal) = if sun.up {
        let cos_inc = (normal[0] * sun.x + normal[1] * sun.y + normal[2] * sun.z).max(0.0);
        (dni * cos_inc, dni * sun.z.max(0.0))
    } else {
        (0.0, 0.0)
    };
    (beam, dhi * f_sky + (horizontal + dhi) * ALBEDO * f_ground)
}

fn incident_solar(sun: &SunVector, normal: [f64; 3], dni: f64, dhi: f64) -> f64 {
    let (beam, diffuse) = incident_parts(sun, normal, dni, dhi);
    beam + diffuse
}

fn sky_temp(horizontal_ir: f64, dry_bulb: f64) -> f64 {
    if horizontal_ir > 1.0 {
        (horizontal_ir / SIGMA).powf(0.25) - KELVIN
    } else {
        dry_bulb
    }
}

fn enclosure_longwave(model: &Model, interior: &[f64]) -> Vec<f64> {
    let mut flux = vec![0.0; model.surfaces.len()];
    for zone in 0..model.zones.len() {
        let mut index = Vec::new();
        let mut specs = Vec::new();
        let mut temps = Vec::new();
        for (si, surface) in model.surfaces.iter().enumerate() {
            if surface.zone != zone || !surface.active {
                continue;
            }
            if !matches!(surface.kind, SurfaceKind::Opaque { .. }) {
                continue;
            }
            let tilt = tilt_from_normal(surface.normal[2]);
            index.push(si);
            specs.push(EnclosureSurface {
                area: surface.area,
                azimuth_deg: surface.azimuth_deg,
                tilt_deg: tilt,
                emissivity: surface.thermal_absorptance_int,
                class: class_from_tilt(tilt),
            });
            temps.push(interior[si]);
        }
        let net = net_longwave(&specs, &temps);
        for (local, si) in index.into_iter().enumerate() {
            flux[si] = net[local];
        }
    }
    flux
}

fn build_surface(
    surface: &Surface,
    timestep_hours: f64,
    t0: f64,
) -> Result<SurfaceRuntime, String> {
    let SurfaceKind::Opaque { layers } = &surface.kind else {
        return Ok(SurfaceRuntime {
            ctf: None,
            t_out: t0,
            t_in: t0,
            t_out_hist: Vec::new(),
            t_in_hist: Vec::new(),
            q_out_hist: Vec::new(),
            q_in_hist: Vec::new(),
            steps_until_update: 0,
        });
    };
    let ctf_layers: Vec<CtfLayer> = layers
        .iter()
        .map(|layer| CtfLayer {
            thickness: layer.dx,
            conductivity: layer.k,
            density: layer.rho,
            specific_heat: layer.cp,
        })
        .collect();
    let ctf = conduction_transfer_functions(&ctf_layers, timestep_hours)
        .map_err(|err| format!("surface '{}': {err}", surface.name))?;
    let n = ctf.num_terms + 1;
    Ok(SurfaceRuntime {
        ctf: Some(ctf),
        t_out: t0,
        t_in: t0,
        t_out_hist: vec![t0; n],
        t_in_hist: vec![t0; n],
        q_out_hist: vec![0.0; n],
        q_in_hist: vec![0.0; n],
        steps_until_update: 0,
    })
}

fn commit_conduction(
    surface: &Surface,
    runtime: &mut SurfaceRuntime,
    t_db: f64,
    sky: f64,
    h_conv_ext: f64,
    h_int: f64,
    incident: f64,
    use_sky: bool,
    t_zone: f64,
    q_to_inside: f64,
) -> Option<f64> {
    let ctf = runtime.ctf.as_ref()?;
    if surface.area <= 0.0 {
        return None;
    }
    let (const_out, const_in) = history_conduction(ctf, runtime);
    let q_inside = q_to_inside / surface.area;
    let h_int = h_int.max(0.0);
    let mut h_film = h_conv_ext.max(0.0);
    let mut env = h_film * t_db;
    if use_sky {
        env += surface.solar_absorptance * incident;
        if surface.thermal_absorptance_ext > 0.0 && h_conv_ext > 0.0 {
            let tsk = sky + KELVIN;
            let tso = runtime.t_out + KELVIN;
            let h_rad =
                surface.thermal_absorptance_ext * SIGMA * (tsk * tsk + tso * tso) * (tsk + tso);
            h_film += h_rad;
            env += h_rad * sky;
        }
    }
    let y0 = ctf.cross[0];
    let x0 = ctf.outside[0];
    let z0 = ctf.inside[0];
    let t_old = runtime.t_in;
    let (t_out, t_in) = if matches!(surface.boundary, Boundary::Adiabatic) {
        let denom = z0 - y0 + h_int + SURFACE_DAMP;
        let t_in = (const_in + q_inside + h_int * t_zone + SURFACE_DAMP * t_old) / denom;
        (t_in, t_in)
    } else if y0 > 0.01 {
        // Quick conduction folds the inside balance into the outside face.
        let f1 = y0 / (z0 + h_int);
        let inside_drive = const_in + q_inside + h_int * t_zone;
        let denom = x0 + h_film - f1 * y0;
        let t_out = (-const_out + env + f1 * inside_drive) / denom;
        let t_in = (inside_drive + SURFACE_DAMP * t_old + y0 * t_out) / (z0 + h_int + SURFACE_DAMP);
        (t_out, t_in)
    } else {
        let t_out = (-const_out + env + y0 * runtime.t_in) / (x0 + h_film);
        let t_in = (const_in + q_inside + h_int * t_zone + SURFACE_DAMP * t_old + y0 * t_out)
            / (z0 + h_int + SURFACE_DAMP);
        (t_out, t_in)
    };
    runtime.t_out = t_out;
    runtime.t_in = t_in;
    Some(t_in)
}

fn history_conduction(ctf: &Ctf, runtime: &SurfaceRuntime) -> (f64, f64) {
    let mut outside = 0.0;
    let mut inside = 0.0;
    for term in 1..=ctf.num_terms {
        outside += ctf.outside[term] * runtime.t_out_hist[term]
            - ctf.cross[term] * runtime.t_in_hist[term]
            + ctf.flux[term] * runtime.q_out_hist[term];
        inside += ctf.cross[term] * runtime.t_out_hist[term]
            - ctf.inside[term] * runtime.t_in_hist[term]
            + ctf.flux[term] * runtime.q_in_hist[term];
    }
    (outside, inside)
}

fn shift_conduction_history(runtime: &mut SurfaceRuntime) {
    let Some(ctf) = runtime.ctf.clone() else {
        return;
    };
    let (const_out, const_in) = history_conduction(&ctf, runtime);
    let q_in = runtime.t_out * ctf.cross[0] - runtime.t_in * ctf.inside[0] + const_in;
    let q_out = runtime.t_out * ctf.outside[0] - runtime.t_in * ctf.cross[0] + const_out;
    for term in (1..runtime.t_in_hist.len()).rev() {
        runtime.t_out_hist[term] = runtime.t_out_hist[term - 1];
        runtime.t_in_hist[term] = runtime.t_in_hist[term - 1];
        runtime.q_out_hist[term] = runtime.q_out_hist[term - 1];
        runtime.q_in_hist[term] = runtime.q_in_hist[term - 1];
    }
    if !runtime.t_in_hist.is_empty() {
        runtime.t_out_hist[0] = runtime.t_out;
        runtime.t_in_hist[0] = runtime.t_in;
        runtime.q_out_hist[0] = q_out;
        runtime.q_in_hist[0] = q_in;
    }
}

fn resolve_period(
    model: &Model,
    weather: &WeatherFile,
) -> Result<((i32, u32, u32), (i32, u32, u32)), String> {
    let first = weather.hours.first().unwrap();
    let last = weather.hours.last().unwrap();
    let Some(period) = &model.run_period else {
        return Ok((
            (first.year, first.month, first.day),
            (last.year, last.month, last.day),
        ));
    };
    let begin_year = period.begin_year.unwrap_or(first.year);
    let mut end_year = period.end_year.unwrap_or(begin_year);
    let start = (begin_year, period.begin_month, period.begin_day);
    let mut end = (end_year, period.end_month, period.end_day);
    if date_index(end) < date_index(start) && period.end_year.is_none() {
        end_year = begin_year + 1;
        end = (end_year, period.end_month, period.end_day);
    }
    if date_index(end) < date_index(start) {
        return Err("RunPeriod ends before it starts".to_string());
    }
    Ok((start, end))
}

fn date_index((y, m, d): (i32, u32, u32)) -> i32 {
    y * 372 + m as i32 * 31 + d as i32
}

fn enumerate_dates(start: (i32, u32, u32), end: (i32, u32, u32)) -> Vec<(i32, u32, u32)> {
    let mut out = Vec::new();
    let mut cursor = start;
    loop {
        out.push(cursor);
        if cursor == end {
            break;
        }
        cursor = next_date(cursor);
        if out.len() > 3660 {
            break;
        }
    }
    out
}

fn next_date((y, m, d): (i32, u32, u32)) -> (i32, u32, u32) {
    let leap = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
    let mdays = [
        0,
        31,
        if leap { 29 } else { 28 },
        31,
        30,
        31,
        30,
        31,
        31,
        30,
        31,
        30,
        31,
    ];
    if d < mdays[m as usize] {
        (y, m, d + 1)
    } else if m < 12 {
        (y, m + 1, 1)
    } else {
        (y + 1, 1, 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::epw::parse_epw;
    use crate::model::{
        Boundary, Gain, Infiltration, Layer, RunPeriod, Schedule, Surface, SurfaceKind, Zone,
    };

    fn constant_weather(tdb: f64, hours: usize) -> WeatherFile {
        let mut rows = Vec::new();
        for h in 0..hours {
            let day = (h / 24) as u32 + 1;
            let hour = (h % 24) as u32 + 1;
            rows.push(format!(
                "2021,1,{day},{hour},60,flags,{tdb:.2},-5,50,101325,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9999,99999,9,999999999,0,0.04,0,88,0.2,0,0"
            ));
        }
        let header = "\
LOCATION,Test,IL,USA,Synthetic,999,41.78,-87.75,-6,182
DESIGN CONDITIONS,0
TYPICAL/EXTREME PERIODS,0
GROUND TEMPERATURES,0
HOLIDAYS/DAYLIGHT SAVINGS,No,0,0,0
COMMENTS 1,test
COMMENTS 2,test
DATA PERIODS,1,1,Data,Friday,1/1,1/2
";
        parse_epw(&format!("{header}{}\n", rows.join("\n"))).unwrap()
    }

    fn bare_model(gain_w: f64, heat: f64, cool: f64) -> Model {
        Model {
            location: None,
            north_axis: 0.0,
            world_coordinates: true,
            timesteps_per_hour: 6,
            run_period: Some(RunPeriod {
                begin_month: 1,
                begin_day: 1,
                begin_year: Some(2021),
                end_month: 1,
                end_day: 1,
                end_year: Some(2021),
            }),
            schedules: vec![
                Schedule::Constant(1.0),
                Schedule::Constant(heat),
                Schedule::Constant(cool),
            ],
            zones: vec![Zone {
                name: "Zone".into(),
                volume: 100.0,
                floor_area: 30.0,
                heating_schedule: Some(1),
                cooling_schedule: Some(2),
                heating_capacity_w: f64::INFINITY,
                cooling_capacity_w: f64::INFINITY,
                equipment: vec![],
                gains: vec![Gain::Power {
                    schedule: 0,
                    convective_w: gain_w,
                    radiant_w: 0.0,
                    latent_w: 0.0,
                }],
                infiltration_m3_s: vec![],
            }],
            surfaces: vec![],
            outputs: vec![],
            warnings: vec![],
        }
    }

    #[test]
    fn ideal_loads_remove_a_steady_internal_gain() {
        let model = bare_model(1000.0, 20.0, 24.0);
        let weather = constant_weather(24.0, 24);
        let sim = simulate(&model, &weather, &SimOptions { warmup_days: 0 }).unwrap();
        let last = sim.steps.last().unwrap();
        assert!(
            (last.zones[0].temperature - 24.0).abs() < 0.02,
            "temp {}",
            last.zones[0].temperature
        );
        assert!(
            (last.zones[0].cooling_w - 1000.0).abs() < 5.0,
            "cooling {}",
            last.zones[0].cooling_w
        );
        assert!(last.zones[0].heating_w < 1.0);
    }

    #[test]
    fn baseboard_takes_its_sequence_share_and_ideal_loads_the_rest() {
        let mut model = bare_model(-2000.0, 20.0, 24.0);
        model.zones[0].equipment = vec![
            ZoneEquipment::ElectricBaseboard {
                name: "BB".into(),
                capacity_w: 500.0,
                efficiency: 0.5,
                availability: None,
                heating_sequence: 1,
            },
            ZoneEquipment::IdealLoads {
                heating_capacity_w: f64::INFINITY,
                cooling_capacity_w: f64::INFINITY,
                heating_sequence: 2,
                cooling_sequence: 1,
            },
        ];
        let weather = constant_weather(0.0, 24);
        let sim = simulate(&model, &weather, &SimOptions { warmup_days: 0 }).unwrap();
        let last = sim.steps.last().unwrap();
        let zone = &last.zones[0];
        assert!(
            (zone.temperature - 20.0).abs() < 0.02,
            "temp {}",
            zone.temperature
        );
        assert!(
            (zone.baseboards[0].heating_w - 500.0).abs() < 1.0,
            "baseboard {}",
            zone.baseboards[0].heating_w
        );
        assert!(
            (zone.baseboards[0].electric_w - 1000.0).abs() < 1.0,
            "electric {}",
            zone.baseboards[0].electric_w
        );
        assert!(
            (zone.ideal_heating_w - 1500.0).abs() < 5.0,
            "ideal {}",
            zone.ideal_heating_w
        );
        assert!(
            (zone.heating_w - 2000.0).abs() < 5.0,
            "total {}",
            zone.heating_w
        );

        model.zones[0].equipment[0] = ZoneEquipment::ElectricBaseboard {
            name: "BB".into(),
            capacity_w: 500.0,
            efficiency: 0.5,
            availability: None,
            heating_sequence: 2,
        };
        if let ZoneEquipment::IdealLoads {
            heating_sequence, ..
        } = &mut model.zones[0].equipment[1]
        {
            *heating_sequence = 1;
        }
        let sim = simulate(&model, &weather, &SimOptions { warmup_days: 0 }).unwrap();
        let zone = &sim.steps.last().unwrap().zones[0];
        assert!(
            zone.baseboards[0].heating_w < 1.0,
            "baseboard after unlimited ideal loads {}",
            zone.baseboards[0].heating_w
        );
        assert!(
            (zone.ideal_heating_w - 2000.0).abs() < 5.0,
            "ideal {}",
            zone.ideal_heating_w
        );
    }

    #[test]
    fn steady_wall_heat_loss_matches_overall_u() {
        let k = 0.1;
        let thickness = 0.1;
        let area = 1.0;
        let mut model = bare_model(0.0, 20.0, 50.0);
        model.surfaces.push(Surface {
            name: "Wall".into(),
            zone: 0,
            area,
            normal: [0.0, -1.0, 0.0],
            azimuth_deg: 180.0,
            roughness: crate::convection::Roughness::MediumRough,
            kind: SurfaceKind::Opaque {
                layers: vec![Layer {
                    k,
                    dx: thickness,
                    n: 1,
                    rho: 10.0,
                    cp: 100.0,
                }],
            },
            boundary: Boundary::Outdoors,
            solar_absorptance: 0.0,
            solar_absorptance_int: 0.0,
            thermal_absorptance_ext: 0.0,
            thermal_absorptance_int: 0.9,
            active: true,
        });
        model.run_period.as_mut().unwrap().end_day = 2;
        let weather = constant_weather(0.0, 48);
        let sim = simulate(&model, &weather, &SimOptions { warmup_days: 0 }).unwrap();
        let last_hour: Vec<_> = sim.steps.iter().rev().take(6).collect();
        let heats: Vec<f64> = last_hour
            .iter()
            .map(|step| step.zones[0].heating_w)
            .collect();
        let mean = heats.iter().sum::<f64>() / heats.len() as f64;
        assert!(
            (4.0..25.0).contains(&mean),
            "steady heating {mean} W is outside the TARP/DOE-2 range for this wall"
        );
        assert!(
            heats.iter().all(|h| (h - mean).abs() / mean < 0.02),
            "heating did not settle: {heats:?}"
        );
        assert!((last_hour[0].zones[0].temperature - 20.0).abs() < 0.05);
        let _ = (k, thickness, area);
    }

    #[test]
    fn third_order_humidity_tracks_outdoor_air_plus_latent_gain() {
        let mut model = bare_model(0.0, 20.0, 20.0);
        model.timesteps_per_hour = 6;
        model.run_period.as_mut().unwrap().end_day = 2;
        let flow = 50.0 / 3600.0;
        model.zones[0].volume = 50.0;
        model.zones[0].infiltration_m3_s = vec![Infiltration {
            schedule: None,
            flow_m3_s: flow,
        }];
        model.zones[0].gains = vec![Gain::Power {
            schedule: 0,
            convective_w: 0.0,
            radiant_w: 0.0,
            latent_w: 100.0,
        }];
        let weather = constant_weather(20.0, 48);
        let sim = simulate(&model, &weather, &SimOptions { warmup_days: 0 }).unwrap();
        let last = &sim.steps.last().unwrap().zones[0];
        let w_out = crate::psychrometrics::psy_w_fn_tdb_rh_pb(20.0, 0.5, 101325.0);
        let rho = crate::psychrometrics::psy_rho_air_fn_pb_tdb_w(
            101325.0,
            last.temperature,
            last.humidity_ratio,
        );
        let hg = crate::psychrometrics::psy_hg_air_fn_w_tdb(last.humidity_ratio, last.temperature);
        let expected = w_out + (100.0 / hg) / (rho * flow);
        assert!(
            (last.humidity_ratio - expected).abs() < 2.0e-4,
            "humidity {} expected {expected}",
            last.humidity_ratio
        );
        assert!(
            last.humidity_ratio
                < crate::psychrometrics::psy_w_fn_tdb_rh_pb(last.temperature, 1.0, 101325.0)
        );
    }
}
