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
