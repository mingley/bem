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
