//! Rust port of the EnergyPlus whole-building energy simulation engine.
//!
//! Numerics in `psychrometrics`, `solar`, and the zone air predictor are
//! translated from NatLabRockies/EnergyPlus `develop` at `ae0eaba`. The
//! program is not an official EnergyPlus release. See `LICENSE.txt`.

#![forbid(unsafe_code)]

pub mod baseboard;
pub mod convection;
pub mod ctf;
pub mod epjson;
pub mod epw;
pub mod idf;
pub mod model;
pub mod psychrometrics;
pub mod radiosity;
pub mod report;
pub mod simple_glazing;
pub mod simulate;
pub mod solar;

pub const KELVIN: f64 = 273.15;
pub const TRIPLE_POINT_OF_WATER_K: f64 = 273.16;
pub const DEG_TO_RAD: f64 = std::f64::consts::PI / 180.0;
pub const SIGMA: f64 = 5.670374419e-8;
