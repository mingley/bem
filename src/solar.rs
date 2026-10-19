//! Daily solar coefficients and sun direction from EnergyPlus.
//!
//! `daily_solar_coeffs` is `WeatherManager::CalculateDailySolarCoeffs`.
//! `sun4` is `SolarShading::SUN4`.

use crate::DEG_TO_RAD;

/// Coefficients for one day of year. Angles are dimensionless sines and cosines.
#[derive(Clone, Debug)]
pub struct DailySolar {
    pub ashrae_a: f64,
    pub ashrae_b: f64,
    pub ashrae_c: f64,
    pub ann_var_sol_const: f64,
    pub equation_of_time: f64,
    pub sin_solar_declination: f64,
    pub cos_solar_declination: f64,
}

/// Direction cosines of the sun. `x` is east, `y` is north, `z` is up.
#[derive(Clone, Debug)]
pub struct SunVector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub up: bool,
}

fn fourier9(coef: &[f64; 9], sin_x: f64, cos_x: f64) -> f64 {
    let sin2 = sin_x * cos_x * 2.0;
    let cos2_m_sin2 = cos_x * cos_x - sin_x * sin_x;
    coef[0]
        + coef[1] * sin_x
        + coef[2] * cos_x
        + coef[3] * sin2
        + coef[4] * cos2_m_sin2
        + coef[5] * (sin_x * cos2_m_sin2 + cos_x * sin2)
        + coef[6] * (cos_x * cos2_m_sin2 - sin_x * sin2)
        + coef[7] * (2.0 * sin2 * cos2_m_sin2)
        + coef[8] * (cos2_m_sin2 * cos2_m_sin2 - sin2 * sin2)
}

/// ASHRAE sky coefficients, equation of time, and solar declination for `day_of_year` (1-366).
pub fn daily_solar_coeffs(day_of_year: i32, latitude_deg: f64) -> DailySolar {
    const SINE_SOL_DECL: [f64; 9] = [
        0.00561800,
        0.0657911,
        -0.392779,
        0.00064440,
        -0.00618495,
        -0.00010101,
        -0.00007951,
        -0.00011691,
        0.00002096,
    ];
    const EQ_OF_TIME: [f64; 9] = [
        0.00021971,
        -0.122649,
        0.00762856,
        -0.156308,
        -0.0530028,
        -0.00388702,
        -0.00123978,
        -0.00270502,
        -0.00167992,
    ];
    const ASHRAE_A: [f64; 9] = [
        1161.6685, 1.1554, 77.3575, -0.5359, -3.7622, 0.9875, -3.3924, -1.7445, 1.1198,
    ];
    const ASHRAE_B: [f64; 9] = [
        0.171631,
        -0.00400448,
        -0.0344923,
        0.00000209,
        0.00325428,
        -0.00085429,
        0.00229562,
        0.0009034,
        -0.0011867,
    ];
    const ASHRAE_C: [f64; 9] = [
        0.0905151,
        -0.00322522,
        -0.0407966,
        0.000104164,
        0.00745899,
        -0.00086461,
        0.0013111,
        0.000808275,
        -0.00170515,
    ];

    let day_correction = std::f64::consts::PI * 2.0 / 366.0;
    let mut x = day_correction * f64::from(day_of_year);
    let sin_x = x.sin();
    let cos_x = x.cos();
    let sin_decl = fourier9(&SINE_SOL_DECL, sin_x, cos_x);
    let cos_decl = (1.0 - sin_decl * sin_decl).sqrt();
    let equation_of_time = fourier9(&EQ_OF_TIME, sin_x, cos_x);
    let ann_var = 1.000047 + 0.000352615 * sin_x + 0.0334454 * cos_x;
    let a = fourier9(&ASHRAE_A, sin_x, cos_x);
    if latitude_deg < 0.0 {
        x -= std::f64::consts::PI;
    }
    let sin_x = x.sin();
    let cos_x = x.cos();
    let b = fourier9(&ASHRAE_B, sin_x, cos_x);
    let c = fourier9(&ASHRAE_C, sin_x, cos_x);
    DailySolar {
        ashrae_a: a,
        ashrae_b: b,
        ashrae_c: c,
        ann_var_sol_const: ann_var,
        equation_of_time,
        sin_solar_declination: sin_decl,
        cos_solar_declination: cos_decl,
    }
}

/// Sun direction cosines at `current_time_hours` (0 at midnight, 13.5 at 13:30).
pub fn sun4(
    current_time_hours: f64,
    solar: &DailySolar,
    latitude_deg: f64,
    longitude_deg: f64,
    time_zone_hours: f64,
) -> SunVector {
    const SUN_IS_UP: f64 = 0.00001;
    let sin_lat = (latitude_deg * DEG_TO_RAD).sin();
    let cos_lat = (latitude_deg * DEG_TO_RAD).cos();
    let meridian = time_zone_hours * 15.0;
    let hr_angle =
        15.0 * (12.0 - (current_time_hours + solar.equation_of_time)) + (meridian - longitude_deg);
    let h = hr_angle * DEG_TO_RAD;
    let z = solar.sin_solar_declination * sin_lat + solar.cos_solar_declination * cos_lat * h.cos();
    if z < SUN_IS_UP {
        return SunVector {
            x: 0.0,
            y: 0.0,
            z,
            up: false,
        };
    }
    let y = solar.sin_solar_declination * cos_lat - solar.cos_solar_declination * sin_lat * h.cos();
    let x = solar.cos_solar_declination * h.sin();
    SunVector { x, y, z, up: true }
}

/// Day of year for a Gregorian date. January 1 is 1.
pub fn day_of_year(year: i32, month: u32, day: u32) -> i32 {
    let leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
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
    let mut n = day as i32;
    for m in 1..month {
        n += mdays[m as usize];
    }
    n
}

/// 0 = Sunday .. 6 = Saturday. Sakamoto's method.
pub fn weekday(year: i32, month: u32, day: u32) -> u32 {
    let t = [0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4];
    let mut y = year;
    if month < 3 {
        y -= 1;
    }
    let dow = (y + y / 4 - y / 100 + y / 400 + t[(month as usize) - 1] + day as i32).rem_euclid(7);
    dow as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn june_sun_is_high_at_solar_noon_in_chicago() {
        let day = day_of_year(2021, 6, 21);
        let solar = daily_solar_coeffs(day, 41.78);
        // Time zone -6, longitude -87.75, so solar noon is a bit before 12:00 standard time.
        let mut best = SunVector {
            x: 0.0,
            y: 0.0,
            z: -1.0,
            up: false,
        };
        for step in 0..24 * 4 {
            let hour = step as f64 / 4.0;
            let sun = sun4(hour, &solar, 41.78, -87.75, -6.0);
            if sun.z > best.z {
                best = sun;
            }
        }
        assert!(best.up);
        // Chicago June solar altitude is about 71 degrees. sin(71 deg) ~= 0.945.
        assert!(best.z > 0.90, "noon altitude cosine z={}", best.z);
    }

    #[test]
    fn january_night_has_sun_down() {
        let solar = daily_solar_coeffs(day_of_year(2021, 1, 1), 41.78);
        let sun = sun4(2.0, &solar, 41.78, -87.75, -6.0);
        assert!(!sun.up);
    }

    #[test]
    fn known_weekday() {
        // 1 Jan 2021 was a Friday.
        assert_eq!(weekday(2021, 1, 1), 5);
    }
}
