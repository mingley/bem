//! Psychrometric routines translated from EnergyPlus `Psychrometrics.cc` / `.hh`.
//!
//! Numeric clamps match the C++ functions. The recurring error reporter
//! (`EP_psych_errors`) and the result cache are not ported; the raw equations
//! are, which is what the cache stores.

use crate::{KELVIN, TRIPLE_POINT_OF_WATER_K};

/// Iteration state carried by `PsyTwbFnTdbWPb` and `PsyTsatFnPb`.
#[derive(Clone, Debug)]
pub struct PsychState {
    pub iconv_tol: f64,
    pub last_patm: f64,
    pub last_tboil: f64,
    pub press_save: f64,
    pub tsat_save: f64,
}

impl Default for PsychState {
    fn default() -> Self {
        Self {
            iconv_tol: 0.0001,
            last_patm: -99999.0,
            last_tboil: -99999.0,
            press_save: -99999.0,
            tsat_save: -99999.0,
        }
    }
}

/// Linear correction from EnergyPlus `General::Iterate`.
///
/// Solves `Y(X) = 0`. Returns `(result_x, converged)`.
pub fn iterate(x0: f64, y0: f64, tol: f64, iter: i32, x1: &mut f64, y1: &mut f64) -> (f64, bool) {
    const SMALL: f64 = 1.0e-9;
    const PERTURB: f64 = 0.1;
    if iter != 1 && ((x0 - *x1).abs() < tol || y0 == 0.0) {
        return (x0, true);
    }
    let result_x = if iter == 1 {
        if x0.abs() > SMALL {
            x0 * (1.0 + PERTURB)
        } else {
            PERTURB
        }
    } else {
        let mut dy = y0 - *y1;
        if dy.abs() < SMALL {
            dy = SMALL;
        }
        (y0 * *x1 - *y1 * x0) / dy
    };
    *x1 = x0;
    *y1 = y0;
    (result_x, false)
}

fn hum_ratio(dw: f64) -> f64 {
    dw.max(1.0e-5)
}

/// Moist-air density `{kg/m3}` from pressure `{Pa}`, dry-bulb `{C}`, humidity ratio.
pub fn psy_rho_air_fn_pb_tdb_w(pb: f64, tdb: f64, dw: f64) -> f64 {
    pb / (287.0 * (tdb + KELVIN) * (1.0 + 1.6077687 * hum_ratio(dw)))
}

/// Latent heat of vaporization of water in moist air `{J/kg}`.
pub fn psy_hfg_air_fn_w_tdb(_w: f64, t: f64) -> f64 {
    let temperature = t.max(0.0);
    (2_500_940.0 + 1858.95 * temperature) - (4180.0 * temperature)
}

/// Enthalpy of water vapor `{J/kg}`.
pub fn psy_hg_air_fn_w_tdb(_w: f64, t: f64) -> f64 {
    2_500_940.0 + 1858.95 * t
}

/// Moist-air enthalpy `{J/kg}` from dry-bulb `{C}` and humidity ratio.
pub fn psy_h_fn_tdb_w(tdb: f64, dw: f64) -> f64 {
    1.00484e3 * tdb + hum_ratio(dw) * (2.50094e6 + 1.85895e3 * tdb)
}

/// Moist-air specific heat `{J/kg-C}`.
pub fn psy_cp_air_fn_w(dw: f64) -> f64 {
    1.00484e3 + hum_ratio(dw) * 1.85895e3
}

/// Dry-bulb `{C}` from enthalpy `{J/kg}` and humidity ratio.
pub fn psy_tdb_fn_hw(h: f64, dw: f64) -> f64 {
    let w = hum_ratio(dw);
    (h - 2.50094e6 * w) / (1.00484e3 + 1.85895e3 * w)
}

/// Vapor density `{kg/m3}` using the 0 C lower-bound saturation curve.
pub fn psy_rhov_fn_tdb_rh_lbnd0c(tdb: f64, rh: f64) -> f64 {
    rh / (461.52 * (tdb + KELVIN)) * (23.7093 - 4111.0 / ((tdb + KELVIN) - 35.45)).exp()
}

/// Vapor density `{kg/m3}` from dry-bulb, humidity ratio, and pressure.
pub fn psy_rhov_fn_tdb_w_pb(tdb: f64, dw: f64, pb: f64) -> f64 {
    let w = hum_ratio(dw);
    w * pb / (461.52 * (tdb + KELVIN) * (w + 0.62198))
}

/// Relative humidity `{0-1}` from dry-bulb and vapor density, 0 C curve.
pub fn psy_rh_fn_tdb_rhov_lbnd0c(tdb: f64, rhovapor: f64) -> f64 {
    let rh = if rhovapor > 0.0 {
        rhovapor * 461.52 * (tdb + KELVIN) * (-23.7093 + 4111.0 / ((tdb + KELVIN) - 35.45)).exp()
    } else {
        0.0
    };
    if !(0.0..=1.0).contains(&rh) {
        rh.clamp(0.01, 1.0)
    } else {
        rh
    }
}

/// Saturation pressure `{Pa}` from dry-bulb `{C}`, Hyland & Wexler (ASHRAE HOF).
pub fn psy_psat_fn_temp(t: f64) -> f64 {
    let tkel = t + KELVIN;
    if tkel < 173.15 {
        0.001405102123874164
    } else if tkel < TRIPLE_POINT_OF_WATER_K {
        const C1: f64 = -5674.5359;
        const C2: f64 = 6.3925247;
        const C3: f64 = -0.9677843e-2;
        const C4: f64 = 0.62215701e-6;
        const C5: f64 = 0.20747825e-8;
        const C6: f64 = -0.9484024e-12;
        const C7: f64 = 4.1635019;
        (C1 / tkel + C2 + tkel * (C3 + tkel * (C4 + tkel * (C5 + C6 * tkel))) + C7 * tkel.ln())
            .exp()
    } else if tkel <= 473.15 {
        const C8: f64 = -5800.2206;
        const C9: f64 = 1.3914993;
        const C10: f64 = -0.048640239;
        const C11: f64 = 0.41764768e-4;
        const C12: f64 = -0.14452093e-7;
        const C13: f64 = 6.5459673;
        (C8 / tkel + C9 + tkel * (C10 + tkel * (C11 + tkel * C12)) + C13 * tkel.ln()).exp()
    } else {
        1_555_073.745636215
    }
}

/// Specific volume `{m3/kg}`. Returns 0.83 when the ideal-gas result is negative.
pub fn psy_v_fn_tdb_w_pb(tdb: f64, dw: f64, pb: f64) -> f64 {
    let w = hum_ratio(dw);
    let v = 1.59473e2 * (1.0 + 1.6078 * w) * (1.8 * tdb + 492.0) / pb;
    if v < 0.0 {
        0.83
    } else {
        v
    }
}

/// Humidity ratio from dry-bulb and enthalpy. Floors at `1e-5`.
pub fn psy_w_fn_tdb_h(tdb: f64, h: f64) -> f64 {
    let w = (h - 1.00484e3 * tdb) / (2.50094e6 + 1.85895e3 * tdb);
    if w < 0.0 {
        1.0e-5
    } else {
        w
    }
}

/// Relative humidity `{0-1}` from dry-bulb, humidity ratio, and barometric pressure.
pub fn psy_rh_fn_tdb_w_pb(state: &mut PsychState, tdb: f64, dw: f64, pb: f64) -> f64 {
    let pws = psy_psat_fn_temp(tdb);
    let _ = state;
    let w = hum_ratio(dw);
    let u = w / (0.62198 * pws / (pb - pws));
    let rh = u / (1.0 - (1.0 - u) * (pws / pb));
    if !(0.0..=1.0).contains(&rh) {
        rh.clamp(0.01, 1.0)
    } else {
        rh
    }
}

/// Humidity ratio from dew-point `{C}` and barometric pressure `{Pa}`.
pub fn psy_w_fn_tdp_pb(tdp: f64, pb: f64) -> f64 {
    let pdew = psy_psat_fn_temp(tdp);
    let w = pdew * 0.62198 / (pb - pdew);
    if w < 0.0 {
        let mut delta_t = 0.0;
        let mut pdew1 = pdew;
        while pdew1 >= pb {
            delta_t += 1.0;
            pdew1 = psy_psat_fn_temp(tdp - delta_t);
        }
        pdew1 * 0.62198 / (pb - pdew1)
    } else {
        w
    }
}

/// Humidity ratio from dry-bulb, relative humidity `{0-1}`, and pressure.
pub fn psy_w_fn_tdb_rh_pb(tdb: f64, rh: f64, pb: f64) -> f64 {
    let pdew = rh * psy_psat_fn_temp(tdb);
    let w = pdew * 0.62198 / (pb - pdew).max(1000.0);
    if w < 1.0e-5 {
        1.0e-5
    } else {
        w
    }
}

/// Humidity ratio from dry-bulb, wet-bulb, and pressure.
pub fn psy_w_fn_tdb_twb_pb(tdb: f64, twb_in: f64, pb: f64) -> f64 {
    let twb = if twb_in > tdb { tdb } else { twb_in };
    let pwet = psy_psat_fn_temp(twb);
    let wet = 0.62198 * pwet / (pb - pwet);
    let w = ((2501.0 - 2.381 * twb) * wet - (tdb - twb)) / (2501.0 + 1.805 * tdb - 4.186 * twb);
    if w < 0.0 {
        psy_w_fn_tdb_rh_pb(tdb, 0.0001, pb)
    } else {
        w
    }
}

/// Enthalpy from dry-bulb, relative humidity, and pressure.
pub fn psy_h_fn_tdb_rh_pb(tdb: f64, rh: f64, pb: f64) -> f64 {
    psy_h_fn_tdb_w(tdb, psy_w_fn_tdb_rh_pb(tdb, rh, pb).max(1.0e-5))
}

/// Saturation temperature `{C}` at pressure `{Pa}`.
pub fn psy_tsat_fn_pb(state: &mut PsychState, press: f64) -> f64 {
    const ITMAX: i32 = 50;
    const CONV_TOL: f64 = 0.0001;
    if press == state.press_save {
        return state.tsat_save;
    }
    state.press_save = press;
    let t_sat = if press >= 1_555_000.0 {
        200.0
    } else if press <= 0.0017 {
        -100.0
    } else if (611.000..611.25).contains(&press) {
        0.0
    } else {
        let mut t_sat = 100.0;
        let mut x1 = 0.0;
        let mut y1 = 0.0;
        for iter in 1..=ITMAX {
            let p_sat = psy_psat_fn_temp(t_sat);
            let error = press - p_sat;
            let (result_x, converged) = iterate(t_sat, error, CONV_TOL, iter, &mut x1, &mut y1);
            t_sat = result_x;
            if converged {
                break;
            }
        }
        t_sat
    };
    state.tsat_save = t_sat;
    t_sat
}

/// Wet-bulb temperature `{C}` from dry-bulb, humidity ratio, and pressure.
pub fn psy_twb_fn_tdb_w_pb(state: &mut PsychState, tdb: f64, dw: f64, patm: f64) -> f64 {
    const ITMAX: i32 = 100;
    let mut w = dw;
    if w < 0.0 {
        w = 1.0e-5;
    }
    let t_boil = if patm != state.last_patm {
        let t_boil = psy_tsat_fn_pb(state, patm);
        state.last_patm = patm;
        state.last_tboil = t_boil;
        t_boil
    } else {
        state.last_tboil
    };
    let mut wbt = tdb;
    let mut x1 = 0.0;
    let mut y1 = 0.0;
    for iter in 1..=ITMAX {
        if wbt >= t_boil - 0.09 {
            wbt = t_boil - 0.1;
        }
        let p_sat_star = psy_psat_fn_temp(wbt);
        let w_star = 0.62198 * p_sat_star / (patm - p_sat_star);
        let new_w = if wbt >= 0.0 {
            ((2501.0 - 2.326 * wbt) * w_star - 1.006 * (tdb - wbt))
                / (2501.0 + 1.86 * tdb - 4.186 * wbt)
        } else {
            ((2830.0 - 0.24 * wbt) * w_star - 1.006 * (tdb - wbt))
                / (2830.0 + 1.86 * tdb - 2.1 * wbt)
        };
        let error = w - new_w;
        let (result_x, converged) = iterate(wbt, error, state.iconv_tol, iter, &mut x1, &mut y1);
        wbt = result_x;
        if converged {
            break;
        }
    }
    if wbt > tdb {
        tdb
    } else {
        wbt
    }
}

fn f6(x: f64, a0: f64, a1: f64, a2: f64, a3: f64, a4: f64, a5: f64) -> f64 {
    a0 + x * (a1 + x * (a2 + x * (a3 + x * (a4 + x * a5))))
}

fn f7(x: f64, a0: f64, a1: f64, a2: f64, a3: f64, a4: f64, a5: f64, a6: f64) -> f64 {
    (a0 + x * (a1 + x * (a2 + x * (a3 + x * (a4 + x * (a5 + x * a6)))))) / 1.0e10
}

/// Saturation temperature `{C}` from moist-air enthalpy `{J/kg}` and pressure `{Pa}`.
pub fn psy_tsat_fn_h_pb(state: &mut PsychState, h: f64, pb: f64) -> f64 {
    let mut hh = h + 1.78637e4;
    let hloc = if h >= 0.0 {
        h.max(0.00001)
    } else {
        h.min(-0.00001)
    };
    let case_range = [
        -4.24e4, -2.2138e4, -6.7012e2, 2.7297e4, 7.5222e4, 1.8379e5, 4.7577e5, 1.5445e6, 3.8353e6,
        4.5866e7,
    ];
    let mut beg = 0i32;
    let mut end = 9i32;
    while beg + 1 < end {
        let mid = (beg + end) >> 1;
        if hh > case_range[mid as usize] {
            beg = mid;
        } else {
            end = mid;
        }
    }
    let case_index = beg + 1;
    let mut t = match case_index {
        1 => {
            if hh < -4.24e4 {
                hh = -4.24e4;
            }
            f6(
                hh,
                -19.44,
                8.53675e-4,
                -5.12637e-9,
                -9.85546e-14,
                -1.00102e-18,
                -4.2705e-24,
            )
        }
        2 => f6(
            hh,
            -1.94224e1,
            8.5892e-4,
            -4.50709e-9,
            -6.19492e-14,
            8.71734e-20,
            8.73051e-24,
        ),
        3 => f6(
            hh,
            -1.94224e1,
            8.59061e-4,
            -4.4875e-9,
            -5.76696e-14,
            7.72217e-19,
            3.97894e-24,
        ),
        4 => f6(
            hh,
            -2.01147e1,
            9.04936e-4,
            -6.83305e-9,
            2.3261e-14,
            7.27237e-20,
            -6.31939e-25,
        ),
        5 => f6(
            hh,
            -1.82124e1,
            8.31683e-4,
            -6.16461e-9,
            3.06411e-14,
            -8.60964e-20,
            1.03003e-25,
        ),
        6 => f6(
            hh,
            -1.29419,
            3.88538e-4,
            -1.30237e-9,
            2.78254e-15,
            -3.27225e-21,
            1.60969e-27,
        ),
        7 => f6(
            hh,
            2.39214e1,
            1.27519e-4,
            -1.52089e-10,
            1.1043e-16,
            -4.33919e-23,
            7.05296e-30,
        ),
        8 => f6(
            hh,
            4.88446e1,
            3.85534e-5,
            -1.78805e-11,
            4.87224e-18,
            -7.15283e-25,
            4.36246e-32,
        ),
        _ => {
            if hh > 4.5866e7 {
                hh = 4.5866e7;
            }
            f7(
                hh,
                7.60565e11,
                5.80534e4,
                -7.36433e-3,
                5.11531e-10,
                -1.93619e-17,
                3.70511e-25,
                -2.77313e-33,
            )
        }
    };
    if (pb - 1.0133e5).abs() / 1.0133e5 > 0.01 {
        let mut t1 = t;
        let h1 = psy_h_fn_tdb_w(t1, psy_w_fn_tdb_twb_pb(t1, t1, pb));
        let mut y1 = h1 - hloc;
        if (y1 / hloc).abs() <= 0.1e-4 {
            t = t1;
        } else {
            let mut t2 = t1 * 0.9;
            let mut iter_count = 0;
            while iter_count <= 30 {
                iter_count += 1;
                let h2 = psy_h_fn_tdb_w(t2, psy_w_fn_tdb_twb_pb(t2, t2, pb));
                let y2 = h2 - hloc;
                if (y2 / hloc).abs() <= 0.1e-4 || y2 == y1 {
                    t = t2;
                    break;
                }
                let tn = t2 - y2 / (y2 - y1) * (t2 - t1);
                t1 = t2;
                t2 = tn;
                y1 = y2;
            }
        }
    }
    let _ = state;
    t
}

/// Dew-point `{C}` from humidity ratio and pressure.
pub fn psy_tdp_fn_w_pb(state: &mut PsychState, w: f64, pb: f64) -> f64 {
    let w0 = hum_ratio(w);
    let pdew = pb * w0 / (0.62198 + w0);
    psy_tsat_fn_pb(state, pdew)
}

/// Dew-point `{C}` from dry-bulb, wet-bulb, and pressure.
pub fn psy_tdp_fn_tdb_twb_pb(state: &mut PsychState, tdb: f64, twb: f64, pb: f64) -> f64 {
    let w = psy_w_fn_tdb_twb_pb(tdb, twb, pb).max(1.0e-5);
    let tdp = psy_tdp_fn_w_pb(state, w, pb);
    if tdp > twb {
        twb
    } else {
        tdp
    }
}

/// Vapor density from dry-bulb and relative humidity, continuous across 0 C.
pub fn psy_rhov_fn_tdb_rh(tdb: f64, rh: f64) -> f64 {
    psy_psat_fn_temp(tdb) * rh / (461.52 * (tdb + KELVIN))
}

/// Relative humidity from dry-bulb and vapor density, continuous across 0 C.
pub fn psy_rh_fn_tdb_rhov(tdb: f64, rhovapor: f64) -> f64 {
    let rh = if rhovapor > 0.0 {
        rhovapor * 461.52 * (tdb + KELVIN) / psy_psat_fn_temp(tdb)
    } else {
        0.0
    };
    if !(0.0..=1.0).contains(&rh) {
        rh.clamp(0.01, 1.0)
    } else {
        rh
    }
}

/// Specific heat of chilled or hot water `{J/kg-K}`. EnergyPlus returns a constant.
pub fn cp_water(_temperature: f64) -> f64 {
    4180.0
}

/// Density of liquid water `{kg/m3}` from dry-bulb `{C}`.
pub fn rho_h2o(tb: f64) -> f64 {
    1000.1207 + 8.3215874e-04 * tb - 4.929976e-03 * tb * tb + 8.4791863e-06 * tb * tb * tb
}

/// Sensible enthalpy difference `{J/kg}` at a single humidity ratio.
pub fn psy_delta_h_sen_fn_tdb2_tdb1_w(tdb2: f64, tdb1: f64, w: f64) -> f64 {
    (1.00484e3 + hum_ratio(w) * 1.85895e3) * (tdb2 - tdb1)
}

/// Sensible enthalpy difference `{J/kg}` using the smaller of two humidity ratios.
pub fn psy_delta_h_sen_fn_tdb2_w2_tdb1_w1(tdb2: f64, w2: f64, tdb1: f64, w1: f64) -> f64 {
    psy_delta_h_sen_fn_tdb2_tdb1_w(tdb2, tdb1, w1.min(w2))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn saturation_pressure_matches_hyland_wexler_reference_points() {
        let cases = [
            (-10.0, 259.902865),
            (0.0, 611.1535709),
            (0.01, 611.6570244),
            (20.0, 2338.8037),
            (25.0, 3169.21647),
            (100.0, 101418.7168),
        ];
        for (t, expected) in cases {
            let p = psy_psat_fn_temp(t);
            assert!(
                (p - expected).abs() / expected < 1.0e-6,
                "psat({t}) = {p}, expected {expected}"
            );
        }
    }

    #[test]
    fn saturation_temperature_inverts_saturation_pressure() {
        let mut state = PsychState::default();
        for t in [-20.0, 0.5, 20.0, 35.0, 60.0, 99.0] {
            let p = psy_psat_fn_temp(t);
            let back = psy_tsat_fn_pb(&mut state, p);
            assert!((back - t).abs() < 0.02, "T {t} -> P {p} -> T {back}");
        }
    }

    #[test]
    fn enthalpy_round_trip_and_density() {
        let h = psy_h_fn_tdb_w(25.0, 0.009);
        assert!((h - 48047.72375).abs() < 1.0e-3);
        let t = psy_tdb_fn_hw(h, 0.009);
        assert!((t - 25.0).abs() < 1.0e-9);
        let rho = psy_rho_air_fn_pb_tdb_w(101325.0, 20.0, 0.007);
        assert!((rho - 1.19092497).abs() < 1.0e-6);
    }

    #[test]
    fn humidity_and_relative_humidity_round_trip() {
        let mut state = PsychState::default();
        let pb = 101325.0;
        let w = psy_w_fn_tdb_rh_pb(22.0, 0.45, pb);
        let rh = psy_rh_fn_tdb_w_pb(&mut state, 22.0, w, pb);
        assert!((rh - 0.45).abs() < 1.0e-6, "rh {rh}");
        let twb = psy_twb_fn_tdb_w_pb(&mut state, 22.0, w, pb);
        assert!(twb < 22.0 && twb > 10.0, "twb {twb}");
        let w_sat = psy_w_fn_tdb_rh_pb(22.0, 1.0, pb);
        let twb_sat = psy_twb_fn_tdb_w_pb(&mut state, 22.0, w_sat, pb);
        assert!((twb_sat - 22.0).abs() < 0.05, "saturated twb {twb_sat}");
    }

    #[test]
    fn sensible_delta_matches_enthalpy_difference() {
        let w = 0.008;
        let direct = psy_h_fn_tdb_w(30.0, w) - psy_h_fn_tdb_w(18.0, w);
        let delta = psy_delta_h_sen_fn_tdb2_tdb1_w(30.0, 18.0, w);
        assert!((direct - delta).abs() < 1.0e-6);
    }
}
