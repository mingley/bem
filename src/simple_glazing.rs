//! `WindowMaterial:SimpleGlazingSystem` block model from EnergyPlus `MaterialGlass::SetupSimpleWindowGlazingSystem`.
//!
//! U-factor and SHGC become a normal solar transmittance and the inward-flowing
//! portion of the layer absorptance. Their sum is the original SHGC.

#[derive(Clone, Copy, Debug)]
pub struct SimpleGlazing {
    pub transmittance: f64,
    pub solar_absorptance: f64,
    pub inward_fraction: f64,
    pub reflectance: f64,
    /// Glazing resistance with the standard film coefficients removed, m²K/W.
    pub resistance: f64,
}

impl SimpleGlazing {
    /// Absorptance that flows into the zone. Transmittance plus this value is SHGC.
    pub fn inward_absorptance(self) -> f64 {
        self.inward_fraction * self.solar_absorptance
    }
}

pub fn simple_glazing(u_factor: f64, shgc: f64) -> SimpleGlazing {
    if u_factor <= 0.0 {
        return SimpleGlazing {
            transmittance: 0.0,
            solar_absorptance: 0.0,
            inward_fraction: 0.0,
            reflectance: 1.0,
            resistance: 1.0,
        };
    }
    let riw = if u_factor < 5.85 {
        1.0 / (0.359073 * u_factor.ln() + 6.949915)
    } else {
        1.0 / (1.788041 * u_factor - 2.886625)
    };
    let row = 1.0 / (0.025342 * u_factor + 29.163853);
    let rlw = (1.0 / u_factor - riw - row).max(0.001);

    let transmittance = solar_transmittance(u_factor, shgc).max(0.0);
    let delta = shgc - transmittance;
    let (ris, ros) = summer_film_resistances(u_factor, delta);
    let denom = ros + rlw + ris;
    let inward_fraction = if denom.abs() > 1.0e-12 {
        (ros + 0.5 * rlw) / denom
    } else {
        0.0
    };
    let solar_absorptance = if inward_fraction.abs() > 1.0e-12 {
        delta / inward_fraction
    } else {
        0.0
    };
    let solar_absorptance = solar_absorptance.max(0.0);
    let reflectance = (1.0 - transmittance - solar_absorptance).max(0.0);
    SimpleGlazing {
        transmittance,
        solar_absorptance,
        inward_fraction,
        reflectance,
        resistance: rlw,
    }
}

/// Inside and outside face temperatures of a single simple-glazing layer.
///
/// Absorbed solar is split equally between the two faces, matching
/// `CalcWindowHeatBalance`. `q_to_inside` is extra heat gain on the room face,
/// W/m², positive into the glass.
pub fn window_face_temperatures(
    resistance: f64,
    h_outside: f64,
    t_outside: f64,
    q_absorbed_outside: f64,
    h_inside: f64,
    t_zone: f64,
    q_to_inside: f64,
) -> (f64, f64) {
    let g = 1.0 / resistance.max(0.001);
    let h_out = h_outside.max(0.0);
    let h_in = h_inside.max(0.0);
    let rhs_out = h_out * t_outside + q_absorbed_outside;
    let rhs_in = h_in * t_zone + q_to_inside;
    let det = (h_out + g) * (g + h_in) - g * g;
    if det.abs() < 1.0e-12 {
        return (t_outside, t_zone);
    }
    let t_so = (rhs_out * (g + h_in) + g * rhs_in) / det;
    let t_si = ((h_out + g) * rhs_in + g * rhs_out) / det;
    (t_so, t_si)
}

/// Beam solar properties at incidence cosine `cos_inc`, from the simple-glazing
/// angular curves in `WindowManager::TransAndReflAtPhi`.
///
/// `cos_inc` is 1 at normal incidence. The returned transmittance and reflectance
/// are normalized to the normal-incidence values. Diffuse solar keeps those
/// normal values.
pub fn simple_glazing_beam(
    u_factor: f64,
    shgc: f64,
    normal_transmittance: f64,
    normal_reflectance: f64,
    cos_inc: f64,
) -> (f64, f64, f64) {
    let cs = cos_inc.clamp(0.0, 1.0);
    let (mut trans_factor, mut reflect_factor) = angular_curve_factors(u_factor, shgc, cs);
    if cs == 1.0 {
        trans_factor = 1.0;
        reflect_factor = 0.0;
    }
    let transmittance = (normal_transmittance * trans_factor).clamp(0.0, 1.0);
    let reflectance = (normal_reflectance * (1.0 - reflect_factor) + reflect_factor)
        .clamp(0.0, (0.9999 - transmittance).max(0.0));
    let absorptance = (1.0 - transmittance - reflectance).max(0.0);
    (transmittance, reflectance, absorptance)
}

#[allow(unused_assignments)]
fn angular_curve_factors(u: f64, shgc: f64, cs: f64) -> (f64, f64) {
    let cs2 = cs * cs;
    let cs3 = cs2 * cs;
    let cs4 = cs2 * cs2;
    let trans = |b: f64, c: f64, d: f64, e: f64| b * cs + c * cs2 + d * cs3 + e * cs4;
    let trans_a = trans(3.36, -3.85, 1.49, 0.01);
    let trans_b = trans(2.83, -2.42, 0.04, 0.55);
    let trans_c = trans(2.45, -1.58, -0.64, 0.77);
    let trans_d = trans(2.85, -2.58, 0.40, 0.35);
    let trans_e = trans(1.51, 2.49, -5.87, 2.88);
    let trans_f = trans(1.21, 3.14, -6.37, 3.03);
    let trans_g = trans(1.09, 3.54, -6.84, 3.23);
    let trans_h = trans(0.98, 3.83, -7.13, 3.33);
    let trans_i = trans(0.79, 3.93, -6.86, 3.15);
    let trans_j = trans(0.08, 6.02, -8.84, 3.74);
    let trans_fghi = (trans_f + trans_g + trans_h + trans_i) / 4.0;
    let trans_fh = (trans_f + trans_h) / 2.0;
    let trans_bdcd = (trans_b + trans_d + trans_c + trans_d) / 4.0;
    let reflect = |b: f64, c: f64, d: f64, e: f64, curve: f64| {
        1.0 + b * cs + c * cs2 + d * cs3 + e * cs4 - curve
    };
    let reflect_a = reflect(-0.70, 2.57, -3.20, 1.33, trans_a);
    let reflect_b = reflect(-1.87, 6.50, -7.86, 3.23, trans_b);
    let reflect_c = reflect(-2.52, 8.40, -9.86, 3.99, trans_c);
    let reflect_d = reflect(-1.85, 6.40, -7.64, 3.11, trans_d);
    let reflect_e = reflect(-1.57, 5.60, -6.82, 2.80, trans_e);
    let reflect_f = reflect(-3.15, 10.98, -13.14, 5.32, trans_f);
    let reflect_g = reflect(-3.25, 11.32, -13.54, 5.49, trans_g);
    let reflect_h = reflect(-3.39, 11.70, -13.94, 5.64, trans_h);
    let reflect_i = reflect(-4.06, 13.55, -15.74, 6.27, trans_i);
    let reflect_j = reflect(-4.35, 14.27, -16.32, 6.39, trans_j);
    let reflect_fghi = (reflect_f + reflect_g + reflect_h + reflect_i) / 4.0;
    let reflect_fh = (reflect_f + reflect_h) / 2.0;
    let reflect_bdcd = (reflect_b + reflect_d + reflect_c + reflect_d) / 4.0;

    let mut trans_tmp = 0.0;
    let mut reflect_tmp = 0.0;
    if u < 1.4195 {
        if shgc > 0.45 {
            trans_tmp = trans_e;
            reflect_tmp = reflect_e;
        } else if (0.35..0.45).contains(&shgc) || shgc == 0.45 {
            trans_tmp = lerp(shgc, 0.35, 0.45, trans_j, trans_e);
            reflect_tmp = lerp(shgc, 0.35, 0.45, reflect_j, reflect_e);
        } else {
            trans_tmp = trans_j;
            reflect_tmp = reflect_j;
        }
    } else if (1.4195..=1.7034).contains(&u) {
        if shgc > 0.55 {
            trans_tmp = trans_e;
            reflect_tmp = reflect_e;
        } else if shgc > 0.5 && shgc <= 0.55 {
            trans_tmp = bilerp(
                u, shgc, 1.4195, 1.7034, 0.50, 0.55, trans_e, trans_e, trans_fghi, trans_e,
            );
            reflect_tmp = bilerp(
                u,
                shgc,
                1.4195,
                1.7034,
                0.50,
                0.55,
                reflect_e,
                reflect_e,
                reflect_fghi,
                reflect_e,
            );
        } else if shgc > 0.45 && shgc <= 0.5 {
            trans_tmp = lerp(u, 1.4195, 1.7034, trans_e, trans_fghi);
            reflect_tmp = lerp(u, 1.4195, 1.7034, reflect_e, reflect_fghi);
        } else if shgc > 0.35 && shgc <= 0.45 {
            trans_tmp = bilerp(
                u, shgc, 1.4195, 1.7034, 0.35, 0.45, trans_j, trans_e, trans_fghi, trans_fghi,
            );
            reflect_tmp = bilerp(
                u,
                shgc,
                1.4195,
                1.7034,
                0.35,
                0.45,
                reflect_j,
                reflect_e,
                reflect_fghi,
                reflect_fghi,
            );
        } else if shgc > 0.3 && shgc <= 0.35 {
            trans_tmp = lerp(u, 1.4195, 1.7034, trans_j, trans_fghi);
            reflect_tmp = lerp(u, 1.4195, 1.7034, reflect_j, reflect_fghi);
        } else if shgc > 0.25 && shgc <= 0.3 {
            trans_tmp = bilerp(
                u, shgc, 1.4195, 1.7034, 0.25, 0.3, trans_j, trans_j, trans_fh, trans_fghi,
            );
            reflect_tmp = bilerp(
                u,
                shgc,
                1.4195,
                1.7034,
                0.25,
                0.3,
                reflect_j,
                reflect_j,
                reflect_fh,
                reflect_fghi,
            );
        } else {
            trans_tmp = lerp(u, 1.4195, 1.7034, trans_j, trans_fh);
            reflect_tmp = lerp(u, 1.4195, 1.7034, reflect_j, reflect_fh);
        }
    } else if u > 1.7034 && u < 3.4068 {
        if shgc > 0.55 {
            trans_tmp = trans_e;
            reflect_tmp = reflect_e;
        } else if (0.5..=0.55).contains(&shgc) {
            trans_tmp = lerp(shgc, 0.5, 0.55, trans_fghi, trans_e);
            reflect_tmp = lerp(shgc, 0.5, 0.55, reflect_fghi, reflect_e);
        } else if shgc > 0.3 && shgc < 0.5 {
            trans_tmp = trans_fghi;
            reflect_tmp = reflect_fghi;
        } else if (0.25..=0.3).contains(&shgc) {
            trans_tmp = lerp(shgc, 0.25, 0.30, trans_fh, trans_fghi);
            reflect_tmp = lerp(shgc, 0.25, 0.30, reflect_fh, reflect_fghi);
        } else {
            trans_tmp = trans_fh;
            reflect_tmp = reflect_fh;
        }
    } else if (3.4068..=4.5424).contains(&u) {
        if shgc > 0.65 {
            trans_tmp = lerp(u, 3.4068, 4.5424, trans_e, trans_a);
            reflect_tmp = lerp(u, 3.4068, 4.5424, reflect_e, reflect_a);
        } else if shgc > 0.6 && shgc <= 0.65 {
            trans_tmp = bilerp(
                u, shgc, 3.4068, 4.5424, 0.6, 0.65, trans_e, trans_e, trans_bdcd, trans_a,
            );
            reflect_tmp = bilerp(
                u,
                shgc,
                3.4068,
                4.5424,
                0.6,
                0.65,
                reflect_e,
                reflect_e,
                reflect_bdcd,
                reflect_a,
            );
        } else if shgc > 0.55 && shgc <= 0.6 {
            trans_tmp = lerp(u, 3.4068, 4.5424, trans_e, trans_bdcd);
            reflect_tmp = lerp(u, 3.4068, 4.5424, reflect_e, reflect_bdcd);
        } else if shgc > 0.5 && shgc <= 0.55 {
            trans_tmp = bilerp(
                u, shgc, 3.4068, 4.5424, 0.5, 0.55, trans_fghi, trans_e, trans_bdcd, trans_bdcd,
            );
            reflect_tmp = bilerp(
                u,
                shgc,
                3.4068,
                4.5424,
                0.5,
                0.55,
                reflect_fghi,
                reflect_e,
                reflect_bdcd,
                reflect_bdcd,
            );
        } else if shgc > 0.45 && shgc <= 0.5 {
            trans_tmp = lerp(u, 3.4068, 4.5424, trans_fghi, trans_bdcd);
            reflect_tmp = lerp(u, 3.4068, 4.5424, reflect_fghi, reflect_bdcd);
        } else if shgc > 0.3 && shgc <= 0.45 {
            trans_tmp = bilerp(
                u, shgc, 3.4068, 4.5424, 0.3, 0.45, trans_fghi, trans_fghi, trans_d, trans_bdcd,
            );
            reflect_tmp = bilerp(
                u,
                shgc,
                3.4068,
                4.5424,
                0.3,
                0.45,
                reflect_fghi,
                reflect_fghi,
                reflect_d,
                reflect_bdcd,
            );
        } else if shgc > 0.25 && shgc <= 0.3 {
            trans_tmp = bilerp(
                u, shgc, 3.4068, 4.5424, 0.25, 0.3, trans_fh, trans_fghi, trans_d, trans_d,
            );
            reflect_tmp = bilerp(
                u,
                shgc,
                3.4068,
                4.5424,
                0.25,
                0.3,
                reflect_fh,
                reflect_fghi,
                reflect_d,
                reflect_d,
            );
        } else {
            trans_tmp = lerp(u, 3.4068, 4.5424, trans_fh, trans_d);
            reflect_tmp = lerp(u, 3.4068, 4.5424, reflect_fh, reflect_d);
        }
    } else if shgc > 0.65 {
        trans_tmp = trans_a;
        reflect_tmp = reflect_a;
    } else if (0.6..=0.65).contains(&shgc) {
        trans_tmp = lerp(shgc, 0.6, 0.65, trans_bdcd, trans_a);
        reflect_tmp = lerp(shgc, 0.6, 0.65, reflect_bdcd, reflect_a);
    } else if shgc > 0.45 && shgc < 0.6 {
        trans_tmp = trans_bdcd;
        reflect_tmp = reflect_bdcd;
    } else if (0.3..=0.45).contains(&shgc) {
        trans_tmp = lerp(shgc, 0.3, 0.45, trans_d, trans_bdcd);
        reflect_tmp = lerp(shgc, 0.3, 0.45, reflect_d, reflect_bdcd);
    } else {
        trans_tmp = trans_d;
        reflect_tmp = reflect_d;
    }
    (trans_tmp, reflect_tmp)
}

fn lerp(x: f64, x0: f64, x1: f64, f0: f64, f1: f64) -> f64 {
    f0 + ((x - x0) / (x1 - x0)) * (f1 - f0)
}

fn bilerp(
    x: f64,
    y: f64,
    x1: f64,
    x2: f64,
    y1: f64,
    y2: f64,
    fx1y1: f64,
    fx1y2: f64,
    fx2y1: f64,
    fx2y2: f64,
) -> f64 {
    let denom = (x2 - x1) * (y2 - y1);
    fx1y1 / denom * (x2 - x) * (y2 - y)
        + fx2y1 / denom * (x - x1) * (y2 - y)
        + fx1y2 / denom * (x2 - x) * (y - y1)
        + fx2y2 / denom * (x - x1) * (y - y1)
}

fn solar_transmittance(u_factor: f64, shgc: f64) -> f64 {
    let high = || {
        if shgc < 0.7206 {
            0.939998 * shgc * shgc + 0.20332 * shgc
        } else {
            1.30415 * shgc - 0.30515
        }
    };
    let low = || {
        if shgc <= 0.15 {
            0.41040 * shgc
        } else {
            0.085775 * shgc * shgc + 0.963954 * shgc - 0.084958
        }
    };
    if u_factor > 4.5 {
        high()
    } else if u_factor < 3.4 {
        low()
    } else {
        ((u_factor - 3.4) / (4.5 - 3.4)) * (high() - low()) + low()
    }
}

fn summer_film_resistances(u_factor: f64, delta: f64) -> (f64, f64) {
    let ris_high =
        1.0 / (29.436546 * delta.powi(3) - 21.943415 * delta.powi(2) + 9.945872 * delta + 7.426151);
    let ros_high = 1.0 / (2.225824 * delta + 20.577080);
    let ris_low = 1.0
        / (199.8208128 * delta.powi(3) - 90.639733 * delta.powi(2) + 19.737055 * delta + 6.766575);
    let ros_low = 1.0 / (5.763355 * delta + 20.541528);
    if u_factor > 4.5 {
        (ris_high, ros_high)
    } else if u_factor < 3.4 {
        (ris_low, ros_low)
    } else {
        let weight = (u_factor - 3.4) / (4.5 - 3.4);
        (
            weight * (ris_low - ris_high) + ris_low,
            weight * (ros_low - ros_high) + ros_low,
        )
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SolarReceiver {
    pub area: f64,
    pub inside_absorptance: f64,
    /// Outward normal points downward (CosTilt < -0.5), so the surface is a floor.
    pub floor: bool,
}

/// Watts absorbed by each receiver from transmitted beam and diffuse solar.
///
/// Beam falls on floors in proportion to area times inside absorptance over the
/// floor area (`SurfIntAbsFac`). Whatever the floor does not absorb, plus the
/// transmitted diffuse, is spread over all receivers by area times absorptance.
pub fn absorb_transmitted(receivers: &[SolarReceiver], beam_w: f64, diffuse_w: f64) -> Vec<f64> {
    let mut absorbed = vec![0.0; receivers.len()];
    let floor_area: f64 = receivers.iter().filter(|r| r.floor).map(|r| r.area).sum();
    let mut reflected_beam = beam_w;
    if floor_area > 0.0 {
        for (i, receiver) in receivers.iter().enumerate() {
            if receiver.floor {
                let portion = beam_w * receiver.area * receiver.inside_absorptance / floor_area;
                absorbed[i] += portion;
                reflected_beam -= portion;
            }
        }
    }
    let pool = diffuse_w + reflected_beam.max(0.0);
    let weight: f64 = receivers
        .iter()
        .map(|r| r.area * r.inside_absorptance.max(0.0))
        .sum();
    if weight > 0.0 {
        for (i, receiver) in receivers.iter().enumerate() {
            absorbed[i] += pool * receiver.area * receiver.inside_absorptance.max(0.0) / weight;
        }
    }
    absorbed
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn block_model_reproduces_shgc() {
        let cases = [
            (2.2, 0.4, 0.3143476),
            (5.0, 0.6, 0.46039128),
            (4.0, 0.5, 0.373842795),
        ];
        for (u, shgc, transmittance) in cases {
            let glass = simple_glazing(u, shgc);
            assert!(
                (glass.transmittance - transmittance).abs() < 1.0e-8,
                "U={u} T={} expected {transmittance}",
                glass.transmittance
            );
            let rebuilt = glass.transmittance + glass.inward_absorptance();
            assert!((rebuilt - shgc).abs() < 1.0e-8, "SHGC {rebuilt} vs {shgc}");
        }
    }

    #[test]
    fn rating_films_recover_the_u_factor() {
        let u = 2.2;
        let glass = simple_glazing(u, 0.4);
        let h_in = if u < 5.85 {
            0.359073 * u.ln() + 6.949915
        } else {
            1.788041 * u - 2.886625
        };
        let h_out = 0.025342 * u + 29.163853;
        let (t_so, t_si) =
            window_face_temperatures(glass.resistance, h_out, 0.0, 0.0, h_in, 20.0, 0.0);
        let flux = h_in * (t_si - 20.0);
        let expected = u * (0.0 - 20.0);
        assert!(
            (flux - expected).abs() < 1.0e-6,
            "flux {flux} expected {expected}, faces {t_so} {t_si}"
        );
    }

    #[test]
    fn normal_incidence_keeps_the_block_model_transmittance() {
        let glass = simple_glazing(2.2, 0.4);
        let (transmittance, reflectance, absorptance) =
            simple_glazing_beam(2.2, 0.4, glass.transmittance, glass.reflectance, 1.0);
        assert!((transmittance - glass.transmittance).abs() < 1.0e-12);
        assert!((reflectance - glass.reflectance).abs() < 1.0e-12);
        assert!((absorptance - glass.solar_absorptance).abs() < 1.0e-9);
    }

    #[test]
    fn oblique_beam_transmits_less_than_normal_incidence() {
        let glass = simple_glazing(2.2, 0.4);
        let (transmittance, _, _) =
            simple_glazing_beam(2.2, 0.4, glass.transmittance, glass.reflectance, 0.5);
        assert!(transmittance < glass.transmittance);
        assert!(transmittance > 0.0);
    }

    #[test]
    fn beam_is_absorbed_on_the_floor_and_the_reflection_is_redistributed() {
        let receivers = [
            SolarReceiver {
                area: 10.0,
                inside_absorptance: 0.5,
                floor: true,
            },
            SolarReceiver {
                area: 10.0,
                inside_absorptance: 0.5,
                floor: false,
            },
        ];
        let absorbed = absorb_transmitted(&receivers, 100.0, 0.0);
        // Floor takes 100 * 10 * 0.5 / 10 = 50 W. The other 50 W is split by equal area*alpha.
        assert!((absorbed[0] - 75.0).abs() < 1.0e-9, "{}", absorbed[0]);
        assert!((absorbed[1] - 25.0).abs() < 1.0e-9, "{}", absorbed[1]);
        assert!((absorbed.iter().sum::<f64>() - 100.0).abs() < 1.0e-9);
    }
}
