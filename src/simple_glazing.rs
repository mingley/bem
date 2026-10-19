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
    }
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
