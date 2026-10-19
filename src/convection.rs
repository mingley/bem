//! Surface convection coefficients translated from EnergyPlus `ConvectionCoefficients`.
//!
//! Interior uses `CalcASHRAETARPNatural` (the TARP default). Exterior uses the DOE-2
//! model: TARP natural plus a roughness-adjusted MoWiTT forced term. The coefficient
//! is floored at EnergyPlus `LowHConvLimit` (0.1 W/m²K).

/// EnergyPlus `Material::SurfaceRoughness`, in the order of `RoughnessMultiplier`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Roughness {
    VeryRough = 0,
    Rough = 1,
    MediumRough = 2,
    MediumSmooth = 3,
    Smooth = 4,
    VerySmooth = 5,
}

impl Roughness {
    pub fn parse(text: &str) -> Self {
        match text.trim().to_ascii_lowercase().as_str() {
            "veryrough" => Self::VeryRough,
            "rough" => Self::Rough,
            "mediumrough" => Self::MediumRough,
            "mediumsmooth" => Self::MediumSmooth,
            "smooth" => Self::Smooth,
            "verysmooth" => Self::VerySmooth,
            _ => Self::MediumRough,
        }
    }

    /// `RoughnessMultiplier` in `ConvectionCoefficients.cc`.
    pub fn multiplier(self) -> f64 {
        const MULTIPLIER: [f64; 6] = [2.17, 1.67, 1.52, 1.13, 1.11, 1.0];
        MULTIPLIER[self as usize]
    }
}

pub const LOW_H_CONV: f64 = 0.1;

/// `CalcASHRAEVerticalWall`.
pub fn ashrae_vertical_wall(delta_temp: f64) -> f64 {
    1.31 * delta_temp.abs().powf(1.0 / 3.0)
}

/// `CalcWaltonUnstableHorizontalOrTilt`.
pub fn walton_unstable(delta_temp: f64, cosine_tilt: f64) -> f64 {
    9.482 * delta_temp.abs().powf(1.0 / 3.0) / (7.238 - cosine_tilt.abs())
}

/// `CalcWaltonStableHorizontalOrTilt`.
pub fn walton_stable(delta_temp: f64, cosine_tilt: f64) -> f64 {
    1.810 * delta_temp.abs().powf(1.0 / 3.0) / (1.382 + cosine_tilt.abs())
}

/// `CalcASHRAETARPNatural`. `cos_tilt` is the cosine of the outward tilt.
///
/// Interior calls pass the negated cosine, matching EnergyPlus.
pub fn ashrae_tarp_natural(t_surf: f64, t_air: f64, cos_tilt: f64) -> f64 {
    let delta_temp = t_surf - t_air;
    if delta_temp == 0.0 || cos_tilt == 0.0 {
        return ashrae_vertical_wall(delta_temp);
    }
    if (delta_temp < 0.0 && cos_tilt < 0.0) || (delta_temp > 0.0 && cos_tilt > 0.0) {
        walton_unstable(delta_temp, cos_tilt)
    } else {
        walton_stable(delta_temp, cos_tilt)
    }
}

/// Interior film. `cos_tilt_outward` is the outward normal's z component.
pub fn tarp_interior_hc(t_surf: f64, t_air: f64, cos_tilt_outward: f64) -> f64 {
    ashrae_tarp_natural(t_surf, t_air, -cos_tilt_outward).max(LOW_H_CONV)
}

/// `CalcMoWITTNatural`.
pub fn mowitt_natural(delta_temp: f64) -> f64 {
    0.84 * delta_temp.abs().powf(1.0 / 3.0)
}

/// `CalcMoWITTForcedWindward`.
pub fn mowitt_forced_windward(wind_at_z: f64) -> f64 {
    3.26 * wind_at_z.max(0.0).powf(0.89)
}

/// `CalcMoWITTForcedLeeward`.
pub fn mowitt_forced_leeward(wind_at_z: f64) -> f64 {
    3.55 * wind_at_z.max(0.0).powf(0.617)
}

/// `Windward`. Leeward is more than about 90° from the outward normal.
/// A nearly horizontal surface is always windward.
pub fn windward(cos_tilt: f64, azimuth_deg: f64, wind_direction_deg: f64) -> bool {
    if cos_tilt.abs() >= 0.98 {
        return true;
    }
    let mut diff = (wind_direction_deg - azimuth_deg).abs();
    if diff - 180.0 > 0.001 {
        diff -= 360.0;
    }
    diff.abs() - 90.0 <= 0.001
}

/// `CalcDOE2Forced`: roughness-adjusted forced coefficient, not yet added to natural.
pub fn doe2_forced(
    t_surf: f64,
    t_air: f64,
    cos_tilt: f64,
    hf_smooth: f64,
    roughness: Roughness,
) -> f64 {
    let hn = ashrae_tarp_natural(t_surf, t_air, cos_tilt);
    let hc_smooth = (hn * hn + hf_smooth * hf_smooth).sqrt();
    roughness.multiplier() * (hc_smooth - hn)
}

/// DOE-2 outside convection, `Hn + Hf`, floored at `LowHConvLimit`.
pub fn doe2_exterior_hc(
    t_surf: f64,
    t_air: f64,
    cos_tilt: f64,
    wind_speed: f64,
    wind_direction_deg: f64,
    azimuth_deg: f64,
    roughness: Roughness,
) -> f64 {
    let hf_smooth = if windward(cos_tilt, azimuth_deg, wind_direction_deg) {
        mowitt_forced_windward(wind_speed)
    } else {
        mowitt_forced_leeward(wind_speed)
    };
    let hf = doe2_forced(t_surf, t_air, cos_tilt, hf_smooth, roughness);
    let hn = ashrae_tarp_natural(t_surf, t_air, cos_tilt);
    (hn + hf).max(LOW_H_CONV)
}

/// Azimuth of an outward normal, degrees clockwise from north. East is 90.
pub fn azimuth_from_normal(normal: [f64; 3]) -> f64 {
    let mut deg = normal[0].atan2(normal[1]).to_degrees();
    if deg < 0.0 {
        deg += 360.0;
    }
    deg
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tarp_matches_ashrae_and_walton_equations() {
        let vertical = ashrae_tarp_natural(28.0, 20.0, 0.0);
        assert!((vertical - 2.62).abs() < 1.0e-12, "{vertical}");
        let unstable = ashrae_tarp_natural(24.0, 20.0, 1.0);
        assert!(
            (unstable - 2.4129106724531044).abs() < 1.0e-12,
            "{unstable}"
        );
        let stable = ashrae_tarp_natural(16.0, 20.0, 1.0);
        assert!((stable - 1.20621154662571).abs() < 1.0e-12, "{stable}");
    }

    #[test]
    fn interior_negates_outward_cosine_and_floors_zero_delta() {
        let hc = tarp_interior_hc(20.0, 20.0, 0.0);
        assert!((hc - LOW_H_CONV).abs() < 1.0e-12);
        // Outward floor normal points down. Interior cosine is +1, surface warmer: unstable.
        let floor = tarp_interior_hc(24.0, 20.0, -1.0);
        let expected = ashrae_tarp_natural(24.0, 20.0, 1.0);
        assert!((floor - expected).abs() < 1.0e-12);
    }

    #[test]
    fn doe2_windward_very_smooth_is_the_smooth_combination() {
        let total = doe2_exterior_hc(28.0, 20.0, 0.0, 3.0, 180.0, 180.0, Roughness::VerySmooth);
        assert!((total - 9.05409826063257).abs() < 1.0e-9, "{total}");
        let rough = doe2_exterior_hc(28.0, 20.0, 0.0, 3.0, 180.0, 180.0, Roughness::Rough);
        assert!((rough - 13.36494409525639).abs() < 1.0e-9, "{rough}");
    }

    #[test]
    fn windward_follows_the_energyplus_angle_test() {
        assert!(windward(0.0, 180.0, 180.0));
        assert!(!windward(0.0, 180.0, 0.0));
        assert!(windward(0.99, 180.0, 0.0));
    }

    #[test]
    fn south_facing_normal_is_azimuth_180() {
        let azimuth = azimuth_from_normal([0.0, -1.0, 0.0]);
        assert!((azimuth - 180.0).abs() < 1.0e-9, "{azimuth}");
    }
}
