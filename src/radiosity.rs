//! Interior long-wave exchange from EnergyPlus `HeatBalanceIntRadExchange`.
//!
//! View factors are `CalcApproximateViewFactors` plus the reciprocity and
//! closure loop from `FixViewFactors`. The exchange coefficients are Hottel's
//! Script F (`CalcScriptF`), multiplied by the EnergyPlus Stefan-Boltzmann
//! constant. Net flux on surface `i` is sum_j ScriptF_ij * (T_j^4 - T_i^4),
//! in W/m², positive when the surface gains heat.

use crate::KELVIN;

/// EnergyPlus `Constant::StefanBoltzmann`, W/(m²·K⁴).
pub const STEFAN_BOLTZMANN: f64 = 5.6697e-8;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurfaceClass {
    Wall,
    Floor,
    Roof,
    Ceiling,
}

#[derive(Clone, Copy, Debug)]
pub struct EnclosureSurface {
    pub area: f64,
    pub azimuth_deg: f64,
    pub tilt_deg: f64,
    pub emissivity: f64,
    pub class: SurfaceClass,
}

/// Tilt of an outward unit normal, degrees from up. Floors are near 180.
pub fn tilt_from_normal(normal_z: f64) -> f64 {
    normal_z.clamp(-1.0, 1.0).acos().to_degrees()
}

pub fn class_from_tilt(tilt_deg: f64) -> SurfaceClass {
    if tilt_deg > 170.0 {
        SurfaceClass::Floor
    } else if tilt_deg < 10.0 {
        SurfaceClass::Roof
    } else {
        SurfaceClass::Wall
    }
}

/// Net long-wave flux for one enclosure, W/m², aligned with `surfaces`.
pub fn net_longwave(surfaces: &[EnclosureSurface], temps_c: &[f64]) -> Vec<f64> {
    let n = surfaces.len();
    if n == 0 || temps_c.len() != n {
        return vec![0.0; temps_c.len()];
    }
    if n == 1 {
        return vec![0.0];
    }
    let Some(script_f) = script_f(surfaces) else {
        return vec![0.0; n];
    };
    let mut flux = vec![0.0; n];
    for i in 0..n {
        let ti = (temps_c[i] + KELVIN).powi(4);
        let mut sum = 0.0;
        for j in 0..n {
            if i == j {
                continue;
            }
            let tj = (temps_c[j] + KELVIN).powi(4);
            sum += script_f[i][j] * (tj - ti);
        }
        flux[i] = sum;
    }
    flux
}

fn script_f(surfaces: &[EnclosureSurface]) -> Option<Vec<Vec<f64>>> {
    let n = surfaces.len();
    let areas: Vec<f64> = surfaces.iter().map(|s| s.area).collect();
    let mut f = approximate_view_factors(surfaces);
    fix_view_factors(&areas, &mut f);
    let mut emiss: Vec<f64> = surfaces
        .iter()
        .map(|s| s.emissivity.clamp(0.0, 0.99999))
        .collect();
    for value in &mut emiss {
        if *value > 0.99999 {
            *value = 0.99999;
        }
    }

    let mut c = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            c[i][j] = areas[i] * f[i][j];
        }
    }
    let mut excite = vec![0.0; n];
    for i in 0..n {
        let reflect = (1.0 - emiss[i]).max(1.0e-8);
        let emiss_fac = areas[i] / reflect;
        excite[i] = -emiss[i] * emiss_fac;
        c[i][i] -= emiss_fac;
    }
    let mut inverse = invert(&c)?;
    for j in 0..n {
        for i in 0..n {
            inverse[i][j] *= excite[j];
        }
    }
    let mut script = vec![vec![0.0; n]; n];
    for i in 0..n {
        let reflect = (1.0 - emiss[i]).max(1.0e-8);
        let emiss_fac = emiss[i] / reflect;
        for j in 0..n {
            let delta = if i == j { emiss[i] } else { 0.0 };
            script[i][j] = emiss_fac * (inverse[i][j] - delta) * STEFAN_BOLTZMANN;
        }
    }
    Some(script)
}

fn sees(i: &EnclosureSurface, j: &EnclosureSurface) -> bool {
    const SAME_ANGLE: f64 = 10.0;
    if i.class == SurfaceClass::Floor && j.class == SurfaceClass::Floor {
        return false;
    }
    if i.class == SurfaceClass::Floor
        || j.class == SurfaceClass::Floor
        || i.class == SurfaceClass::Roof
        || j.class == SurfaceClass::Roof
        || i.class == SurfaceClass::Ceiling
        || j.class == SurfaceClass::Ceiling
    {
        // Floors are seen by every other class. Roofs and ceilings are seen
        // when the tilt test below also passes; floors always see them.
        if i.class == SurfaceClass::Floor || j.class == SurfaceClass::Floor {
            return true;
        }
    }
    let az = (i.azimuth_deg - j.azimuth_deg).abs();
    let tilt = (i.tilt_deg - j.tilt_deg).abs();
    (az > SAME_ANGLE && az < 360.0 - SAME_ANGLE) || tilt > SAME_ANGLE
}

fn approximate_view_factors(surfaces: &[EnclosureSurface]) -> Vec<Vec<f64>> {
    let n = surfaces.len();
    let mut seen = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            if i != j && sees(&surfaces[i], &surfaces[j]) {
                seen[i] += surfaces[j].area;
            }
        }
    }
    let mut f = vec![vec![0.0; n]; n];
    for i in 0..n {
        if seen[i] <= 0.0 {
            continue;
        }
        for j in 0..n {
            if i != j && sees(&surfaces[i], &surfaces[j]) {
                f[i][j] = surfaces[j].area / seen[i];
            }
        }
    }
    f
}

fn fix_view_factors(areas: &[f64], f: &mut [Vec<f64>]) {
    let n = areas.len();
    if n == 0 {
        return;
    }
    let total: f64 = areas.iter().sum();
    let largest = areas.iter().cloned().fold(0.0, f64::max);
    if n > 3 && largest > 0.99 * (total - largest) {
        if let Some(i) = areas.iter().position(|area| *area == largest) {
            f[i][i] = (1.2 * largest / total).min(0.9);
        }
    }
    let mut af = vec![vec![0.0; n]; n];
    for from in 0..n {
        for to in 0..n {
            af[from][to] = f[from][to] * areas[from];
        }
    }
    if n <= 3 {
        enforce_reciprocity(&mut af);
        for from in 0..n {
            for to in 0..n {
                f[from][to] = if areas[from] > 0.0 {
                    af[from][to] / areas[from]
                } else {
                    0.0
                };
            }
        }
        let mut max_row = 0.0_f64;
        for row in f.iter() {
            max_row = max_row.max(row.iter().sum());
        }
        if max_row > 1.01 {
            for row in f.iter_mut() {
                for value in row.iter_mut() {
                    *value /= max_row;
                }
            }
        }
        return;
    }
    let mut previous = 10.0;
    for _ in 0..400 {
        for from in 0..n {
            let sum: f64 = af[from].iter().sum();
            let coefficient = if sum.abs() > 1.0e-10 {
                areas[from] / sum
            } else {
                1.0
            };
            for value in &mut af[from] {
                *value *= coefficient;
            }
        }
        enforce_reciprocity(&mut af);
        for from in 0..n {
            for to in 0..n {
                let view = if areas[from] > 0.0 {
                    af[from][to] / areas[from]
                } else {
                    0.0
                };
                if view.abs() < 1.0e-10 {
                    f[from][to] = 0.0;
                    af[from][to] = 0.0;
                } else {
                    f[from][to] = view;
                }
            }
        }
        let sum_f: f64 = f.iter().map(|row| row.iter().sum::<f64>()).sum();
        let gap = (sum_f - n as f64).abs();
        if (previous - gap).abs() < 0.00001 || gap <= 0.001 {
            break;
        }
        previous = gap;
    }
}

fn enforce_reciprocity(af: &mut [Vec<f64>]) {
    let n = af.len();
    for i in 0..n {
        for j in i..n {
            let average = 0.5 * (af[i][j] + af[j][i]);
            af[i][j] = average;
            af[j][i] = average;
        }
    }
}

fn invert(matrix: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = matrix.len();
    let mut rows = vec![vec![0.0; n * 2]; n];
    for i in 0..n {
        for j in 0..n {
            rows[i][j] = matrix[i][j];
        }
        rows[i][n + i] = 1.0;
    }
    for col in 0..n {
        let mut pivot = col;
        for row in col + 1..n {
            if rows[row][col].abs() > rows[pivot][col].abs() {
                pivot = row;
            }
        }
        if rows[pivot][col].abs() < 1.0e-14 {
            return None;
        }
        rows.swap(col, pivot);
        let scale = rows[col][col];
        for value in &mut rows[col] {
            *value /= scale;
        }
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = rows[row][col];
            if factor == 0.0 {
                continue;
            }
            for j in 0..n * 2 {
                rows[row][j] -= factor * rows[col][j];
            }
        }
    }
    Some(rows.into_iter().map(|row| row[n..].to_vec()).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plate(
        area: f64,
        emissivity: f64,
        class: SurfaceClass,
        azimuth: f64,
        tilt: f64,
    ) -> EnclosureSurface {
        EnclosureSurface {
            area,
            azimuth_deg: azimuth,
            tilt_deg: tilt,
            emissivity,
            class,
        }
    }

    #[test]
    fn two_parallel_plates_match_the_grey_body_formula() {
        let surfaces = [
            plate(2.0, 0.9, SurfaceClass::Floor, 0.0, 180.0),
            plate(2.0, 0.9, SurfaceClass::Roof, 0.0, 0.0),
        ];
        let flux = net_longwave(&surfaces, &[30.0, 10.0]);
        let t1 = 30.0 + KELVIN;
        let t2 = 10.0 + KELVIN;
        let expected = STEFAN_BOLTZMANN * (t2.powi(4) - t1.powi(4)) / (1.0 / 0.9 + 1.0 / 0.9 - 1.0);
        assert!(
            (flux[0] - expected).abs() / expected.abs() < 1.0e-6,
            "flux {} expected {expected}",
            flux[0]
        );
        let net_power = flux[0] * 2.0 + flux[1] * 2.0;
        assert!(net_power.abs() < 1.0e-6, "enclosure power {net_power}");
    }

    #[test]
    fn six_sided_box_conserves_longwave() {
        let walls = [
            (0.0, 90.0, SurfaceClass::Wall),
            (90.0, 90.0, SurfaceClass::Wall),
            (180.0, 90.0, SurfaceClass::Wall),
            (270.0, 90.0, SurfaceClass::Wall),
            (0.0, 180.0, SurfaceClass::Floor),
            (0.0, 0.0, SurfaceClass::Roof),
        ];
        let surfaces: Vec<_> = walls
            .into_iter()
            .map(|(azimuth, tilt, class)| plate(10.0, 0.9, class, azimuth, tilt))
            .collect();
        let temps = [18.0, 19.0, 21.0, 22.0, 16.0, 24.0];
        let flux = net_longwave(&surfaces, &temps);
        let power: f64 = flux.iter().map(|q| q * 10.0).sum();
        assert!(power.abs() < 1.0e-4, "box power {power}");
    }
}
