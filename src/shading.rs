//! Shadows from a `Shading:Overhang` on the window it references.
//!
//! The overhang is a rectangle built in the window's local axes, matching the
//! EnergyPlus dimensions: height above the window, tilt from the window, left
//! and right extensions, and depth. Beam sunlit fraction is the window area the
//! overhang's shadow does not cover. Diffuse ratios follow `SkyDifSolarShading`:
//! the sky dome is split into equal-radiance patches, and each patch is shaded
//! with the same projection. Incoming rays travel opposite the sun direction.

/// One overhang attached to a window or door.
#[derive(Clone, Debug)]
pub struct Overhang {
    pub window: usize,
    pub vertices: Vec<[f64; 3]>,
}

/// Rectangle for an overhang on `window_vertices`.
///
/// `tilt_from_window_deg` of 90 puts the plate perpendicular to the window.
/// Local X is `cross(up, outward normal)` and left is the negative-X side.
pub fn overhang_rectangle(
    window_vertices: &[[f64; 3]],
    normal: [f64; 3],
    height_above: f64,
    tilt_from_window_deg: f64,
    left: f64,
    right: f64,
    depth: f64,
) -> Vec<[f64; 3]> {
    if window_vertices.is_empty() || depth <= 0.0 {
        return Vec::new();
    }
    let origin = window_vertices[0];
    let (axis_u, axis_v) = window_axes(normal);
    let locals: Vec<(f64, f64)> = window_vertices
        .iter()
        .map(|point| local(*point, origin, axis_u, axis_v))
        .collect();
    let min_u = locals.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
    let max_u = locals.iter().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max);
    let max_v = locals.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);
    let tilt = tilt_from_window_deg.to_radians();
    let depth_dir = normalize(add(scale(normal, tilt.sin()), scale(axis_v, tilt.cos())));
    let top = max_v + height_above;
    let inner_left = add(add(origin, scale(axis_u, min_u - left)), scale(axis_v, top));
    let inner_right = add(
        add(origin, scale(axis_u, max_u + right)),
        scale(axis_v, top),
    );
    let outer_right = add(inner_right, scale(depth_dir, depth));
    let outer_left = add(inner_left, scale(depth_dir, depth));
    vec![inner_left, inner_right, outer_right, outer_left]
}

/// Fraction of the window that still receives beam sun. `sun` points toward the sun.
pub fn beam_sunlit_fraction(
    window: &[[f64; 3]],
    normal: [f64; 3],
    shade: &[[f64; 3]],
    sun: [f64; 3],
) -> f64 {
    if window.len() < 3 || shade.len() < 3 {
        return 1.0;
    }
    let cos_inc = dot(normal, sun);
    if cos_inc <= 1.0e-8 {
        return 1.0;
    }
    let origin = window[0];
    let projected: Vec<[f64; 3]> = shade
        .iter()
        .filter_map(|point| project_to_plane(*point, origin, normal, sun, cos_inc))
        .collect();
    if projected.len() < 3 {
        return 1.0;
    }
    let (axis_u, axis_v) = window_axes(normal);
    let window_2d: Vec<(f64, f64)> = window
        .iter()
        .map(|point| local(*point, origin, axis_u, axis_v))
        .collect();
    let shade_2d: Vec<(f64, f64)> = projected
        .iter()
        .map(|point| local(*point, origin, axis_u, axis_v))
        .collect();
    let window_area = polygon_area_2d(&window_2d).abs();
    if window_area < 1.0e-8 {
        return 1.0;
    }
    let overlap = clip_polygon(&shade_2d, &window_2d);
    let shaded = polygon_area_2d(&overlap).abs();
    (1.0 - (shaded / window_area).clamp(0.0, 1.0)).clamp(0.0, 1.0)
}

/// Shaded sky diffuse divided by the unshaded value.
///
/// `isotropic` multiplies sky diffuse. `horizon` is the same ratio for the
/// lowest altitude ring, kept for a later anisotropic sky.
#[derive(Clone, Copy, Debug)]
pub struct DiffuseShadeRatio {
    pub isotropic: f64,
    pub horizon: f64,
}

/// Isotropic and horizon-band shade ratios for every overhang on one surface.
///
/// Patch altitudes are 7.5° through 82.5° and azimuths are 0° through 345°,
/// the `NPhi = 6`, `NTheta = 24` grid in `SolarShading.cc`. A patch behind the
/// surface is left out of both sums. With no overhang the ratios are 1.
pub fn diffuse_shade_ratios(
    window: &[[f64; 3]],
    normal: [f64; 3],
    shades: &[&[[f64; 3]]],
) -> DiffuseShadeRatio {
    if window.len() < 3 || shades.is_empty() {
        return DiffuseShadeRatio {
            isotropic: 1.0,
            horizon: 1.0,
        };
    }
    const N_PHI: usize = 6;
    const N_THETA: usize = 24;
    const EPS: f64 = 1.0e-10;
    let d_phi = std::f64::consts::FRAC_PI_2 / N_PHI as f64;
    let d_theta = std::f64::consts::TAU / N_THETA as f64;
    let weight = d_theta * d_phi;
    let phi_min = 0.5 * d_phi;
    let mut with_iso = 0.0;
    let mut wo_iso = 0.0;
    let mut with_horiz = 0.0;
    let mut wo_horiz = 0.0;
    for i_phi in 0..N_PHI {
        let phi = phi_min + i_phi as f64 * d_phi;
        let (sin_phi, cos_phi) = phi.sin_cos();
        for i_theta in 0..N_THETA {
            let theta = i_theta as f64 * d_theta;
            let (sin_theta, cos_theta) = theta.sin_cos();
            let sun = [cos_phi * cos_theta, cos_phi * sin_theta, sin_phi];
            let cos_inc = dot(normal, sun);
            if cos_inc < 0.0 {
                continue;
            }
            let fac = cos_phi * weight * cos_inc;
            let mut fraction = 1.0_f64;
            for shade in shades {
                fraction = fraction.min(beam_sunlit_fraction(window, normal, shade, sun));
            }
            let shaded = fac * fraction;
            with_iso += shaded;
            wo_iso += fac;
            if i_phi == 0 {
                with_horiz += shaded;
                wo_horiz += fac;
            }
        }
    }
    DiffuseShadeRatio {
        isotropic: shade_ratio(with_iso, wo_iso, EPS),
        horizon: shade_ratio(with_horiz, wo_horiz, EPS),
    }
}

fn shade_ratio(with_shade: f64, without_shade: f64, eps: f64) -> f64 {
    if without_shade.abs() > eps {
        with_shade / without_shade
    } else {
        with_shade / (without_shade + eps)
    }
}

fn project_to_plane(
    point: [f64; 3],
    origin: [f64; 3],
    normal: [f64; 3],
    sun: [f64; 3],
    cos_inc: f64,
) -> Option<[f64; 3]> {
    let t = dot(sub(point, origin), normal) / cos_inc;
    if t < -1.0e-8 {
        return None;
    }
    Some(sub(point, scale(sun, t)))
}

fn window_axes(normal: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let mut axis_u = cross([0.0, 0.0, 1.0], normal);
    if norm(axis_u) < 1.0e-8 {
        axis_u = [1.0, 0.0, 0.0];
    }
    let axis_u = normalize(axis_u);
    let axis_v = normalize(cross(normal, axis_u));
    (axis_u, axis_v)
}

fn local(point: [f64; 3], origin: [f64; 3], axis_u: [f64; 3], axis_v: [f64; 3]) -> (f64, f64) {
    let delta = sub(point, origin);
    (dot(delta, axis_u), dot(delta, axis_v))
}

fn clip_polygon(subject: &[(f64, f64)], clip: &[(f64, f64)]) -> Vec<(f64, f64)> {
    if subject.is_empty() || clip.len() < 3 {
        return Vec::new();
    }
    let centroid = polygon_centroid(clip);
    let mut output = subject.to_vec();
    for i in 0..clip.len() {
        if output.is_empty() {
            break;
        }
        let a = clip[i];
        let b = clip[(i + 1) % clip.len()];
        let input = output;
        output = Vec::new();
        let mut previous = *input.last().unwrap();
        for &current in &input {
            let current_in = inside_half_plane(current, a, b, centroid);
            let previous_in = inside_half_plane(previous, a, b, centroid);
            if current_in {
                if !previous_in {
                    if let Some(hit) = segment_hit(previous, current, a, b) {
                        output.push(hit);
                    }
                }
                output.push(current);
            } else if previous_in {
                if let Some(hit) = segment_hit(previous, current, a, b) {
                    output.push(hit);
                }
            }
            previous = current;
        }
    }
    output
}

fn inside_half_plane(
    point: (f64, f64),
    a: (f64, f64),
    b: (f64, f64),
    centroid: (f64, f64),
) -> bool {
    let edge = (b.0 - a.0, b.1 - a.1);
    let to_point = (point.0 - a.0, point.1 - a.1);
    let to_center = (centroid.0 - a.0, centroid.1 - a.1);
    let point_side = edge.0 * to_point.1 - edge.1 * to_point.0;
    let center_side = edge.0 * to_center.1 - edge.1 * to_center.0;
    point_side * center_side >= -1.0e-10
}

fn segment_hit(p: (f64, f64), q: (f64, f64), a: (f64, f64), b: (f64, f64)) -> Option<(f64, f64)> {
    let r = (q.0 - p.0, q.1 - p.1);
    let s = (b.0 - a.0, b.1 - a.1);
    let denom = r.0 * s.1 - r.1 * s.0;
    if denom.abs() < 1.0e-12 {
        return None;
    }
    let qp = (a.0 - p.0, a.1 - p.1);
    let t = (qp.0 * s.1 - qp.1 * s.0) / denom;
    Some((p.0 + t * r.0, p.1 + t * r.1))
}

fn polygon_centroid(points: &[(f64, f64)]) -> (f64, f64) {
    let n = points.len() as f64;
    let x = points.iter().map(|p| p.0).sum::<f64>() / n;
    let y = points.iter().map(|p| p.1).sum::<f64>() / n;
    (x, y)
}

fn polygon_area_2d(points: &[(f64, f64)]) -> f64 {
    if points.len() < 3 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 0..points.len() {
        let j = (i + 1) % points.len();
        sum += points[i].0 * points[j].1 - points[j].0 * points[i].1;
    }
    0.5 * sum
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: [f64; 3], factor: f64) -> [f64; 3] {
    [a[0] * factor, a[1] * factor, a[2] * factor]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn normalize(a: [f64; 3]) -> [f64; 3] {
    let length = norm(a);
    if length < 1.0e-12 {
        a
    } else {
        scale(a, 1.0 / length)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn forty_five_degree_sun_shades_half_of_a_two_meter_window() {
        let window = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 0.0, 2.0],
            [0.0, 0.0, 2.0],
        ];
        let normal = [0.0, -1.0, 0.0];
        let shade = overhang_rectangle(&window, normal, 0.0, 90.0, 0.0, 0.0, 1.0);
        let s = (0.5_f64).sqrt();
        let sun = [0.0, -s, s];
        let fraction = beam_sunlit_fraction(&window, normal, &shade, sun);
        assert!(
            (fraction - 0.5).abs() < 0.02,
            "sunlit {fraction}, shade {shade:?}"
        );
    }

    #[test]
    fn sun_behind_the_wall_does_not_cast_an_overhang_shadow() {
        let window = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 0.0, 2.0],
            [0.0, 0.0, 2.0],
        ];
        let normal = [0.0, -1.0, 0.0];
        let shade = overhang_rectangle(&window, normal, 0.0, 90.0, 0.0, 0.0, 1.0);
        let fraction = beam_sunlit_fraction(&window, normal, &shade, [0.0, 1.0, 0.5]);
        assert!((fraction - 1.0).abs() < 1.0e-9);
    }

    #[test]
    fn a_window_with_no_overhang_keeps_the_full_sky() {
        let window = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 0.0, 2.0],
            [0.0, 0.0, 2.0],
        ];
        let ratio = diffuse_shade_ratios(&window, [0.0, -1.0, 0.0], &[]);
        assert!((ratio.isotropic - 1.0).abs() < 1.0e-12);
        assert!((ratio.horizon - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn a_horizontal_overhang_reduces_isotropic_sky_diffuse() {
        let window = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 0.0, 2.0],
            [0.0, 0.0, 2.0],
        ];
        let normal = [0.0, -1.0, 0.0];
        let shade = overhang_rectangle(&window, normal, 0.0, 90.0, 0.0, 0.0, 1.0);
        let ratio = diffuse_shade_ratios(&window, normal, &[&shade]);
        assert!(
            (ratio.isotropic - 0.706).abs() < 0.02,
            "isotropic {}",
            ratio.isotropic
        );
        assert!(
            (ratio.horizon - 0.933).abs() < 0.02,
            "horizon {}",
            ratio.horizon
        );
    }
}
