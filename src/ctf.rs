//! Conduction transfer functions from `Construction::calculateTransferFunction`.
//!
//! Coefficients use Seem's state-space method: a finite-difference grid, the
//! exponential of `A * dt`, and the flux history series. EnergyPlus converts
//! layer properties to English units before that algebra and converts the
//! finished coefficients back to SI. This module follows that path for a
//! one-dimensional wall with no internal source.

/// `HVAC` surface temperature iteration is elsewhere. CTF history uses hours.
const PHYS_PROP_LIMIT: f64 = 1.0e-6;
const MIN_NODES: usize = 6;
const MAX_CTF_TERMS: usize = 19;
const MAX_SUM_ERROR: f64 = 0.01;
const MAX_TIME_STEP_HOURS: f64 = 7.0;
const SERIES_LIMIT: f64 = 1.0e-13;
const DP_LIMIT: f64 = 1.0e-20;

const CFC: f64 = 4.1868;
const CFL: f64 = 0.3048;
const CFM: f64 = 0.45359237;
const CFT: f64 = 5.0 / 9.0;
const CFE: f64 = CFC * CFM * CFT / 3.6;
const CFD: f64 = CFM / (CFL * CFL * CFL);
const CFK: f64 = CFE / (CFL * CFT);
const CFU: f64 = CFK / CFL;

#[derive(Clone, Debug)]
pub struct CtfLayer {
    pub thickness: f64,
    pub conductivity: f64,
    pub density: f64,
    pub specific_heat: f64,
}

/// SI conduction transfer functions. Index 0 is the current timestep.
/// `flux[0]` is unused. Flux history is dimensionless.
#[derive(Clone, Debug)]
pub struct Ctf {
    pub outside: Vec<f64>,
    pub cross: Vec<f64>,
    pub inside: Vec<f64>,
    pub flux: Vec<f64>,
    pub u_value: f64,
    pub num_terms: usize,
    /// Zone timesteps between history updates. `1` when the CTF step is the zone step.
    pub zone_steps_per_update: usize,
}

/// Build CTFs for `layers`, outside layer first. `zone_timestep_hours` is `TimeStepZone`.
pub fn conduction_transfer_functions(
    layers_in: &[CtfLayer],
    zone_timestep_hours: f64,
) -> Result<Ctf, String> {
    if layers_in.is_empty() {
        return Err("construction has no layers".to_string());
    }
    if zone_timestep_hours <= 0.0 {
        return Err("zone timestep must be positive".to_string());
    }

    let mut thickness = Vec::new();
    let mut conductivity = Vec::new();
    let mut density = Vec::new();
    let mut specific_heat = Vec::new();
    let mut resistance = Vec::new();
    let mut resistive = Vec::new();

    for layer in layers_in {
        if layer.thickness > 3.0 {
            return Err(format!(
                "material thickness {:.4} m is above the 3 m CTF limit",
                layer.thickness
            ));
        }
        let mut thick = layer.thickness;
        let mut k = layer.conductivity;
        let mut rho = layer.density;
        let mut cp = layer.specific_heat;
        let mut r = if k > PHYS_PROP_LIMIT { thick / k } else { 0.0 };
        let is_resistive = if k <= PHYS_PROP_LIMIT {
            true
        } else {
            thick * (rho * cp / k).sqrt() < PHYS_PROP_LIMIT
        };
        if is_resistive {
            if r < 1.0e-3 {
                return Err(
                    "resistive layer R-value is below the 0.001 m²K/W CTF limit".to_string()
                );
            }
            // Outer resistive layers become an equivalent still-air film.
            cp = 1.007e3;
            rho = 1.1614;
            k = 0.0263;
            thick = k * r;
            r = thick / k;
        }
        thickness.push(thick);
        conductivity.push(k);
        density.push(rho);
        specific_heat.push(cp);
        resistance.push(r);
        resistive.push(is_resistive);
    }

    for i in 0..thickness.len() {
        resistance[i] *= CFU;
        thickness[i] /= CFL;
        conductivity[i] /= CFK;
        density[i] /= CFD;
        specific_heat[i] /= CFC * 1000.0;
    }

    let r_total: f64 = resistance.iter().sum();
    if r_total <= 0.0 {
        return Err("construction resistance is zero".to_string());
    }
    let cnd = 1.0 / r_total;
    let n_layers = thickness.len();
    let any_massive = resistive.iter().any(|flag| !flag);
    if !any_massive {
        return Ok(resistive_ctf(cnd));
    }

    let mut nodes = vec![0usize; n_layers];
    let mut dx = vec![0.0; n_layers];
    let mut ctf_step = 0.0;
    for layer in 0..n_layers {
        if resistive[layer] && layer > 0 && layer + 1 < n_layers {
            nodes[layer] = 1;
            dx[layer] = thickness[layer];
            continue;
        }
        let alpha = conductivity[layer] / density[layer] / specific_heat[layer];
        let spacing = (2.0 * alpha * zone_timestep_hours).sqrt();
        let estimate = (thickness[layer] / spacing) as usize;
        nodes[layer] = estimate.clamp(MIN_NODES, MAX_CTF_TERMS);
        dx[layer] = thickness[layer] / nodes[layer] as f64;
    }

    let mut rcmax = nodes.iter().sum::<usize>();
    if rcmax == 0 {
        return Err("construction produced no CTF nodes".to_string());
    }
    rcmax -= 1;

    for layer in 0..n_layers {
        if nodes[layer] >= MAX_CTF_TERMS {
            let step =
                density[layer] * specific_heat[layer] * dx[layer] * dx[layer] / conductivity[layer];
            if step > ctf_step {
                ctf_step = step;
            }
        }
    }
    let mut histories = 1usize;
    if ((zone_timestep_hours - ctf_step) / zone_timestep_hours).abs() > 0.1 {
        if ctf_step > zone_timestep_hours {
            histories = (ctf_step / zone_timestep_hours).round().max(1.0) as usize;
            ctf_step = zone_timestep_hours * histories as f64;
        } else {
            ctf_step = zone_timestep_hours;
            histories = 1;
        }
    }

    let mut outside = vec![0.0; MAX_CTF_TERMS];
    let mut cross = vec![0.0; MAX_CTF_TERMS];
    let mut inside = vec![0.0; MAX_CTF_TERMS];
    let mut flux = vec![0.0; MAX_CTF_TERMS];

    let num_terms = loop {
        let (s0, s_hist, e, terms) = calculate_coefficients(
            &conductivity,
            &density,
            &specific_heat,
            &dx,
            &nodes,
            rcmax,
            ctf_step,
        )?;
        if terms > MAX_CTF_TERMS - 1 {
            histories += 1;
            ctf_step += zone_timestep_hours;
            if ctf_step >= MAX_TIME_STEP_HOURS {
                return Err(
                    "CTF series did not converge before the 7 hour timestep limit".to_string(),
                );
            }
            continue;
        }
        // Sum(X) = Sum(Y) = Sum(Z) in the English coefficients, before the SI conversion.
        let mut sum_x = s0[1][1];
        let mut sum_y = s0[0][1];
        let mut sum_z = s0[0][0];
        for term in 1..=terms {
            sum_x += s_hist[term][1][1];
            sum_y += s_hist[term][0][1];
            sum_z += s_hist[term][0][0];
        }
        sum_x = sum_x.abs();
        sum_y = sum_y.abs();
        sum_z = sum_z.abs();
        let biggest = sum_x.max(sum_y).max(sum_z);
        if biggest == 0.0 {
            return Err("CTF series summed to zero".to_string());
        }
        if (sum_x - sum_y).abs() / biggest > MAX_SUM_ERROR
            || (sum_z - sum_y).abs() / biggest > MAX_SUM_ERROR
        {
            histories += 1;
            ctf_step += zone_timestep_hours;
            if ctf_step >= MAX_TIME_STEP_HOURS {
                return Err(
                    "CTF series did not converge before the 7 hour timestep limit".to_string(),
                );
            }
            continue;
        }

        outside[0] = s0[0][0] * CFU;
        cross[0] = s0[0][1] * CFU;
        inside[0] = -s0[1][1] * CFU;
        for term in 1..=terms {
            outside[term] = s_hist[term][0][0] * CFU;
            cross[term] = s_hist[term][0][1] * CFU;
            inside[term] = -s_hist[term][1][1] * CFU;
            flux[term] = -e[term];
        }
        break terms;
    };

    let keep = num_terms + 1;
    Ok(Ctf {
        outside: outside[..keep].to_vec(),
        cross: cross[..keep].to_vec(),
        inside: inside[..keep].to_vec(),
        flux: flux[..keep].to_vec(),
        u_value: cnd * CFU,
        num_terms,
        zone_steps_per_update: histories.max(1),
    })
}

fn resistive_ctf(cnd_english: f64) -> Ctf {
    let u = cnd_english * CFU;
    Ctf {
        outside: vec![u, 0.0],
        cross: vec![u, 0.0],
        inside: vec![u, 0.0],
        flux: vec![0.0, 0.0],
        u_value: u,
        num_terms: 1,
        zone_steps_per_update: 1,
    }
}

fn calculate_coefficients(
    conductivity: &[f64],
    density: &[f64],
    specific_heat: &[f64],
    dx: &[f64],
    nodes: &[usize],
    rcmax: usize,
    ctf_step: f64,
) -> Result<([[f64; 2]; 2], Vec<[[f64; 2]; 2]>, Vec<f64>, usize), String> {
    let n_layers = conductivity.len();
    let mut amat = zeros(rcmax);
    let mut bmat = [0.0; 2];
    let cap = density[0] * specific_heat[0] * dx[0] * 1.5;
    let dxtmp = 1.0 / dx[0] / cap;
    amat[0][0] = -2.0 * conductivity[0] * dxtmp;
    if rcmax > 1 {
        amat[1][0] = conductivity[0] * dxtmp;
    }
    bmat[0] = conductivity[0] * dxtmp;

    let mut layer = 0usize;
    let mut node_in_layer = 2usize;
    for node in 1..rcmax.saturating_sub(1) {
        if node_in_layer == nodes[layer] && n_layers != 1 {
            let next = layer + 1;
            let cap = 0.5
                * (density[layer] * specific_heat[layer] * dx[layer]
                    + density[next] * specific_heat[next] * dx[next]);
            if cap.abs() < 1.0e-12 {
                return Err("CTF interface node has no thermal mass".to_string());
            }
            amat[node - 1][node] = conductivity[layer] / dx[layer] / cap;
            amat[node][node] =
                -(conductivity[layer] / dx[layer] + conductivity[next] / dx[next]) / cap;
            amat[node + 1][node] = conductivity[next] / dx[next] / cap;
            node_in_layer = 0;
            layer += 1;
        } else {
            let cap = density[layer] * specific_heat[layer] * dx[layer];
            let dxtmp = 1.0 / dx[layer] / cap;
            amat[node - 1][node] = conductivity[layer] * dxtmp;
            amat[node][node] = -2.0 * conductivity[layer] * dxtmp;
            amat[node + 1][node] = conductivity[layer] * dxtmp;
        }
        node_in_layer += 1;
    }

    let last = n_layers - 1;
    let cap = density[last] * specific_heat[last] * dx[last] * 1.5;
    let dxtmp = 1.0 / dx[last] / cap;
    amat[rcmax - 1][rcmax - 1] = -2.0 * conductivity[last] * dxtmp;
    if rcmax > 1 {
        amat[rcmax - 2][rcmax - 1] = conductivity[last] * dxtmp;
    }
    bmat[1] = conductivity[last] * dxtmp;
    let cmat = [-conductivity[0] / dx[0], conductivity[last] / dx[last]];
    let dmat = [conductivity[0] / dx[0], -conductivity[last] / dx[last]];

    let aexp = matrix_exponential(&amat, ctf_step);
    let ainv = invert_matrix(&amat)?;
    let (gamma1, gamma2) = gammas(&ainv, &aexp, &bmat, ctf_step);

    let mut s0 = [[0.0; 2]; 2];
    s0[0][0] = cmat[0] * gamma2[0][0] + dmat[0];
    s0[1][0] = cmat[0] * gamma2[1][0];
    s0[0][1] = cmat[1] * gamma2[0][rcmax - 1];
    s0[1][1] = cmat[1] * gamma2[1][rcmax - 1] + dmat[1];
    let (left, right) = (s0[1][0], s0[0][1]);
    let (left, right) = symmetrized(left, right);
    s0[1][0] = left;
    s0[0][1] = right;

    let mut e = vec![0.0; rcmax + 1];
    let mut s_hist = vec![[[0.0; 2]; 2]; rcmax + 1];
    let mut rnew = identity(rcmax);
    let mut num_terms = rcmax;
    let mut converged = false;
    let mut inum = 1usize;
    while !converged && inum < rcmax {
        let (phi_r0, trace) = phi_r(&aexp, &rnew);
        e[inum] = -trace / inum as f64;
        let rold = rnew.clone();
        rnew = phi_r0;
        for i in 0..rcmax {
            rnew[i][i] += e[inum];
        }
        for j in 0..2 {
            for state in 0..rcmax {
                s_hist[inum][j][0] += cmat[0]
                    * (rold[state][0] * gamma1[j][state] + rnew[state][0] * gamma2[j][state]);
                s_hist[inum][j][1] += cmat[1]
                    * (rold[state][rcmax - 1] * gamma1[j][state]
                        + rnew[state][rcmax - 1] * gamma2[j][state]);
            }
            s_hist[inum][j][j] += e[inum] * dmat[j];
        }
        let (left, right) = symmetrized(s_hist[inum][1][0], s_hist[inum][0][1]);
        s_hist[inum][1][0] = left;
        s_hist[inum][0][1] = right;
        if e[1] == 0.0 {
            num_terms = 1;
            converged = true;
        } else if (e[inum] / e[1]).abs() < SERIES_LIMIT {
            num_terms = inum;
            converged = true;
        }
        inum += 1;
    }
    if !converged {
        let mut trace = 0.0;
        for ir in 0..rcmax {
            for is in 0..rcmax {
                trace += aexp[is][ir] * rnew[ir][is];
            }
        }
        e[rcmax] = -trace / rcmax as f64;
        for j in 0..2 {
            for state in 0..rcmax {
                s_hist[rcmax][j][0] += cmat[0] * rnew[state][0] * gamma1[j][state];
                s_hist[rcmax][j][1] += cmat[1] * rnew[state][rcmax - 1] * gamma1[j][state];
            }
            s_hist[rcmax][j][j] += e[rcmax] * dmat[j];
        }
        let (left, right) = symmetrized(s_hist[rcmax][1][0], s_hist[rcmax][0][1]);
        s_hist[rcmax][1][0] = left;
        s_hist[rcmax][0][1] = right;
        num_terms = rcmax;
    }
    Ok((s0, s_hist, e, num_terms))
}

fn symmetrized(mut left: f64, mut right: f64) -> (f64, f64) {
    if left.abs() != right.abs() && left.abs() > 0.0 && right.abs() > 0.0 {
        let average = 0.5 * (left.abs() + right.abs());
        left *= average / left.abs();
        right *= average / right.abs();
    }
    (left, right)
}

fn gammas(
    ainv: &[Vec<f64>],
    aexp: &[Vec<f64>],
    bmat: &[f64; 2],
    step: f64,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = ainv.len();
    let mut atemp = zeros(n);
    for ir in 0..n {
        for ic in 0..n {
            let eye = if ic == ir { 1.0 } else { 0.0 };
            atemp[ic][ir] = aexp[ic][ir] - eye;
        }
    }
    let mut gamma1 = vec![vec![0.0; n]; 2];
    for i in 0..n {
        for state in 0..n {
            gamma1[0][i] += ainv[state][i] * atemp[0][state] * bmat[0];
            gamma1[1][i] += ainv[state][i] * atemp[n - 1][state] * bmat[1];
        }
    }
    let mut gamma2 = vec![vec![0.0; n]; 2];
    for i in 0..n {
        for j in 0..2 {
            for state in 0..n {
                let mut value = ainv[state][i] * (gamma1[j][state] / step);
                if (j == 0 && state == 0) || (j == 1 && state + 1 == n) {
                    value -= ainv[state][i] * bmat[j];
                }
                gamma2[j][i] += value;
            }
        }
    }
    // Seem stores Gamma1 as Gamma1 - Gamma2 before the history coefficients.
    for i in 0..n {
        for j in 0..2 {
            gamma1[j][i] -= gamma2[j][i];
        }
    }
    (gamma1, gamma2)
}

fn phi_r(aexp: &[Vec<f64>], rnew: &[Vec<f64>]) -> (Vec<Vec<f64>>, f64) {
    let n = aexp.len();
    let tiny = f64::EPSILON;
    let mut phi = zeros(n);
    let mut trace = 0.0;
    for ir in 0..n {
        for ic in 0..n {
            for state in 0..n {
                if rnew[ic][state].abs() > tiny
                    && aexp[state][ir].abs() > (tiny / rnew[ic][state]).abs()
                {
                    phi[ic][ir] += aexp[state][ir] * rnew[ic][state];
                }
            }
        }
        trace += phi[ir][ir];
    }
    (phi, trace)
}

fn matrix_exponential(amat: &[Vec<f64>], step: f64) -> Vec<Vec<f64>> {
    let n = amat.len();
    let tiny = f64::EPSILON;
    let mut scaled = amat.to_vec();
    let mut row_norm = 1.0e-12_f64;
    for column in 0..n {
        let mut sum = 0.0;
        for row in 0..n {
            sum += scaled[row][column].abs();
        }
        row_norm = row_norm.max(sum * step);
    }
    let k = (row_norm.ln() / 2.0_f64.ln()) as i32 + 1;
    let factor = step / 2.0_f64.powi(k);
    for row in &mut scaled {
        for value in row {
            *value *= factor;
        }
    }
    let limit = (3.0 * row_norm + 6.0).min(100.0) as i32;
    let mut aexp = identity(n);
    for i in 0..n {
        for j in 0..n {
            aexp[j][i] += scaled[j][i];
        }
    }
    let mut previous = scaled.clone();
    let mut power = 1i32;
    while power < limit {
        power += 1;
        let mut next = zeros(n);
        let mut significant = true;
        for ir in 0..n {
            for ic in 0..n {
                for mid in 0..n {
                    if scaled[ic][mid].abs() > tiny
                        && previous[mid][ir].abs() > (power as f64 * tiny / scaled[ic][mid]).abs()
                    {
                        next[ic][ir] += previous[mid][ir] * scaled[ic][mid] / power as f64;
                    }
                }
            }
        }
        previous = next;
        for ir in 0..n {
            for ic in 0..n {
                aexp[ic][ir] += previous[ic][ir];
            }
        }
        for ir in 0..n {
            for ic in 0..n {
                if aexp[ic][ir].abs() > tiny {
                    if (previous[ic][ir] / aexp[ic][ir]).abs() > DP_LIMIT {
                        significant = false;
                    }
                } else {
                    significant = false;
                }
            }
        }
        if significant {
            break;
        }
    }
    let mut squares = k.max(0);
    while squares > 0 {
        let old = aexp.clone();
        aexp = zeros(n);
        let mut backed_up = true;
        for ir in 0..n {
            for ic in 0..n {
                for mid in 0..n {
                    if (old[mid][ir] * old[ic][mid]).abs() > tiny {
                        aexp[ic][ir] += old[mid][ir] * old[ic][mid];
                        backed_up = false;
                    }
                }
            }
        }
        if backed_up {
            aexp = old;
            break;
        }
        squares -= 1;
    }
    aexp
}

fn invert_matrix(amat: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = amat.len();
    let mut work = amat.to_vec();
    let mut inverse = identity(n);
    for pivot in 0..n.saturating_sub(1) {
        let diag = work[pivot][pivot];
        if diag.abs() < 1.0e-18 {
            return Err("CTF state matrix is singular".to_string());
        }
        for col in pivot + 1..n {
            work[col][pivot] /= diag;
        }
        for col in 0..=pivot {
            inverse[col][pivot] /= diag;
        }
        work[pivot][pivot] = 1.0;
        for row in pivot + 1..n {
            let factor = work[pivot][row];
            for col in pivot + 1..n {
                work[col][row] -= factor * work[col][pivot];
            }
            for col in 0..=pivot {
                inverse[col][row] -= factor * inverse[col][pivot];
            }
            work[pivot][row] = 0.0;
        }
    }
    let diag = work[n - 1][n - 1];
    if diag.abs() < 1.0e-18 {
        return Err("CTF state matrix is singular".to_string());
    }
    for col in 0..n {
        inverse[col][n - 1] /= diag;
    }
    work[n - 1][n - 1] = 1.0;
    for pivot in (1..n).rev() {
        for row in 0..pivot {
            let factor = work[pivot][row];
            for col in 0..n {
                inverse[col][row] -= factor * inverse[col][pivot];
            }
            work[pivot][row] = 0.0;
        }
    }
    Ok(inverse)
}

fn zeros(n: usize) -> Vec<Vec<f64>> {
    vec![vec![0.0; n]; n]
}

fn identity(n: usize) -> Vec<Vec<f64>> {
    let mut matrix = zeros(n);
    for i in 0..n {
        matrix[i][i] = 1.0;
    }
    matrix
}

#[cfg(test)]
mod tests {
    use super::*;

    fn layer(thickness: f64, k: f64) -> CtfLayer {
        CtfLayer {
            thickness,
            conductivity: k,
            density: 2000.0,
            specific_heat: 1000.0,
        }
    }

    fn sums_match_u(ctf: &Ctf, u: f64) {
        // Steady state: q = U (T_out - T_in), with q = (Ysum T_out - Zsum T_in) / (1 - Fsum).
        let flux_sum: f64 = ctf.flux.iter().sum();
        let scale = (1.0 - flux_sum).abs().max(1.0e-12);
        let cross: f64 = ctf.cross.iter().sum::<f64>() / scale;
        let outside: f64 = ctf.outside.iter().sum::<f64>() / scale;
        let inside: f64 = ctf.inside.iter().sum::<f64>() / scale;
        assert!((cross - u).abs() / u < 0.02, "cross {cross} vs U {u}");
        assert!((outside - u).abs() / u < 0.02, "outside {outside}");
        assert!((inside - u).abs() / u < 0.02, "inside {inside}");
    }

    #[test]
    fn single_layer_series_sums_to_conductance() {
        let ctf = conduction_transfer_functions(&[layer(0.1, 0.1)], 1.0 / 6.0).unwrap();
        sums_match_u(&ctf, 1.0);
    }

    #[test]
    fn two_layers_sum_to_the_series_u_value() {
        let ctf = conduction_transfer_functions(&[layer(0.1, 0.04), layer(0.15, 1.4)], 1.0 / 6.0)
            .unwrap();
        let u = 1.0 / (0.1 / 0.04 + 0.15 / 1.4);
        sums_match_u(&ctf, u);
    }
}
