use crate::orbital::Orbital;

use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

const ZERO_PHI_PLANE: f64 = 0.0;
const R_BOUND_MAX: f64 = 40e-10; // 40 angstrom
const PROB_THRESHOLD: f64 = 0.08;

pub fn sample(n: u64, l: u64, m: i64, grid_size: usize, sample_amount: u64) -> Vec<Vec<bool>> {
    let orbital = Orbital::new(n, l, m);
    let grid_isize = grid_size as isize;

    let r_bound = discover_r_bound(grid_isize, &orbital, sample_amount).unwrap();

    // each pixel represents a length in space. we take the volume as a cube with side length of that
    // this volume is multiplied by |psi(r, theta, phi)|^2 at each point in space to get the probability
    // at the volume around that point
    let delta_volume = (r_bound / (grid_size as f64 / 2.0 - 1.0)).powi(3);

    // x and z are currently in range of [-grid_size/2, grid_size/2], we normalise that to [-r_bound, r_bound]
    // with the normalisation factor
    let norm_factor = r_bound as f64 / (grid_size as f64 / 2.0 - 1.0);

    let mut grid = vec![vec![false; grid_size]; grid_size];
    grid.par_iter_mut().enumerate().for_each(|(i, row)| {
        let mut rng = SmallRng::from_entropy();
        let z = ((grid_isize - 1) / 2 - i as isize) as f64 * norm_factor;
        row.iter_mut().enumerate().for_each(|(j, hit)| {
            let x = (j as isize - (grid_isize - 1) / 2) as f64 * norm_factor;

            let (r, theta, _) = spherical(x, 0.0, z);

            let mut p = orbital.probability(r, theta, ZERO_PHI_PLANE, delta_volume);

            // This calculates the probability at this point after sample_amount of sampling,
            // so later on in sample() we can sample each pixel only once, rather than
            // actually sampling sample_amount of them, which is very large.
            p = 1.0 - (1.0 - p).powf(sample_amount as f64);

            *hit |= rng.gen_bool(p);
        });
    });
    grid
}

// Find the maximum r_bound value so that we can see the entire orbital shape at a reasonable scale.
// Instead of a zoomed-in sub picture or a lot of empty space
fn discover_r_bound(
    grid_size: isize,
    orbital: &Orbital,
    sample_amount: u64,
) -> Result<f64, &'static str> {
    // We start as if the edges represent R_BOUND_MAX, which is very zoomed out.
    // We then calculate the average probability at each row and column.
    let delta_volume = (R_BOUND_MAX / (grid_size as f64 / 2.0 - 1.0)).powi(3);
    let norm_factor = R_BOUND_MAX as f64 / (grid_size as f64 / 2.0 - 1.0);

    let mut probs = vec![vec![0.0; grid_size as usize]; grid_size as usize];
    // ith row
    probs.par_iter_mut().enumerate().for_each(|(i, row_probs)| {
        let z = ((grid_size - 1) / 2 - i as isize) as f64 * norm_factor;
        // jth column
        row_probs.iter_mut().enumerate().for_each(|(j, prob)| {
            let x = (j as isize - (grid_size - 1) / 2) as f64 * norm_factor;
            let (r, theta, _) = spherical(x, 0.0, z);

            let mut p = orbital.probability(r, theta, ZERO_PHI_PLANE, delta_volume);

            p = 1.0 - (1.0 - p).powf(sample_amount as f64);

            *prob = p;
        });
    });

    let row_cum = probs
        .iter()
        .fold(Vec::new(), |mut acc: Vec<f64>, row: &Vec<f64>| {
            acc.push(row.into_iter().sum());
            acc
        });

    let col_cum = probs.iter().fold(probs[0].clone(), |acc, row| {
        let zipped = acc.into_iter().zip(row);
        zipped.map(|(a, b)| a + b).collect()
    });
    // At the current R_BOUND_MAX scale, the edges have a very low probability. We progress
    // from the edges to the centre to find the first row and column such that the average
    // probability is greater than PROB_THRESHOLD
    let row_bound: f64;
    let col_bound: f64;
    if let Some((row, _)) = row_cum
        .iter()
        .map(|&x| x / grid_size as f64)
        .enumerate()
        .find(|&(_, prob)| prob >= PROB_THRESHOLD)
    {
        let z = ((grid_size - 1) / 2 - row as isize) as f64 * norm_factor;
        row_bound = spherical(0.0, 0.0, z).0;
    } else {
        return Err("Max R bound is too small");
    }

    if let Some((col, _)) = col_cum
        .iter()
        .map(|&x| x / grid_size as f64)
        .enumerate()
        .find(|&(_, prob)| prob >= PROB_THRESHOLD)
    {
        let x = (col as isize - (grid_size - 1) / 2) as f64 * norm_factor;
        col_bound = spherical(x, 0.0, 0.0).0;
    } else {
        return Err("Max R bound is too small");
    }

    // Some of the orbitals are wider, some of them are taller. We need to fit according to the
    // longest axis
    Ok(if row_bound > col_bound {
        row_bound
    } else {
        col_bound
    })
}

fn spherical(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let r = (x * x + y * y + z * z).sqrt();
    let theta = if r == 0.0 { 0.0 } else { (z / r).acos() };
    let phi = if x == 0.0 { 0.0 } else { (y / x).acos() };

    (r, theta, phi)
}
