use crate::orbital::{Orbital, Phase};

use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use sphrs::Coordinates;

const ZERO_Y_PLANE: f64 = 0.0;

pub fn sample(n: u64, l: u64, m: i64, grid_size: usize) -> Vec<Vec<Phase>> {
    let orbital = Orbital::new(n, l, m);
    let grid_isize = grid_size as isize;

    let (r_bound, sample_amount) = discover_r_bound_and_iter(grid_isize, &orbital).unwrap();

    // each pixel represents a length in space. we take the volume as a cube with side length of that
    // this volume is multiplied by |psi(r, theta, phi)|^2 at each point in space to get the probability
    // at the volume around that point
    let delta_volume = (r_bound / (grid_size as f64 / 2.0 - 1.0)).powi(3);

    // x and z are currently in range of [-grid_size/2, grid_size/2], we normalise that to [-r_bound, r_bound]
    // with the normalisation factor
    let norm_factor = r_bound as f64 / (grid_size as f64 / 2.0 - 1.0);

    let mut grid = vec![vec![Phase::Zero; grid_size]; grid_size];
    grid.par_iter_mut().enumerate().for_each(|(i, row)| {
        let mut rng = SmallRng::from_entropy();
        let z = ((grid_isize - 1) / 2 - i as isize) as f64 * norm_factor;
        row.iter_mut().enumerate().for_each(|(j, cell)| {
            let x = (j as isize - (grid_isize - 1) / 2) as f64 * norm_factor;

            let coord = Coordinates::cartesian(x, ZERO_Y_PLANE, z);

            let (mut p, phase) = orbital.probability_with_phase(&coord, delta_volume);

            // This calculates the probability at this point after sample_amount of sampling
            p = 1.0 - (1.0 - p).powf(sample_amount as f64);

            if rng.gen_bool(p) {
                *cell = phase;
            }
        });
    });
    grid
}

// Find the maximum r_bound value so that we can see the entire orbital shape at a reasonable scale.
// Instead of a zoomed-in sub picture or a lot of empty space. Also finds the appropirate amount of
// "iterations" so that the regions with the highest probability are the right brightness
fn discover_r_bound_and_iter(
    grid_size: isize,
    orbital: &Orbital,
) -> Result<(f64, u64), &'static str> {
    // Predicted r bound using an empirical quadratic regression
    let pred_r_bound: f64 = 1.0446e-10 * (orbital.n() * orbital.n()) as f64
        + 1.7629e-10 * orbital.n() as f64
        - 7.5457e-12;
    let r_bound_max = pred_r_bound * 2.;

    // We want to determine the amount of "iterations" we sample the entire grid.
    // We want, at this amount of iterations, the average probability at top 12.5% of rows
    // and columns to be AVG_PROB_AT_MAX_RUN
    const AVG_PROB_AT_MAX_RUN: f64 = 0.999;

    // We want to find an r_bound such that
    // the minimum of (average probability at each row, average probability at each column)
    // after iter_amount of samplings greater than PROB_THRESHOLD
    const PROB_THRESHOLD: f64 = 0.08;

    // We start as if the edges represent r_bound_max, which is very zoomed out.
    // We then calculate the average probability at each row and column.
    let delta_volume = (r_bound_max / (grid_size as f64 / 2.0 - 1.0)).powi(3);
    // norm_factor is multiplied onto x and z coordinates of the grid such that the result
    // at the edges equal to r_bound_max metre
    let norm_factor = r_bound_max as f64 / (grid_size as f64 / 2.0 - 1.0);

    let mut probs = vec![vec![0.0; grid_size as usize]; grid_size as usize];
    let sample_interval = 2;
    // ith row
    probs
        .par_iter_mut()
        .enumerate()
        .filter(|(x, _)| x % sample_interval == 0)
        .for_each(|(i, row_probs)| {
            let z = ((grid_size - 1) / 2 - i as isize) as f64 * norm_factor;
            // jth column
            row_probs
                .iter_mut()
                .enumerate()
                .filter(|(x, _)| x % sample_interval == 0)
                .for_each(|(j, prob)| {
                    let x = (j as isize - (grid_size - 1) / 2) as f64 * norm_factor;

                    let coord = Coordinates::cartesian(x, ZERO_Y_PLANE, z);

                    let p = orbital.probability(&coord, delta_volume);
                    *prob = p;
                });
        });

    // Here we calculate the average single-sample probability at each row and each column
    let row_avg: Vec<f64> = probs
        .iter()
        .fold(Vec::new(), |mut acc: Vec<f64>, row: &Vec<f64>| {
            acc.push(row.into_iter().sum());
            acc
        })
        .iter()
        .map(|&x| x / (grid_size / sample_interval as isize) as f64)
        .collect();

    let col_avg: Vec<f64> = probs
        .iter()
        .fold(probs[0].clone(), |acc, row| {
            let zipped = acc.into_iter().zip(row);
            zipped.map(|(a, b)| a + b).collect()
        })
        .iter()
        .map(|&x| x / (grid_size / sample_interval as isize) as f64)
        .collect();

    // We then calculate the average probabilites of top grid_size/8 rows and columns
    let mut row_avg_sorted = row_avg.clone();
    row_avg_sorted.sort_unstable_by(|a, b| b.partial_cmp(a).unwrap());
    let mut col_avg_sorted = col_avg.clone();
    col_avg_sorted.sort_unstable_by(|a, b| b.partial_cmp(a).unwrap());

    let limit = (grid_size / 8) as usize;
    let row_top_sum: f64 = row_avg_sorted[..limit].iter().sum();
    let col_top_sum: f64 = col_avg_sorted[..limit].iter().sum();
    let top_avg = (row_top_sum + col_top_sum) / ((limit * 2) as f64);

    // 1-(1-top_avg) ^ (sample_amount) = AVG_PROB_AT_MAX_RUN
    // Solve for sample_amount, get sample_amount = ln(1-top_avg) / ln(1-AVG_PROB_AT_MAX_RUN)
    // which is log base (1-top_avg) of (1-AVG_PROB_AT_MAX_RUN)
    let sample_amount: u64 = (1.0 - AVG_PROB_AT_MAX_RUN).log(1.0 - top_avg).round() as u64;

    // Now we know how many iterations we want to sample, we can then work out where the edges of
    // our grid should be
    //
    // At the current R_BOUND_MAX scale, the edges have a very low probability. We progress
    // from the edges to the centre to find the first row and column such that the average
    // probability is greater than PROB_THRESHOLD
    let row_bound: f64;
    let col_bound: f64;
    if let Some((row, _)) = row_avg
        .iter()
        .enumerate()
        .find(|&(_, prob)| 1.0 - (1.0 - prob).powf(sample_amount as f64) >= PROB_THRESHOLD)
    {
        if row == grid_size as usize || row == 0 {
            println!("Max R bound is too small: row limit at {}", row)
        }
        let z = ((grid_size - 1) / 2 - row as isize) as f64 * norm_factor;
        row_bound = z.abs();
    } else {
        return Err("No row satisfies minimum probability threshold");
    }

    if let Some((col, _)) = col_avg
        .iter()
        .enumerate()
        .find(|&(_, prob)| 1.0 - (1.0 - prob).powf(sample_amount as f64) >= PROB_THRESHOLD)
    {
        if col == grid_size as usize || col == 0 {
            println!("Max R bound is too small: row limit at {}", col)
        }
        let x = (col as isize - (grid_size - 1) / 2) as f64 * norm_factor;
        col_bound = x.abs();
    } else {
        return Err("No column satisfies minimum probability threshold");
    }

    // Some of the orbitals are wider, some of them are taller. We need to fit according to the
    // longest axis
    Ok((f64::max(row_bound, col_bound), sample_amount))
}
