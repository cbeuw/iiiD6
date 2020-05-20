use crate::orbital::Orbital;

use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

const ZERO_PHI_PLANE: f64 = 0.0;
const R_BOUND_MAX: f64 = 40e-10; // 40 angstrom
const PROB_THRESHOLD_AT_R_MAX: f64 = 0.10; // 10%
pub struct Sampler {
    grid: Vec<Vec<bool>>,
    probs: Vec<Vec<f64>>,
}

impl Sampler {
    pub fn new(n: u64, l: u64, m: i64, grid_size: usize, sample_amount: u64) -> Self {
        let orbital = Orbital::new(n, l, m);
        let grid_isize = grid_size as isize;

        let r_bound = Self::discover_R_bound(grid_isize, &orbital, sample_amount).unwrap();
        let delta_volume = (r_bound / (grid_size as f64 / 2.0 - 1.0)).powi(3);
        let norm_factor = r_bound as f64 / (grid_size as f64 / 2.0 - 1.0);

        let mut probs = vec![vec![0.0; grid_size]; grid_size];
        probs.par_iter_mut().enumerate().for_each(|(i, row_probs)| {
            let z = ((grid_isize - 1) / 2 - i as isize) as f64 * norm_factor;
            row_probs.iter_mut().enumerate().for_each(|(j, prob)| {
                let x = (j as isize - (grid_isize - 1) / 2) as f64 * norm_factor;

                let (r, theta, _) = Self::spherical(x, 0.0, z);

                let mut p = orbital.probability(r, theta, ZERO_PHI_PLANE, delta_volume);

                p = 1.0 - (1.0 - p).powf(sample_amount as f64);

                if p < 0.0 || p > 1.0 {
                    panic!(format!("bad prob {}", p));
                }

                *prob = p;
            });
        });

        Sampler {
            grid: vec![vec![false; grid_size]; grid_size],

            probs,
        }
    }

    fn discover_R_bound(
        grid_size: isize,
        orbital: &Orbital,
        sample_amount: u64,
    ) -> Result<f64, &'static str> {
        let delta_volume = (R_BOUND_MAX / (grid_size as f64 / 2.0 - 1.0)).powi(3);
        let norm_factor = R_BOUND_MAX as f64 / (grid_size as f64 / 2.0 - 1.0);

        let mut row_cum = vec![0.0; grid_size as usize];
        let mut col_cum = vec![0.0; grid_size as usize];
        // ith row
        for i in 0..grid_size {
            let z = ((grid_size - 1) / 2 - i) as f64 * norm_factor;
            // jth column
            for j in 0..grid_size {
                let x = (j - (grid_size - 1) / 2) as f64 * norm_factor;
                let (r, theta, _) = Self::spherical(x, 0.0, z);

                let mut prob = orbital.probability(r, theta, ZERO_PHI_PLANE, delta_volume);

                prob = 1.0 - (1.0 - prob).powf(sample_amount as f64);

                row_cum[i as usize] += prob;
                col_cum[j as usize] += prob;
            }
        }

        let row_bound: f64;
        let col_bound: f64;
        if let Some((row, _)) = row_cum
            .iter()
            .map(|&x| x / grid_size as f64)
            .enumerate()
            .find(|&(_, prob)| prob >= PROB_THRESHOLD_AT_R_MAX)
        {
            let z = ((grid_size - 1) / 2 - row as isize) as f64 * norm_factor;
            let (r, _, _) = Self::spherical(0.0, 0.0, z);
            row_bound = r;
        } else {
            return Err("Max R bound is too small");
        }
        if let Some((col, _)) = col_cum
            .iter()
            .map(|&x| x / grid_size as f64)
            .enumerate()
            .find(|&(_, prob)| prob >= PROB_THRESHOLD_AT_R_MAX)
        {
            let x = (col as isize - (grid_size - 1) / 2) as f64 * norm_factor;
            let (r, _, _) = Self::spherical(x, 0.0, 0.0);
            col_bound = r;
        } else {
            return Err("Max R bound is too small");
        }

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

    pub fn sample(&mut self) -> &Vec<Vec<bool>> {
        let probs = &self.probs;
        self.grid.par_iter_mut().enumerate().for_each(|(i, row)| {
            let mut rng = SmallRng::from_entropy();
            for (j, hit) in row.iter_mut().enumerate() {
                *hit |= rng.gen_bool(probs[i][j]);
            }
        });
        &self.grid
    }
}
