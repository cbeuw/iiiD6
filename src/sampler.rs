use crate::orbital;

use rand::{Rng, SeedableRng};
use rand::rngs::SmallRng;
use rayon::prelude::*;

const DEFAULT_R_BOUND:f64 = 10e-10; // 4  angstrom
//const DEFAULT_DELTA_VOLUME:f64 = 1e-30; // 1 cubic angstrom
const ZERO_PHI_PLANE:f64 = 0.0;
const PROB_AMPLIFICATION:f64 = 50.0;
pub struct Sampler{
    grid: Vec<Vec<bool>>,   
    probs: Vec<Vec<f64>>
}

impl Sampler{
    pub fn new(n: u64, l:u64, m:i64, grid_size:usize) -> Self{
        let mut probs = Vec::new();
        let delta_volume =  (DEFAULT_R_BOUND / (grid_size as f64 / 2.0 - 1.0)).powi(3);

        let orbital = orbital::Orbital::new(n,l,m);

        let grid_isize = grid_size as isize;
        let norm_factor =DEFAULT_R_BOUND as f64/ (grid_size as f64 / 2.0 - 1.0);
        for i in 0..grid_isize{
            let z = (grid_isize/2 - i) as f64 * norm_factor;
            let mut row_probs = Vec::new();
            for j in 0..grid_isize{
                let x = (j - grid_isize/2) as f64 * norm_factor ;

                let r = (x*x+z*z).sqrt();
                //let theta = if z != 0.0 {(x/z).atan()} else {0.0};
                let theta = if (j - grid_isize/2) == 0 || (grid_isize/2 - i) == 0 {0.0} else {(z/r).acos()};
                let mut prob = orbital.probability(r,theta,ZERO_PHI_PLANE,delta_volume);

                prob *= PROB_AMPLIFICATION;

                if prob < 0.0 || prob > 1.0 {
                    panic!(format!("bad prob {}", prob));
                }

                row_probs.push(prob);
            }
            //println!("({}): {}",i,row_probs.iter().cloned().fold(0./0., f64::max));
            probs.push(row_probs);
        }

        Sampler{
            grid: vec![vec![false;grid_size];grid_size],

            probs
        }
    }

    pub fn sample(&mut self, iter_per_row: usize) -> &Vec<Vec<bool>>{
        let probs = &self.probs;
        self.grid.par_iter_mut().enumerate().for_each(|(i,row)|{
            let mut rng = SmallRng::from_entropy();
            for _ in 0..iter_per_row {
                for (j, hit) in row.iter_mut().enumerate() {
                    *hit |= rng.gen_bool(probs[i][j]);
                };
            }
        });
        &self.grid
    }

    pub fn sample_mono(&mut self, iter_per_row: usize) -> &Vec<Vec<bool>>{
        let probs = &self.probs;
        self.grid.par_iter_mut().enumerate().for_each(|(i,row)|{
            let mut rng = SmallRng::from_entropy();
                for (j, hit) in row.iter_mut().enumerate() {
                    let prob = 1.0-(1.0-probs[i][j]).powi(iter_per_row as i32);
                    *hit |= rng.gen_bool(prob);
                };
        });
        &self.grid
    }
}