use crate::orbital;

use rand::{Rng, SeedableRng};
use rand::rngs::SmallRng;
use rayon::prelude::*;

const DEFAULT_R_BOUND:f64 = 4e-10; // 4 angstrom
//const DEFAULT_DELTA_VOLUME:f64 = 1e-30; // 1 cubic angstrom
const ZERO_PHI_PLANE:f64 = 0.0;
pub struct Sampler{
    delta_volume: f64,
    phi_plane: f64,

    orbital: orbital::Orbital,
    grid: Vec<Vec<bool>>,

    coords: Vec<Vec<(f64,f64)>>
}

impl Sampler{
    pub fn new(n: u64, l:u64, m:i64, grid_size:usize) -> Self{
        let norm_factor =DEFAULT_R_BOUND as f64/ (grid_size as f64 / 2.0);

        let mut coords = Vec::new();

        let grid_isize = grid_size as isize;

        for i in 0..grid_isize{
            let z = (grid_isize/2 - i) as f64 * norm_factor;
            let mut row_coords = Vec::new();
            for j in 0..grid_isize{
                let x = (j - grid_isize/2) as f64 * norm_factor ;

                let r = (x*x+z*z).sqrt();
                let theta = (x as f64/z as f64).atan();
                row_coords.push((r,theta));
            }
            coords.push(row_coords);
        }

        Sampler{
            delta_volume: (DEFAULT_R_BOUND / (grid_size as f64 / 2.0)).powi(3),
            phi_plane: ZERO_PHI_PLANE,

            orbital: orbital::Orbital::new(n,l,m),
            grid: vec![vec![false;grid_size];grid_size],

            coords
        }
    }

    pub fn sample(&mut self, iter_per_row: usize) -> &Vec<Vec<bool>>{
        let orbital = &self.orbital;
        let coords = &self.coords;
        let phi = self.phi_plane;
        //let volume = DEFAULT_DELTA_VOLUME;
        let volume = self.delta_volume;
        //let mut val = String::new();
        self.grid.par_iter_mut().enumerate().for_each(|(i,row)|{
            let mut rng = SmallRng::from_entropy();
            //let mut largest = 0.0;
            for _ in 0..iter_per_row {
                for (j, hit) in row.iter_mut().enumerate() {
                    let (r,theta) = coords[i][j];
                    let prob = orbital.probability(r,theta,phi,volume);
                    //val = format!("{} {1:+e}",val, prob);
                    //if prob > largest { largest = prob};
                    //if prob < smallest{smallest = prob};
                    *hit |= rng.gen_bool(prob);
                };
            }
            //println!("{}",val);
        });
        //println!("{}", largest);
        //println!("{}", smallest);
        &self.grid
    }
}