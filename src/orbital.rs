
use crate::laguerre;
use sphrs::{ComplexSHType, Coordinates, SHEval};
use num_complex::Complex64;
use std::f64::consts;

pub const REDUCED_BOHR_RADIUS:f64 = 5.294651e-11;

pub struct Orbital{
    n:u64,
    l:u64,
    m:i64,

    laguerre: laguerre::Laguerre,
    root_term: f64,
}

impl Orbital{
    pub fn new(n: u64, l:u64, m:i64) -> Self{
        let facts: Vec<_> = (0..=n+l).scan(1, |prod, x|{
            *prod = *prod * if x==0 {1} else {x};
            Some(*prod)
        }).collect();
    
        let root_term = (
            (8 * facts[(n-l-1) as usize]) as f64
            /
            ((n*n*n*n * 2 * facts[(n+l) as usize]) as f64 * REDUCED_BOHR_RADIUS*REDUCED_BOHR_RADIUS*REDUCED_BOHR_RADIUS)
        ).sqrt();
    
        let laguerre = laguerre::Laguerre::new(n-l-1, 2*l+1);

        Orbital{
            n,l,m,
            laguerre,
            root_term
        }
    }

    pub fn psi(&self,r:f64, theta:f64, phi:f64) -> Complex64{
        let rho = r * 2.0 / (self.n as f64 * REDUCED_BOHR_RADIUS);
        let unit_sphere: Coordinates<f64> = Coordinates::spherical(1.0, theta, phi);

        Complex64::new(
            self.root_term              *
            consts::E.powf(-rho / 2.0)  *
            rho.powi(self.l as i32)     *
            self.laguerre.L(rho)
            ,
            0.0) 
        *
        ComplexSHType::Spherical.eval(self.l as i64,self.m as i64,&unit_sphere)
    }

    pub fn probability(&self,r:f64, theta:f64, phi:f64, delta_V:f64) -> f64 {
        self.psi(r,theta,phi).norm_sqr() * delta_V
    }
}