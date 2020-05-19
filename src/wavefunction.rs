use crate::laguerre;
use sphrs::{ComplexSHType, Coordinates, SHEval};
use num_complex::Complex64;
use std::f64::consts;

const FRAC_1_REDUCED_BOHR_RADIUS_ATTO:f64 = 1.888697507170487e-8;
const FRAC_1_REDUCED_BOHR_RADIUS_CUBED_ATTO:f64 = 6.737320712965956e-24;
pub const REDUCED_BOHR_RADIUS_ATTO:u64 = 52_946_541;

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
            FRAC_1_REDUCED_BOHR_RADIUS_CUBED_ATTO
            *
            (8 * facts[(n-l-1) as usize]) as f64
            /
            (n*n*n*n * 2 * facts[(n+l) as usize]) as f64
        ).sqrt();
    
        let laguerre = laguerre::Laguerre::new(n-l-1, 2*l+1);

        Orbital{
            n,l,m,
            laguerre,
            root_term
        }
    }

    pub fn psi(&self,r:f64, theta:f64, phi:f64) -> Complex64{
        let rho = (r * 2.0 / self.n as f64) * FRAC_1_REDUCED_BOHR_RADIUS_ATTO;
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

    pub fn probability(&self,r:f64, theta:f64, phi:f64) -> f64 {
        self.psi(r,theta,phi).norm_sqr() * (r * r) * theta.sin()
    }
}