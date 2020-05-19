use crate::laguerre;
use sphrs::{ComplexSHType, Coordinates, SHEval};
use num_complex::Complex64;
use std::f64::consts;

const FRAC_1_REDUCED_BOHR_RADIUS:f64 = 1.888706223265207e10;
const FRAC_1_REDUCED_BOHR_RADIUS_CUBED:f64 = 6.73741398899e30;
/*
struct Orbital {
    n: u64,
    l: u64,
    m: i64,

    root_term: f64,
    laguerre_L: Box<dyn Fn(f64) -> f64>,
}

impl Orbital{
    fn new(n:u64, l:u64, m:i64) -> &Self {
        let facts: Vec<_> = (0..=n+l).scan(1, |prod, x|{
            *prod = *prod * if x==0 {1} else {x};
            Some(*prod)
        }).collect();
    
        let root_term = (
            FRAC_1_REDUCED_BOHR_RADIUS_CUBED
            *
            (8 * facts[(n-l-1) as usize]) as f64
            /
            (n*n*n*n * 2 * facts[(n+l) as usize]) as f64
        ).sqrt();
    
        let laguerre_L = laguerre::laguerre_polynomial(n-l-1, 2*l+1);
    
        &Orbital{
            n,
            l,
            m,

            root_term,
            laguerre_L: Box::new(laguerre_L),
        }
    }
}
*/

pub fn wavefuncion(n: u64, l: u64, m:i64) -> impl Fn(f64,f64,f64) -> Complex64{
    let facts: Vec<_> = (0..=n+l).scan(1, |prod, x|{
        *prod = *prod * if x==0 {1} else {x};
        Some(*prod)
    }).collect();

    let root_term = (
        FRAC_1_REDUCED_BOHR_RADIUS_CUBED
        *
        (8 * facts[(n-l-1) as usize]) as f64
        /
        (n*n*n*n * 2 * facts[(n+l) as usize]) as f64
    ).sqrt();

    let laguerre_L = laguerre::laguerre_polynomial(n-l-1, 2*l+1);
    
    move |r:f64, theta: f64, phi:f64| {
        let rho = (r * 2.0 / n as f64) * FRAC_1_REDUCED_BOHR_RADIUS;
        let unit_sphere: Coordinates<f64> = Coordinates::spherical(1.0, theta, phi);
        print!("rou: {}\n", rho);
        print!("laguerre: {}\n",laguerre::laguerre_polynomial(n-l-1, 2*l+1)(rho));
        print!("harmonics: {}\n", ComplexSHType::Spherical.eval(l as i64,m as i64,&unit_sphere));


        Complex64::new(
        root_term                   *
        consts::E.powf(-rho / 2.0)  *
        rho.powi(l as i32)          *
        laguerre_L(rho),
        0.0) * 
        ComplexSHType::Spherical.eval(l as i64,m as i64,&unit_sphere)
    }
}