use crate::laguerre::Laguerre;
use num_complex::Complex64;
use sphrs::{ComplexSHType, Coordinates, SHCoordinates, SHEval};
use std::f64::consts::SQRT_2;

const REDUCED_BOHR_RADIUS: f64 = 5.294651e-11;

#[derive(Copy, Clone, Debug)]
pub enum Phase {
    Positive,
    Zero,
    Negative,
}

pub struct Orbital {
    n: u64,
    l: u64,
    m: i64,

    rho_over_r: f64,
    laguerre: Laguerre,
    root_term: f64,
}

impl Orbital {
    // wavefunction as given in https://en.wikipedia.org/wiki/Hydrogen_atom#Wavefunction
    pub fn new(n: u64, l: u64, m: i64) -> Self {
        // precompute all factorials up to n+l (which is the largest factorial we need)
        let facts: Vec<_> = (0..=n + l)
            .scan(1, |prod, x| {
                *prod = *prod * if x == 0 { 1 } else { x };
                Some(*prod)
            })
            .collect();

        // some terms do not require any of r, theta or phi, so we precompute them here
        let root_term = ((8 * facts[(n - l - 1) as usize]) as f64
            / ((n * n * n * n * 2 * facts[(n + l) as usize]) as f64
                * REDUCED_BOHR_RADIUS
                * REDUCED_BOHR_RADIUS
                * REDUCED_BOHR_RADIUS))
            .sqrt();

        let laguerre = Laguerre::new(n - l - 1, 2 * l + 1);
        let rho_over_r = 2.0 / (n as f64 * REDUCED_BOHR_RADIUS);

        Orbital {
            n,
            l,
            m,
            laguerre,
            root_term,
            rho_over_r,
        }
    }

    #[inline(always)]
    pub fn psi(&self, coord: Coordinates<f64>) -> (Complex64, Phase) {
        let rho = coord.r() * self.rho_over_r;
        let unit_sphere: Coordinates<f64> = Coordinates::spherical(1.0, coord.theta(), coord.phi());

        let radial =
            self.root_term * (-rho / 2.0).exp() * rho.powi(self.l as i32) * self.laguerre.L(rho);

        let c_sph_harm = ComplexSHType::Spherical.eval(self.l as i64, self.m as i64, &unit_sphere);
        //let r_sph_harm: f64 = RealSHType::Spherical.eval(self.l as i64, self.m as i64, &unit_sphere);

        let sign = if self.m % 2 == 0 { 1. } else { -1. };
        let r_sph_harm = if self.m < 0 {
            sign * SQRT_2
                * ComplexSHType::Spherical
                    .eval(self.l as i64, -self.m as i64, &unit_sphere)
                    .im
        } else if self.m == 0 {
            c_sph_harm.re
        } else {
            sign * SQRT_2 * c_sph_harm.re
        };

        let phase = if radial * r_sph_harm > 0.0 {
            Phase::Positive
        } else if radial * r_sph_harm == 0.0 {
            Phase::Zero
        } else {
            Phase::Negative
        };

        (radial * c_sph_harm, phase)
    }

    // |psi(r, theta, phi)|^2 is the probability per unit volume at (r, theta, phi). Multiply it by
    // the volume to get the probability to detect an electron in that region
    #[inline(always)]
    pub fn probability(&self, coord: Coordinates<f64>, delta_volume: f64) -> (Phase, f64) {
        let (psi_val, phase) = self.psi(coord);

        let prob = psi_val.norm_sqr() * delta_volume;
        (phase, if !prob.is_nan() { prob } else { 0. })
    }
}
