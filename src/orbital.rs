use crate::laguerre::Laguerre;
use num_complex::Complex64;
use sphrs::{ComplexSHType, Coordinates, RealSHType, SHCoordinates, SHEval};
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

    pub fn n(&self) -> u64 {
        self.n
    }

    #[inline(always)]
    fn phase(&self, sph_harmonics: Complex64, coord: &Coordinates<f64>) -> Phase {
        // Phase calculation
        // This is the sign of R(r)Y_(m, l)(theta, phi), but Y_(m, l) is in its real form
        // (since we can't take the sign of a the complex number psi)

        // We use the Condon-Shortley phase convention for our definition of Y
        let condon_shortley_sign = if self.m.abs() % 2 == 0 { 1. } else { -1. };
        let r_sph_harm = if self.m < 0 {
            // Since we need to calculate Y_(-m,l) anyway, it's cheaper to calculate the real value directly
            condon_shortley_sign * RealSHType::Spherical.eval(self.l as i64, self.m, coord)
        } else if self.m == 0 {
            sph_harmonics.re
        } else {
            condon_shortley_sign * SQRT_2 * sph_harmonics.re
        };

        if r_sph_harm > 0.0 {
            Phase::Positive
        } else if r_sph_harm == 0.0 {
            Phase::Zero
        } else {
            Phase::Negative
        }
    }

    #[inline(always)]
    fn psi_with_phase(&self, coord: &Coordinates<f64>) -> (Complex64, Phase) {
        // let unit_sphere: Coordinates<f64> = Coordinates::spherical(1.0, coord.theta(), coord.phi());
        let rho = coord.r() * self.rho_over_r;
        let radial =
            self.root_term * (-rho / 2.0).exp() * rho.powi(self.l as i32) * self.laguerre.L(rho);
        let sph_harmonics = ComplexSHType::Spherical.eval(self.l as i64, self.m, coord);
        let psi = radial * sph_harmonics;

        let phase = self.phase(sph_harmonics, coord);

        (psi, phase)
    }

    #[inline(always)]
    fn psi(&self, coord: &Coordinates<f64>) -> Complex64 {
        // psi(r, theta, phi) = R(r)Y_(m, l)(theta, phi)
        // where R(r) is the real radial component, and Y_(m, l)(theta, phi) is the complex spherical harmonic
        let rho = coord.r() * self.rho_over_r;
        let radial =
            self.root_term * (-rho / 2.0).exp() * rho.powi(self.l as i32) * self.laguerre.L(rho);
        let psi = radial * ComplexSHType::Spherical.eval(self.l as i64, self.m, coord);

        psi
    }

    #[inline(always)]
    pub fn probability_with_phase(
        &self,
        coord: &Coordinates<f64>,
        delta_volume: f64,
    ) -> (f64, Phase) {
        let (psi_val, phase) = self.psi_with_phase(coord);

        let prob = psi_val.norm_sqr() * delta_volume;
        (if !prob.is_nan() { prob } else { 0. }, phase)
    }

    // |psi(r, theta, phi)|^2 is the probability per unit volume at (r, theta, phi). Multiply it by
    // the volume to get the probability to detect an electron in that region
    #[inline(always)]
    pub fn probability(&self, coord: &Coordinates<f64>, delta_volume: f64) -> f64 {
        let psi_val = self.psi(coord);

        let prob = psi_val.norm_sqr() * delta_volume;
        if !prob.is_nan() {
            prob
        } else {
            0.
        }
    }
}
