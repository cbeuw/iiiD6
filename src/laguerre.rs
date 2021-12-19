use crate::factorial::fact;
pub struct Laguerre {
    n: u64,
    coeffs: Vec<f64>,
}

// Generalised Laguerre polynomial
// https://en.wikipedia.org/wiki/Laguerre_polynomials#Generalized_Laguerre_polynomials
//
// Calculated using the closed summation form
impl Laguerre {
    pub fn new(n: u64, α: u64) -> Self {
        // precompute the coefficient terms ((n + alpha) Choose (n - i)) / i!
        // TODO: divide_product may be better? But this is not a hot spot
        let coeffs: Vec<f64> = (0..=n)
            .rev()
            .map(|j| {
                (fact(n + α) / (fact(j) * fact(α + (n - j)))) as f64 / fact(n - j) as f64
            })
            .collect();

        Laguerre { n, coeffs }
    }

    pub fn L(&self, x: f64) -> f64 {
        // do the summation for a particular x
        (0..=self.n).fold(0.0, |sum: f64, i| {
            let delta: f64 = self.coeffs[i as usize] * x.powi(i as i32);
            if i % 2 == 0 {
                sum + delta
            } else {
                sum - delta
            }
        })
    }
}
