pub struct Laguerre {
    n: u64,
    coeffs: Vec<f64>,
}

// all factorials representible in u64, starting from 0!
const FACTS: [u64; 21] = [
    1,
    1,
    2,
    6,
    24,
    120,
    720,
    5040,
    40320,
    362880,
    3628800,
    39916800,
    479001600,
    6227020800,
    87178291200,
    1307674368000,
    20922789888000,
    355687428096000,
    6402373705728000,
    121645100408832000,
    2432902008176640000,
];

#[inline(always)]
fn fact(i: u64) -> u64 {
    FACTS[i as usize]
}
// Generalised Laguerre polynomial
// https://en.wikipedia.org/wiki/Laguerre_polynomials#Generalized_Laguerre_polynomials
//
// Calculated using the closed summation form
impl Laguerre {
    pub fn new(n: u64, alpha: u64) -> Self {
        // precompute the coefficient terms ((n + alpha) Choose (n - i)) / i!
        // TODO: divide_product may be better? But this is not a hot spot
        let coeffs: Vec<f64> = (0..=n)
            .rev()
            .map(|j| {
                (fact(n + alpha) / (fact(j) * fact(alpha + (n - j)))) as f64 / fact(n - j) as f64
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
