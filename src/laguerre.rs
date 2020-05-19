use num_integer::Integer;

/// Calculate r / (a * b), avoiding overflows and fractions.
///
/// Assumes that a * b divides r evenly.
/// Adapted from https://docs.rs/num-integer/0.1.42/src/num_integer/lib.rs.html#1077
fn divide_product<T: Integer + Clone>(r: T, a: T, b: T) -> T {
    // See http://blog.plover.com/math/choose-2.html for the idea.
    let g = num_integer::gcd(r.clone(), b.clone());
    (r / g.clone()) / (a * (b / g))
}

pub struct Laguerre {
    n: u64,
    coeffs: Vec<f64>,
}

impl Laguerre {
    pub fn new(n: u64, alpha: u64) -> Self {
        let facts: Vec<u64> = (0..=n + alpha)
            .scan(1, |prod, x| {
                *prod = *prod * if x == 0 { 1 } else { x };
                Some(*prod)
            })
            .collect();

        let fact = |i: u64| facts[i as usize];

        // TODO: divide_product may be better?
        let coeffs: Vec<f64> = (0..=n)
            .rev()
            .map(|j| {
                (fact(n + alpha) / (fact(j) * fact(alpha + (n - j)))) as f64 / fact(n - j) as f64
            })
            .collect();

        Laguerre { n, coeffs }
    }

    pub fn L(&self, x: f64) -> f64 {
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
