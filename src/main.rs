mod laguerre;
mod orbital;
mod sampler;

use image::{ImageBuffer, Rgb};
use sampler::Sampler;
use std::env;

const IMG_SIZE: usize = 4096;
const ITERS: u64 = 6_000_000_000;
fn main() {
    let args: Vec<String> = env::args().collect();

    let n: u64 = args[1].parse().unwrap();
    let l: u64 = args[2].parse().unwrap();
    let m: i64 = args[3].parse().unwrap();

    if n < 1 {
        panic!("n cannot be less than 1");
    }
    if l > n - 1 {
        panic!("l cannot be greater than n-1");
    }
    if m < -(l as i64) || m > l as i64 {
        panic!("l must be in range [-l,l]");
    }

    let container = vec![0; 3 * IMG_SIZE * IMG_SIZE];
    let mut image =
        ImageBuffer::<Rgb<u8>, _>::from_raw(IMG_SIZE as u32, IMG_SIZE as u32, container).unwrap();

    let mut sampler = Sampler::new(n, l, m, IMG_SIZE, ITERS);
    let grid = sampler.sample();

    let blue = Rgb([83, 202, 236]);
    let black = Rgb([0, 0, 0]);

    for (x, y, pixel) in image.enumerate_pixels_mut() {
        if grid[y as usize][x as usize] {
            *pixel = blue;
        } else {
            *pixel = black;
        }
    }

    image.save(format!("{}{}{}.png",n,l,m)).unwrap();
}
