mod laguerre;
mod orbital;
mod sampler;

use image::{imageops, FilterType, ImageBuffer, Rgb};
use std::env;

const IMG_SIZE: usize = 2048;
const SAMPLING_SIZE: usize = IMG_SIZE * 4;
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

    let container = vec![0; 3 * SAMPLING_SIZE * SAMPLING_SIZE];
    let mut image =
        ImageBuffer::<Rgb<u8>, _>::from_raw(SAMPLING_SIZE as u32, SAMPLING_SIZE as u32, container)
            .unwrap();

    let grid = sampler::sample(n, l, m, SAMPLING_SIZE);

    let blue = Rgb([83, 202, 236]);
    let black = Rgb([0, 0, 0]);
    let red = Rgb([236, 91, 83]);

    for (x, y, pixel) in image.enumerate_pixels_mut() {
        match grid[y as usize][x as usize] {
            orbital::Phase::Positive => *pixel = blue,
            orbital::Phase::Negative => *pixel = red,
            orbital::Phase::Zero => *pixel = black,
        }
    }

    image = imageops::resize(
        &image,
        IMG_SIZE as u32,
        IMG_SIZE as u32,
        FilterType::Gaussian,
    );
    image.save(format!("{}{}{}.png", n, l, m)).unwrap();
}
