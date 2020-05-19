mod laguerre;
mod orbital;
mod sampler;

use image::{ImageBuffer, Rgb};
use sampler::Sampler;

const IMG_SIZE: usize = 4096;
const ITERS: usize = 2000000000;
fn main() {
    let container = vec![0; 3 * IMG_SIZE * IMG_SIZE];
    let mut image =
        ImageBuffer::<Rgb<u8>, _>::from_raw(IMG_SIZE as u32, IMG_SIZE as u32, container).unwrap();

    let mut sampler = Sampler::new(4, 2, 0, IMG_SIZE, ITERS);
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

    image.save("render.png").unwrap();
}
