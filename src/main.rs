mod orbital;
mod laguerre;
mod sampler;

use sampler::Sampler;
use image::{ImageBuffer, Rgb};
use std::f64::consts::PI;

const IMG_SIZE:usize = 250;
const ITERS:usize = 10000;
fn main() {
    //let orbital = orbital::Orbital::new(2, 1, 1);
    //print!("{:?}\n", orbital.psi(1.0e-10,PI/4.0, PI));
    //print!("{:+e}\n", orbital.probability(1.0e-10,PI/4.0, PI,1.0));
    
    let container = vec![0;3*IMG_SIZE*IMG_SIZE];
    let mut image = ImageBuffer::<Rgb<u8>,_>::from_raw(IMG_SIZE as u32, IMG_SIZE as u32, container).unwrap();

    let mut sampler = Sampler::new(1,0,0, IMG_SIZE);
    let grid = sampler.sample(ITERS);

    let blue = Rgb([0,0,255]);
    let black = Rgb([0,0,0]);

    for (x, y, pixel) in image.enumerate_pixels_mut() {
        if grid[y as usize][x as usize] {
            *pixel = blue;
        } else {
            *pixel = black;
        }
    }

    image.save("render.png").unwrap();
}