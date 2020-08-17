# 3d6
A hydrogen orbital wave function renderer based on random sampling, written in Rust

## Usage
```
$ cargo build --release
$ target/release/iiiD6 [n] [l] [m]
```

Where n, l, m are quantum numbers described here: https://en.wikipedia.org/wiki/Atomic_orbital#Quantum_numbers. 

They can be any integer that satisfies
```
n >= 1
0 <= l <= n-1
-l <= m <= l
```
Although the scaling currenlty doesn't work very well for large n (>4), so `R_BOUND_MAX` in `sampler.rs` may need to be adjusted manually.

## Renders
These can be found in `render` folder
### 200
![200](https://github.com/cbeuw/iiiD6/blob/master/render/200.png)
### 310
![310](https://github.com/cbeuw/iiiD6/blob/master/render/310.png)
### 320
![320](https://github.com/cbeuw/iiiD6/blob/master/render/320.png)
### 410
![410](https://github.com/cbeuw/iiiD6/blob/master/render/410.png)
### 420 🌳
![420](https://github.com/cbeuw/iiiD6/blob/master/render/420.png)
### 422
![422](https://github.com/cbeuw/iiiD6/blob/master/render/422.png)
### 430
![430](https://github.com/cbeuw/iiiD6/blob/master/render/430.png)

## Prints & posters
I've uploaded some of the renders on Redbubble. They make really pretty posters and clock surfaces! https://www.redbubble.com/people/cbeuw/explore
