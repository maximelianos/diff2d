# Differentiable 2D render
<img src=pictures/fractal0.png width=400>

Signed Distance Field is minimal distance from a point to a shape. It allows efficient compression of fonts and binary images, and also allows interesting blending of shapes.

## Auto differentiation

Dtype - scalar, matrix
OpType - +, -, *
struct D - derivable
struct Dp - pointer to D
graph - computational graph

## References
Interior mutability
https://doc.rust-lang.org/book/ch15-05-interior-mutability.html

Global static variable
https://stackoverflow.com/questions/19605132/is-it-possible-to-use-global-variables-in-rust

Global variables
https://www.sitepoint.com/rust-global-variables/

No class constructors
https://doc.rust-lang.org/nomicon/constructors.html

SDF
https://habr.com/ru/post/438316/

## Changelog

v1
* Draw a circle, optimize its position and color
* 1 sample per pixel
* Add dif.rs executable, implement autodiff for scalars with simple operations

## Features


## Implementation
Written from scratch in Rust. Implements basic 2D geometry: 
* Point struct (which serves as a vector), implemented dot product, cross product, projection of one vector onto another, rotation
* Shape struct for all shapes - both analytical, loaded SDF from PNG or computed from binary image via BFS graph algorithm
* ShapeTransform struct, applies scaling, rotation, translation (in this order). ShapeTransform::apply actually computes inverse transform to go from desired point to local shape space. Scaling is also applied to SDF values
* ShapeType denotes analytical shapes (circle, rectangle) or bitmap - SDF loaded from PNG file or computed from binary image


## Files
* Rust sources are src/main.rs and src/loadsdf.rs
* Cargo.toml - dependecies
* resources/ - SDF and images in png format

## Compile and run
Project was developed on Ubuntu OS, using Rust compiler rustc 1.67.1. Compile and run:
```
cargo build --release && ./target/release/sdf
```
Cargo should install Rust dependencies. Rendered images will be saved to working directory. Implemented scenes:


