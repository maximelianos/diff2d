# Differentiable 2D render

Signed Distance Field is minimal distance from a point to a shape. It allows efficient compression of fonts and binary images, and also allows interesting blending of shapes.


## Auto differentiation

<img src=diagram-20230329.svg width=800>

Dtype - scalar, matrix
OpType - +, -, *
struct D - derivable
struct Dp - pointer to D
graph - computational graph

700 mb per 10 000 000 objects in graph

Bugs
1) Simple BFS results in incorrect gradients
2) const values are deleted on backward and result in dangling references

convert -delay 4 -loop 0 *.jpg myimage.gif

convert *.jpg -coalesce -duplicate 1,-2-1 -quiet -layers OptimizePlus  -loop 0 patrol_cycle.gif

## Derivative

$$ I_{ijc} = 1 \cdot step_1 * RGB_{c1} + \\ (1-step_1) step_2 * RGB_{c2} + \\ (1-(1-step_1)step_2)step_3 * RGB_{c3} + \ldots $$
$$ \frac{dI_{ijc}}{dstep_1} = wstep_1 (1\cdot RGB_{c1} + (-step_2)\cdot S_2 ) $$
$$ S_2 = 1\cdot RGB_{c2} + (-step_3)RGB_{c3} + (-step_3)(-step_4)RGB_{c4} +\ldots + (-step_3)(-step_4)\ldots (-step_n)RGB_{cn} $$

## References
SDF
https://habr.com/ru/post/438316/

~~Interior mutability~~
https://doc.rust-lang.org/book/ch15-05-interior-mutability.html

~~Global static variable~~
https://stackoverflow.com/questions/19605132/is-it-possible-to-use-global-variables-in-rust

~~Global variables~~
https://www.sitepoint.com/rust-global-variables/

No user-defined constructors in Rust
https://doc.rust-lang.org/nomicon/constructors.html

Inheritance struct common variables
https://users.rust-lang.org/t/struct-inheritance-embedding-best-practice/10627

Impl for alias type (inheritance)
https://stackoverflow.com/questions/35568871/is-it-possible-to-implement-methods-on-type-aliases

Deref for accessing parent type methods (inheritance)
https://stackoverflow.com/questions/45086595/is-it-considered-a-bad-practice-to-implement-deref-for-newtypes

Easy macros
https://blog.logrocket.com/macros-in-rust-a-tutorial-with-examples/

Cargo dependency resolution
http://aturon.github.io/tech/2018/07/25/cargo-version-selection/

nalgebra dot product
https://zsiciarz.github.io/24daysofrust/book/vol1/day14.html

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


