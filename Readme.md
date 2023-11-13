# Differentiable 2D SDF render

https://github.com/maximelianos/diff2d

Signed Distance Field is minimal distance from a point to a shape. It allows efficient compression of fonts and binary images, and also allows interesting blending of shapes.

In differentiable rendering the scene parameters are optimized by minimizing Mean Squared Error with a reference image.

## Features
* Analytical SDF forward and backward rendering
 * Drawing and differentiating SDF by shape parameters for circle, rectangle, triangle, bitmap SDF. Geometric and color parameters are optimized
 * More complex SDF: 5-angled star, horseshoe, moon
 * Adam optimizer
 * Differentiated alpha blending
 * Bitmap SDF is optimized per-pixel and is single-color. Multiple bitmaps allow multi-color SDF by using differentiated alpha blending
 * Fast derivatives computed by hand
 * Slow automatic self-written differetiation with forward and backward passes
* Textured mesh edge sampling - differentiation w.r.t. texture pixels and triangle positions using Reynolds theorem
* SDF edge sampling, same approach as with textured mesh

<img src=pictures/sgdmomentum.gif width=300> <img src=pictures/adam.gif width=300>

<img src=pictures/tri.gif width=300> <img src=pictures/tricolor.gif width=300>

<img src=pictures/sdf.webp width=300> <img src=pictures/complex-sdf.webp width=300>

<img src=pictures/mesh-divergence.webp width=300> <img src=pictures/sdf-sampling.webp width=300>

The mesh edge sampling optimization works only if texture learning rate is lowered while positions are optimized. Over time it diverges, also if verticies fall outside the screen, they cannot return and decrease effective texture area.

## Auto differentiation

Using forward and backward pass, construct a computation graph in runtime and traverse verticies in same order backward, computing the derivative of loss function w.r.t. to each intermediate node using chain rule.

Example of such graph:

<img src=diagram-20230329.svg width=800>

Important steps in implementation:

```
OpType: a+b, a-b, a*b, a/b, -a, sqrt(a), a^2, smoothstep(a)
struct D: differentiable
struct Dp: pointer to D
graph: computational graph
```

Memory consumption: 700 MB per 10 000 000 nodes in graph. Problems found:
1) BFS backward traversal results in incorrect gradient (fixed, traverse in order of appending to graph)
2) const values are deleted on backward and result in dangling references from intermediate nodes (not fixed, need to create const values on each forward pass)

## SDF edge sampling

$$ \frac{d I_{ijc}}{dr} = \int_{p\in \Gamma} \langle n, \frac{dp}{dr} \rangle \Delta f d\Gamma $$

$$ \Delta f = \lim_{\epsilon\to 0-} f(p+\epsilon n) - \lim_{\epsilon\to 0+} f(p+\epsilon n) $$

$\Delta f$ can be assumed equal to 1.

$$ p = C + (r \cos\alpha, r \sin\alpha) $$

$$ \frac{dp}{dr} = (\cos\alpha, \sin\alpha) $$

$$ n = (\cos\alpha, \sin\alpha) - \text{normal pointing to external} $$

$$ \langle n, \frac{dp}{dr} \rangle = 1 $$

## Complex SDF

Shapes look like https://iquilezles.org/articles/distfunctions2d/, but formulas was derived by me.

5-angled Star

$$ rs - \text{small radius}, rb - \text{big radius} $$

$$ s = (0, rs), b_0 = (0, rb), b = b_0.rot(\alpha) = (rb(-\sin\alpha), rb\cos\alpha) $$

$$ SDF = -\frac{cross(p-s, b-s)}{len(b-s)} $$

**Derivative of SDF with respect to small radius**

$$ \frac{dSDF}{drs} = -\frac{cross'*len - cross*len'}{len^2} $$

$$ cross(p-s, b-s) = (p-s)_x(b-s)_y - (p-s)_y(b-s)_x $$

$$ dcross / drs = (p-s)_x(-1) - (-1)(b-s)_x = -(p-s)_x + (b-s)_x $$

$$ len(b-s) = \sqrt{ (b_x-s_x)^2 + (b_y-s_y)^2 } $$

$$ dlen(b-s) / drs = \frac{1}{2len} (b-s)_y * (0, -1) $$

Now put $dcross$ and $dlen$ into $dSDF/drs$ and get the derivative for small radius!

**Derivative with respect to big radius**

$$ b_0 = (0, rb) \quad drb = (0, 1) $$

$$ b = (rb(-\sin\alpha), rb\cos\alpha) \quad drb = (-\sin\alpha, \cos\alpha) $$

Follow chain rule to compute $db_x/drb$ using derivatives of $b_x$ by $b_{0x}$ AND $b_{0y}$:

$$ \frac{db_x}{drb} = \frac{db_x}{db_{0x}} \frac{db_{0x}}{drb} + \frac{db_x}{db_{0y}} \frac{db_{0y}}{drb} $$

$$ dcross / db_x = -(p-s)_y $$

$$ dcross / db_y = (p-s)_x $$

$$ dcross / drb = -(p-s)_y(-\sin\alpha) + (p-s)_x\cos\alpha $$

$$ dlen / drb = \frac{1}{2len} (2(b-s)_x\frac{db_x}{drb} + 2(b-s)_y\frac{db_y}{drb} ) $$

As with small radius, put $dcross$ and $dlen$ into $dSDF/drb$.


## Files
* Rust sources are src/main.rs, src/loadsdf.rs, src/taskdif.rs and src/taskdif/*
* Cargo.toml - dependecies
* resources/ - precomputed reference images, textures

## Implementation
Written from scratch in Rust. Implements basic 2D geometry: 
* Point struct (which serves as a vector), implemented dot product, cross product, projection of one vector onto another, rotation
* Shape struct for all shapes - analytical, bitmap SDF. Scaling and translation is not implemented. For storage and diff of vectors Point and nalgebra::Vector3 package are used. For image storage ndarray::Array3 package is used.
* Shape::backward sums dMSE/dS for parameters of shape, Shape::step applies Adam or other optimizer
* TriangleMesh defines triangle positions and texture (uv) coordinates. ::render is used for both forward and backward passes, for differentiation on the edge of triangles. ::step applies the optimizer to positions and texture pixels.

## Compile and run
Project was developed on Ubuntu OS, using Rust compiler rustc 1.67.1. Compile and run:
```
cargo build --release && ./target/release/diff2d
```
Cargo should install Rust dependencies. Rendered images will be saved to working directory. Implemented scenes:
* SDF diff
* Bitmap SDF diff, ACDC logo
* Mesh texture and edge sampling
* SDF edge sampling
* Complex SDF
* Auto diff, simple scene
* Multiple image SDF, Portal logo

To create animation:
```
convert -delay 4 -loop 0 *.jpg myimage.gif # imagemagick

convert *.jpg -coalesce -duplicate 1,-2-1 -quiet -layers OptimizePlus  -loop 0 patrol_cycle.gif

ffmpeg -framerate 10 -f image2 -pattern_type glob -i '*.jpg' -c:v libx264 -preset veryslow -crf 10 -b:v 400K -pix_fmt yuv420p output.mp4 

ffmpeg -framerate 20 -pattern_type glob -i '*.jpg' -c:v libwebp -loop 0  -q:v 80 output.webp
```

## Mesh edge sampling performance

Due to unstable optimization, very large textures harm loss and decrease runtime. Loss is final MSE, lower is better.

| optimized texture size  | loss | time, ms/it|
| - | - | - |
| 32x32  | 0.0049 | 9.6 |
| 128x128  | 0.0046 | 10 |
| 1024x1024  | 0.0055 | 45 |

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

ndarray
https://github.com/rust-ndarray/ndarray/blob/master/README-quick-start.md

ndarray sqrt
https://docs.rs/ndarray/0.13.0/ndarray/doc/ndarray_for_numpy_users/index.html
