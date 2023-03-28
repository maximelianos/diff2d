use std::{ops, f32::consts::PI, collections::VecDeque, time::Instant};

use grid::Grid;
use image::GrayImage;
mod loadsdf;

#[derive(Default, Copy, Clone, Debug, PartialEq)]
struct Point {
    x: f32,
    y: f32
}

impl Point {
    fn new(x: f32, y: f32) -> Point {
        return Point{x: x, y: y};
    }

    fn len(&self) -> f32 {
        return (self.x * self.x + self.y * self.y).sqrt();
    }

    fn dot(&self, b: Point) -> f32 {
        return self.x * b.x + self.y * b.y;
    }

    fn cross(&self, b: Point) -> f32 {
        return self.x * b.y - b.x * self.y;
    }

    fn proj(&self, b: Point) -> f32 {
        return self.dot(b) / self.len();
    }

    fn rotate(&self, alpha: f32) -> Point {
        return Point::new(
            self.x * alpha.cos() - self.y * alpha.sin(),
            self.x * alpha.sin() + self.y * alpha.cos());
    }
}

impl ops::Add<Point> for Point {
    type Output = Point;

    fn add(self, _rhs: Point) -> Point {
        //println!("> Point.add(Point) was called");
        return Point {x: self.x + _rhs.x, y: self.y + _rhs.y};
    }
}

impl ops::Mul<f32> for Point {
    type Output = Point;

    fn mul(self, _rhs: f32) -> Point {
        return Point {x: self.x * _rhs, y: self.y * _rhs};
    }
}

impl ops::Neg for Point {
    type Output = Self;

    fn neg(self) -> Point {
        //println!("> Point.neg() was called");
        return Point{x: -self.x, y: -self.y};
    }
}

impl ops::Sub<Point> for Point {
    type Output = Point;

    fn sub(self, _rhs: Point) -> Point {
        //println!("> Point.sub(Point) was called");
        return Point {x: self.x - _rhs.x, y: self.y - _rhs.y};
    }
}

struct ShapeTransform {
    scale: f32,
    rotate: f32, // in radians
    translate: Point
}

impl Default for ShapeTransform {
    fn default() -> Self { Self{scale: 1., rotate: 0., translate: Point::new(0., 0.)} }
}

impl ShapeTransform {
    fn apply(&self, p: Point) -> Point {
        // apply same transform to translation vector
        let translate = self.translate; //(self.translate * (1./self.scale)).rotate(-self.rotate);
        // transform in reverse order
        let mut q = p - translate;
        q = q.rotate(-self.rotate);
        q = q * (1./self.scale);
        return q;
    }
}

#[derive(Default, Copy, Clone, Debug, PartialEq)]
enum ShapeType {
    #[default]
    Circle,
    
    Rectangle,
    Triangle,
    LineSegment,
    Bitmap
}

#[derive(Default)]
struct Shape {
    stype: ShapeType,
    // circle requires C
    r: f32,
    // rect
    w: f32,
    h: f32,
    // tri
    A: Point,
    B: Point,
    C: Point,
    // segment requires A, B
    input_bitmap: GrayImage,
    bitmap: GrayImage, // can be None or GrayImage
    tosdf_scale: f32, // less makes smoother edges in binary image converted to SDF, doesn't affect absolute SDF values
    field_scale: f32, // for bitmaps: d *= field_scale
    color: ShapeColor,
    gradient_color: ShapeColor,
    gradient_p1: Point,
    gradient_p2: Point,

    local_transform: ShapeTransform,
    transform: ShapeTransform,
    layer: usize
}

#[derive(Default, Copy, Clone, Debug, PartialEq)]
struct ShapeColor {
    r: f32,
    g: f32,
    b: f32,
    layer: i32,
}

impl ShapeColor {
    fn new(r: i32, g: i32, b: i32) -> ShapeColor {
        ShapeColor { r: r as f32 / 255., g: g as f32 / 255., b: b as f32 / 255., layer: 0 }
    }
}

impl ops::Add<ShapeColor> for ShapeColor {
    type Output = ShapeColor;

    fn add(self, _rhs: ShapeColor) -> ShapeColor {
        //println!("> Point.add(Point) was called");
        return ShapeColor {
            r: self.r + _rhs.r,
            g: self.g + _rhs.g,
            b: self.b + _rhs.b,
            layer: self.layer };
    }
}

impl ops::Mul<f32> for ShapeColor {
    type Output = ShapeColor;

    fn mul(self, rhs: f32) -> Self::Output {
        return ShapeColor {
            r: self.r * rhs,
            g: self.g * rhs,
            b: self.b * rhs,
            layer: self.layer };
    }
}

impl Shape {
    fn build(&mut self) {
        use ShapeType::*;
        match self.stype {
            Circle => (),
            Rectangle => (),
            Triangle => {
                // order points in clockwise manner
                let a = self.B - self.A;
                let b = self.C - self.A;
                //println!("> Cross product: {}", a.cross(b));
                if a.cross(b) > 0. {
                    let t = self.B;
                    self.B = self.C;
                    self.C = t;
                }

            },
            LineSegment => (),
            Bitmap => ()
        }
    }

    fn distance_pre(&self, point: Point) -> f32 {
        use ShapeType::*;
        // transform is inverse: global space -> local space -> object space
        let p = self.local_transform.apply(self.transform.apply(point));
        match self.stype {
            Circle => {
                return (p - self.C).len() - self.r;
            },
            Rectangle => {
                let w2 = self.w / 2.0;
                let h2 = self.h / 2.0;
                let d1 = p.y.abs() - h2;
                let d2 = p.x.abs() - w2;
                if d1 < 0. && d2 < 0. {
                    return d1.max(d2); // closest border
                } else if -w2 < p.x && p.x < w2 {
                    return d1;
                } else if -h2 < p.y && p.y < h2 {
                    return d2;
                } else {
                    return Point::new(p.x.abs() - w2, p.y.abs() - h2).len();
                }
            },
            Triangle => {
                let a = self.B - self.A;
                let b = self.C - self.B;
                let c = self.A - self.C;
                let pa = p - self.A;
                let pb = p - self.B;
                let pc = p - self.C;
                // on which side of triangle is our point
                let ca = a.cross(pa);
                let cb = b.cross(pb);
                let cc = c.cross(pc);
                let pja = a.proj(pa);
                let pjb = b.proj(pb);
                let pjc = c.proj(pc);
                // case A
                if ca >= 0. && 0. <= pja && pja <= a.len() {
                    return ca / a.len();
                } else if cb >= 0. && 0. <= pjb && pjb <= b.len() {
                    return cb / b.len();
                } else if cc >= 0. && 0. <= pjc && pjc <= c.len() {
                    return cc / c.len();
                // case B
                } else if cb > 0. && pjb < 0. || ca > 0. && !(pja < 0.)  {
                    // a-b side boundary
                    return pb.len();
                } else if cc > 0. && pjc < 0. || cb > 0. && !(pjb < 0.)  {
                    // b-c side boundary
                    return pc.len();
                } else if ca > 0. && pja < 0. || cc > 0. && !(pjc < 0.)  {
                    // c-a side boundary
                    return pa.len();
                } else {
                    // inside triangle. cross product is negative thus max is needed
                    return f32::max(ca / a.len(), 
                        f32::max(cb / b.len(), cc / c.len()));
                }
            },
            LineSegment => {
                let a: Point = self.B - self.A;
                let pa: Point = p - self.A;
                let proj: f32 = a.proj(pa);
                if 0. <= proj && proj <= a.len() {
                    return a.cross(pa).abs() / a.len();
                } else {
                    return f32::min(pa.len(), (p - self.B).len());
                }
            },
            Bitmap => {
                let (w, h): (u32, u32) = self.bitmap.dimensions();
                let q = Point::new(
                    p.x * w as f32, // TODO: think about scaling bitmap
                    p.y * h as f32
                );
                
                // limits can be inaccurate, there are checks after this
                // if q.x < -0. || q.x > w as f32 - 1. || q.y < -0. || q.y > h as f32 - 1. {
                //     return 1000.;
                // }

                let scale: f32;
                if self.field_scale > 0. {
                    scale = self.field_scale;
                } else {
                    scale = 1.;
                }
                let dist_large: f32 = 1000.;

                let x = q.x.floor() as i32;
                let y = q.y.floor() as i32;

                let get_sdf = |px: i32, py: i32| -> f32 {
                    if 0 <= px && px < w as i32 && 0 <= py && py < h as i32 {
                        0.5 - self.bitmap.get_pixel(px as u32, py as u32).0[0] as f32 / 255. // TODO: very ugly, please fix
                    } else {
                        dist_large
                    }
                };
                let v1 = get_sdf(x, y);
                let v2 = get_sdf(x, y+1);
                let v3 = get_sdf(x+1, y+1);
                let v4 = get_sdf(x+1, y);
                
                let kx = q.x - x as f32;
                let ky = q.y - y as f32;

                let k1 = (1.-kx)*(1.-ky);
                let k2 = (1.-kx)*(ky);
                let k3 = (kx)*(ky);
                let k4 = (kx)*(1.-ky);
                return (k1*v1 + k2*v2 + k3*v3 + k4*v4)/4.*scale; // TODO: this is not true distance, please think about it
            }
        }
    }

    fn distance(&self, point: Point) -> f32 {
        // in case of scaling, we must multiply the sdf by scaling factor
        return self.distance_pre(point) * self.local_transform.scale * self.transform.scale;
    }

    fn gradient(&self, point: Point) -> ShapeColor {
        let p = self.local_transform.apply(self.transform.apply(point));
        let gp = self.gradient_p2 - self.gradient_p1;
        let mut k = 1. - gp.proj(p - self.gradient_p1) / gp.len();
        k = k.max(0.).min(1.);
        return self.color * k + self.gradient_color * (1. - k);
    }

    fn bitmap_to_sdf(&mut self, sdf_resolution: f32, sdf_padding: f32) {
        // sdf_resolution is a small number inside [0, 1], makes the SDF matrix small
        // sdf_padding is fraction of padding on all sides of image, to compute accurate SDF
        println!("Convert binary image to SDF");
        let (w, h): (u32, u32) = self.input_bitmap.dimensions();
        
        let w = (w as f32 * (1.+sdf_padding*2.)) as u32;
        let h = (h as f32 * (1.+sdf_padding*2.)) as u32;
        let mut inp_bitmap = image::ImageBuffer::new(w, h);
        let y_off = (h as f32 * sdf_padding / 2.) as i64;
        let x_off = (w as f32 * sdf_padding / 2.) as i64;
        image::imageops::overlay(&mut inp_bitmap, &self.input_bitmap, x_off, y_off);
        self.input_bitmap = inp_bitmap;
        

        println!("Input dimensions: {w}x{h}");
        let outh = (h as f32 * sdf_resolution) as u32;
        let outw = (w as f32 * sdf_resolution) as u32;
        println!("SDF dimensions: {outw}x{outh}");

        let get_pixel = |px: u32, py: u32| -> f32 {
            self.input_bitmap.get_pixel(px, py).0[0] as f32 / 255. // TODO: very ugly, please fix
        };

        let is_inside = |px: u32, py: u32| -> bool {
            get_pixel(px, py) > 0.5
        };

        let mut visited: Grid<u8> = Grid::new(h as usize, w as usize);
        let mut parent: Grid<(usize, usize)> = Grid::new(h as usize, w as usize);
        let mut result: Grid<f32> = Grid::new(h as usize, w as usize);
        let mut deque: VecDeque<(usize, usize)> = VecDeque::new();

        for y in 0..h {
            for x in 0..w {
                if is_inside(x, y) {
                    visited[y as usize][x as usize] = 1;
                    deque.push_back((x as usize, y as usize));
                }
            }
        }
        while !deque.is_empty() {
            let (x, y) = deque.pop_front().unwrap();
            // ************************** BFS
            let px = x as i32;
            let py = y as i32;
            let ux = x;
            let uy = y;
            if visited[uy][ux] == 2 {
                continue;
            }
            
            // println!("Pix at x={ux} y={uy}: {pix}");
            let mut min_dist: f32 = 1000.;
            let mut min_parent: (usize, usize) = (0, 0);
            if is_inside(ux as u32, uy as u32) {
                result[uy][ux] = -1.;
                parent[uy][ux] = (ux, uy);
                min_dist = -1.;
            }
            let cur_point = Point::new(ux as f32, uy as f32);
            for i in -2..3 {
                for j in -2..3 {
                    let nx = px + j;
                    let ny = py + i;
                    if 0 <= nx && nx < w as i32 && 0 <= ny && ny < h as i32 {
                        if visited[ny as usize][nx as usize] == 0 {
                            if i==0&&(j==-1||j==1) || j==0&&(i==-1||i==1) {
                                visited[ny as usize][nx as usize] = 1;
                                deque.push_back((nx as usize, ny as usize));
                            }
                        } else if visited[ny as usize][nx as usize] == 2 {
                            let (pnx, pny) = parent[ny as usize][nx as usize];
                            let par_point = Point::new(pnx as f32, pny as f32);
                            if (par_point - cur_point).len() < min_dist {
                                min_dist = (par_point - cur_point).len();
                                min_parent = (pnx, pny);
                                // println!("Distance between {:?} {:?} = {:?}", cur_point, par_point, min_dist);
                            }
                        }
                    }
                }
            }

            if !is_inside(ux as u32, uy as u32) {
                result[uy][ux] = min_dist;
                parent[uy][ux] = min_parent;
            }

            visited[uy][ux] = 2;
        }

        // ******************* Now propagate inward

        for y in 0..h {
            for x in 0..w {
                visited[y as usize][x as usize] = 0;
                if !is_inside(x, y) {
                    visited[y as usize][x as usize] = 1;
                    deque.push_back((x as usize, y as usize));
                }
            }
        }

        while !deque.is_empty() {
            let (x, y) = deque.pop_front().unwrap();
            // ************************** BFS
            let px = x as i32;
            let py = y as i32;
            let ux = x;
            let uy = y;
            if visited[uy][ux] == 2 {
                continue;
            }
            
            let mut min_dist: f32 = 1000.;
            let mut min_parent: (usize, usize) = (0, 0);
            if !is_inside(ux as u32, uy as u32) {
                parent[uy][ux] = (ux, uy);
            }
            let cur_point = Point::new(ux as f32, uy as f32);
            for i in -2..3 {
                for j in -2..3 {
                    let nx = px + j;
                    let ny = py + i;
                    if 0 <= nx && nx < w as i32 && 0 <= ny && ny < h as i32 {
                        if visited[ny as usize][nx as usize] == 0 {
                            if i==0&&(j==-1||j==1) || j==0&&(i==-1||i==1) {
                                visited[ny as usize][nx as usize] = 1;
                                deque.push_back((nx as usize, ny as usize));
                            }
                        } else if visited[ny as usize][nx as usize] == 2 {
                            let (pnx, pny) = parent[ny as usize][nx as usize];
                            let par_point = Point::new(pnx as f32, pny as f32);
                            if (par_point - cur_point).len() < min_dist {
                                min_dist = (par_point - cur_point).len();
                                min_parent = (pnx, pny);
                                // println!("Distance between {:?} {:?} = {:?}", cur_point, par_point, min_dist);
                            }
                        }
                    }
                }
            }

            let small_distance: f32 = 1.;
            if is_inside(ux as u32, uy as u32) {
                result[uy][ux] = -min_dist + small_distance; // negative, because inside of shape
                parent[uy][ux] = min_parent;
            }

            visited[uy][ux] = 2;
        }

        let sdf_scale = self.tosdf_scale;

        let mut imgbuf = image::ImageBuffer::new(outw, outh);
        for y in 0..outh {
            for x in 0..outw {
                let y_input = (y as f32 / sdf_resolution) as u32;
                let x_input = (x as f32 / sdf_resolution) as f32;
                let val = result[y_input as usize][x_input as usize];
                // println!("Result at x={x} y={y}: {val:.3} parent={:?}", parent[y as usize][x as usize]);
                let res: f32;
                res = (0.5 + -val * sdf_scale).max(0.).min(1.);

                let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
                *pixel = image::Luma([(res * 255.) as u8]);
            }
        }

        // imgbuf.save("tosdf.png").unwrap();

        self.bitmap = imgbuf;
    }
}

#[derive(Default)]
struct TriangleMesh {
    verticies: Vec<Point>,
    indices: Vec<u32>, // triangles
    colors: Vec<ShapeColor>, // for each face
    tc: Vec<Point>,
    texture: GrayImage,
}

impl TriangleMesh {
    fn render(&self, point: Point, color: &mut ShapeColor) -> bool {
        *color = ShapeColor { r: 1., g: 1., b: 1., layer: 0 };
        return false;
    }
}

fn smin(a: f32, b: f32, k: f32) -> f32 {
    let h = (1. - (a - b).abs() / k).max(0.);
    let m = h * h * h / 2.;
    let s = m * k / 3.;
    return a.min(b) - s;
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum SceneType {
    All,
    Letters,
    Bird,
    Flame,
    Fox,
    Shapes,
    Art
}

fn smoothstep(left: f32, right: f32, x: f32) -> f32 {
    if x < left {
        0.
    } else if x > right {
        1.
    } else {
        let y = (x - left) / (right - left);
        y * y * (3. - 2. * y)
    }
}

fn main_letters(scene: SceneType, save_path: &str) {
    println!("\n=== Scene {:?}", scene);
    let mut start_time = Instant::now();

    let yellow = ShapeColor::new(255, 220, 3);
    let orange = ShapeColor::new(251, 133, 0);
    let cyan = ShapeColor::new(0, 128, 128);
    let pink = ShapeColor::new(255, 0, 102);
    let red = ShapeColor::new(204, 0, 0);
    let sky = ShapeColor::new(142, 202, 230);
    let azul = ShapeColor::new(0, 102, 230);
    let grey = ShapeColor::new(40, 42, 42);
    let leaf = ShapeColor::new(0, 153, 0);
    let grass = ShapeColor::new(128, 255, 128);
    let fur = ShapeColor::new(230, 230, 230);

    let main_aliasing_scale = 0.001;
    let main_border_scale = 0.004;

    let scale = 0.08;


    // ********************** SHAPES
    let shapes_scale = scale; //scale;
    let shapes_pos = Point::new(-scale*0.5, -scale*0.5);
    
    let mut s_circ = Shape {
        stype: ShapeType::Circle,
        r: 1.,
        transform: ShapeTransform {
            scale: shapes_scale, 
            rotate: 0.,
            translate: Point::new(0., 0.) // shapes_pos
        },
        color: yellow,
        ..Default::default()
    };
    s_circ.build();

    // let mut s_tri = Shape {
    //     stype: ShapeType::Triangle,
    //     A: Point::new(-4., 4.),
    //     B: Point::new(-5., -3.),
    //     C: Point::new(-1., -2.),
    //     transform: ShapeTransform {
    //         scale: shapes_scale, 
    //         rotate: 0.,
    //         translate: shapes_pos
    //     },
    //     color: grass,
    //     ..Default::default()
    // };
    // s_tri.build();

    // let mut s_rect = Shape {
    //     stype: ShapeType::Rectangle,
    //     w: 1.,
    //     h: 1.5,
    //     transform: ShapeTransform {
    //         scale: shapes_scale, 
    //         rotate: 0.,
    //         translate: Point::new(1., 1.) * shapes_scale + shapes_pos 
    //     },
    //     color: red,
    //     ..Default::default()
    // };
    // s_rect.build();

    // let mut mesh: TriangleMesh = Default::default();
    // let mut mesh_color: ShapeColor = Default::default();

    
    let mut shapes: Vec<Shape> = Vec::new();

    shapes.push(s_circ);
    // shapes.push(s_tri);
    // shapes.push(s_rect);

    println!("Shapes N={}", shapes.len());


    let mut dist: Vec<f32> = vec![0.; shapes.len()];

    println!("Initialization took {:?}", start_time.elapsed());
    start_time = Instant::now();

    let imgx = 400;
    let imgy = imgx;

    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);

    // Iterate over the coordinates and pixels of the image
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let xf = x as f32 / imgx as f32 - 0.5;
        let yf = y as f32 / imgy as f32 - 0.5;
        //let r = (xf / imgx as f32 * 128.0) as u8;
        //let b = (yf / imgy as f32 * 128.0) as u8;

        let p: Point = Point::new(xf, yf);
        let mut r: f32 = 1000.;
        let mut min_idx: usize = 0;
        for i in 0..shapes.len() {
            let ri = shapes[i].distance(p);
            if ri < r {
                min_idx = i;
            }
            r = r.min(ri);
            // r = smin(r, ri, scale/3.);
            dist[i] = ri;
        }

        let mut aliasing_scale = main_aliasing_scale;
        let mut border_scale = main_border_scale;

        let background = ShapeColor{r: 0., g: 0., b: 0., layer: 0};
        let mut out_color: ShapeColor = background;

        let mut blend: ShapeColor = ShapeColor::new(0, 0, 0);
        let mut th = scale / 3.;
        let mut w = 0.;

        for i in 0..shapes.len() {
            let ri = dist[i];
            w = (1. - w) * smoothstep(0., th, -ri);
            let cur_color = shapes[i].color;
            blend = blend + cur_color * w;
        }

        out_color = blend;

        // if mesh.render(p, &mut mesh_color) {
        //     out_color = mesh_color;
        // }


        out_color = out_color * 255.;
        *pixel = image::Rgb([
            out_color.r as u8,
            out_color.g as u8,
            out_color.b as u8]);
    }

    println!("Rendering took {:?}", start_time.elapsed());
    start_time = Instant::now();

    // Save the image as “fractal.png”, the format is deduced from the path
    imgbuf.save(save_path).unwrap();

    println!("Saving took {:?}", start_time.elapsed());
}


fn main() {
    main_letters(SceneType::Shapes, "fractal.png");
    
    return;
}