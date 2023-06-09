#![allow(dead_code)]
#![allow(warnings)]
mod loadsdf;
mod libdif;
use std::{ops::{self, Deref}, f32::consts::PI, collections::{VecDeque, HashMap}, time::Instant, rc::Rc, cell::RefCell, borrow::Borrow, fs};

use image::{GrayImage, GenericImageView};
use libdif::{DGraph, Graph, Dp, GRAPH};

type df32 = Dp;

#[derive(Default, Debug, Clone, Copy)]
struct Point {
    x: df32,
    y: df32
}

impl Point {
    fn const_f(x: f32, y: f32) -> Point {
        return Point{
            x: Dp::const_f(x), 
            y: Dp::const_f(y)
        };
    }

    fn var_f(x: f32, y: f32) -> Point {
        return Point{
            x: Dp::var_f(x), 
            y: Dp::var_f(y)
        };
    }

    fn len(&self) -> df32 {
        return (self.x * self.x + self.y * self.y).sqrt();
    }

    fn dot(&self, b: &Point) -> df32 {
        return self.x * b.x + self.y * b.y;
    }

    fn cross(&self, b: &Point) -> df32 {
        return self.x * b.y - b.x * self.y;
    }

    fn proj(&self, b: &Point) -> df32 {
        return self.dot(b) / self.len();
    }

    fn rotate(&self, alpha: f32) -> Point {
        let acos = Dp::const_f(alpha.cos());
        let asin = Dp::const_f(alpha.sin());
        return Point {
            x: self.x * acos - self.y * asin,
            y: self.x * asin + self.y * acos // TODO derive wrt angle
        };
    }
}

impl ops::Add<&Point> for &Point {
    type Output = Point;

    fn add(self, _rhs: &Point) -> Point {
        //println!("> Point.add(Point) was called");
        return Point {x: self.x + _rhs.x, y: self.y + _rhs.y};
    }
}

impl ops::Mul<f32> for &Point {
    type Output = Point;

    fn mul(self, _rhs: f32) -> Point {
        let a = Dp::const_f(_rhs);
        return Point {x: self.x * a, y: self.y * a};
    }
}

impl ops::Neg for &Point {
    type Output = Point;

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

// #[derive(Default)]
// struct TriangleMesh {
//     verticies: Vec<Point>,
//     indices: Vec<u32>, // triangles
//     colors: Vec<ShapeColor>, // for each face
//     tc: Vec<Point>,
//     texture: GrayImage,
// }

// impl TriangleMesh {
//     fn render(&self, point: Point, color: &mut ShapeColor) -> bool {
//         *color = ShapeColor { r: 1., g: 1., b: 1., layer: 0 };
//         return false;
//     }
// }

// fn smin(a: f32, b: f32, k: f32) -> f32 {
//     let h = (1. - (a - b).abs() / k).max(0.);
//     let m = h * h * h / 2.;
//     let s = m * k / 3.;
//     return a.min(b) - s;
// }

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


fn circle_sdf(c: Point, r: df32, p: Point) -> df32 {
    return (p - c).len() - r;
}


pub fn task_autodiff(save_path: &str) {
    println!("\n=== Scene: autodiff");
    let mut start_time = Instant::now();

    let yellow = ShapeColor::new(255, 220, 3);
    let black = ShapeColor{r: 0., g: 0., b: 0., layer: 0};
    
    println!("Initialization took {:?}", start_time.elapsed());
    start_time = Instant::now();

    let imgx = 128;
    let imgy = imgx;

    let circ_c: Point = Point::var_f(0.2, 0.2);
    let circ_r = Dp::var_f(0.2);

    let circ2_c: Point = Point::var_f(-0.1, -0.2);
    let circ2_r = Dp::var_f(0.3);
    let th = 0.1;

    fs::create_dir("anim");

    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);

    // ******************** FORWARD PASS

    // Iterate over the coordinates and pixels of the image
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        // if x % 80 == 0 && y % 80 == 0 {
        //     println!("x {} y {}", x, y);
        //     {
        //         let graph = GRAPH();
        //         let objects = &graph.as_ref().borrow().objects;
        //         println!("Object count {} free count {}", objects.len(), graph.as_ref().borrow().free_cnt());
        //     }
        // }
        let xf = x as f32 / imgx as f32 - 0.5;
        let yf = y as f32 / imgy as f32 - 0.5;
        //let r = (xf / imgx as f32 * 128.0) as u8;
        //let b = (yf / imgy as f32 * 128.0) as u8;

        let p: Point = Point::const_f(xf, yf);
        let sdf = circle_sdf(circ_c, circ_r, p);
        let smstep = (-sdf).smoothstep(0., th);

        let sdf2 = circle_sdf(circ2_c, circ2_r, p);
        let smstep2 = (-sdf2).smoothstep(0., th);

        let w = smstep;
        let w2 = (Dp::const_f(1.) - w) * smstep2;

        let mut out_color: ShapeColor;
        {
            let graph = GRAPH();
            let objects = &graph.as_ref().borrow().objects;
            out_color = yellow * (objects[w.id].scalar + objects[w2.id].scalar);
        }
        w.backward();
        w2.backward(); // delete nodes

        out_color = out_color * 255.;
        *pixel = image::Rgb([
            out_color.r as u8,
            out_color.g as u8,
            out_color.b as u8]);
    }

    println!("Rendering took {:?}", start_time.elapsed());
    start_time = Instant::now();

    // Save the image as “fractal.png”, the format is deduced from the path
    imgbuf.save("anim/reference.jpg").unwrap();

    println!("Saving took {:?}", start_time.elapsed());

    {
        let graph = GRAPH();
        let objects = &graph.as_ref().borrow().objects;
        println!("Object count {} free count {}", objects.len(), graph.as_ref().borrow().free_cnt());
    }
    // return;

    // ******************** BACKWARD PASS

    let refimg = loadsdf::loadimage("anim/reference.jpg");

    let mut circ_c: Point = Point::var_f(0.2, 0.);
    let mut circ_r = Dp::var_f(0.2);

    let mut circ2_c: Point = Point::var_f(-0.1, 0.);
    let mut circ2_r = Dp::var_f(0.1);

    for _it in 0..20 {
        //let mut mse = Dp::const_f(0.);
        let mut mse: f32 = 0.;
        for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
            let xf = x as f32 / imgx as f32 - 0.5;
            let yf = y as f32 / imgy as f32 - 0.5;
            //let r = (xf / imgx as f32 * 128.0) as u8;
            //let b = (yf / imgy as f32 * 128.0) as u8;

            // again forward pass
            let p: Point = Point::const_f(xf, yf);
            let sdf = circle_sdf(circ_c, circ_r, p);
            let smstep1 = (-sdf).smoothstep(0., th);

            let sdf2 = circle_sdf(circ2_c, circ2_r, p);
            let smstep2 = (-sdf2).smoothstep(0., th);

            let w1 = smstep1;
            let w2 = (Dp::const_f(1.) - w1) * smstep2;
            let w = w1 + w2;

            // let mut out_color: ShapeColor;
            // {
            //     let graph = GRAPH();
            //     let objects = &graph.as_ref().borrow().objects;
            //     out_color = yellow * (objects[w.id].scalar + objects[w2.id].scalar);
            // }

            // mse, derivative
            let pixref = refimg.get_pixel(x, y);
            let dataref = pixref.0;
            let rd = Dp::const_f((dataref[0] as f32) / 255.);
            let gd = Dp::const_f((dataref[1] as f32) / 255.);
            let bd = Dp::const_f((dataref[2] as f32) / 255.);

            let yrd = Dp::const_f(yellow.r);
            let ygd = Dp::const_f(yellow.g);
            let ybd = Dp::const_f(yellow.b);

            let pixmse_d = ((w * yrd - rd).square()
                + (w * ygd - gd).square()
                + (w * ybd - bd).square()
            );

            {
                let graph = GRAPH();
                let freeze;
                {
                    let objects = &graph.as_ref().borrow().objects;

                    if objects[sdf.id].scalar > -th {
                        freeze = false;
                    } else {
                        freeze = true;
                    }
                }
                graph.borrow_mut().freeze = freeze;
            }
            pixmse_d.backward();

            // mse = mse + pixmse_d;
            

            // dr += -pixmse_d * smoothstep_d(0., th, -r) * (-1.) * (-1.);

            let mut out_color: ShapeColor;
            {
                let graph = GRAPH();
                let objects = &graph.as_ref().borrow().objects;
                out_color = yellow * objects[w.id].scalar;
                mse += objects[pixmse_d.id].scalar;
            }

            out_color = out_color * 255.;
            *pixel = image::Rgb([
                out_color.r as u8,
                out_color.g as u8,
                out_color.b as u8]);
        }
        // let imgx_d = Dp::const_f(imgx as f32);
        // let imgy_d = Dp::const_f(imgy as f32);
        // mse = mse / (imgx_d * imgy_d);
        // mse.backward();
        mse = mse / imgx as f32  / imgy as f32;


        let new_r;
        let new_c: (f32, f32);

        let new2_r;
        let new2_c: (f32, f32);
        {
            let graph = GRAPH();
            let objects = &graph.as_ref().borrow().objects;
            println!("Object count {} free count {}", objects.len(), graph.as_ref().borrow().free_cnt());

            let dr: f32 = objects[circ_r.id].scalar_d;
            println!("mse={:.3} circ_r={:.3} dr={:.3}", 
                mse, 
                objects[circ_r.id].scalar,
                dr
            );

            let step = |id: usize, lr: f32| -> f32 {
                let deriv = objects[id].scalar_d;
                return objects[id].scalar - deriv * lr;
            };

            // new_r = objects[circ_r.id].scalar - dr * 0.0000001;
            // new_C = objects[circ_c.x.id].scalar - dr * 0.0000001;
            let lr = 0.000001;
            new_r = step(circ_r.id, lr);
            new_c = (step(circ_c.x.id, lr), step(circ_c.y.id, lr));

            new2_r = step(circ2_r.id, lr);
            new2_c = (step(circ2_c.x.id, lr), step(circ2_c.y.id, lr));
        }
        circ_r = Dp::var_f(new_r);
        circ_c = Point::var_f(new_c.0, new_c.1);

        circ2_r = Dp::var_f(new2_r);
        circ2_c = Point::var_f(new2_c.0, new2_c.1);
        
        
        let filename: String = format!("anim/{}{:0>2}.jpg", save_path, _it);
        imgbuf.save(filename).unwrap();

    }
}




fn test_autodif() {
    let x = Dp::var_f(2.);
    let y = Dp::var_f(3.);

    // *** Test 1
    // let mut z = (x + y) * y * x;

    // *** Test 2
    // let mut z = (x / y).sqrt() - y;

    // *** Test 3
    // let z = x.square();
    // let z1 = -z;
    // let z2 = z1*z;
    // z2.backward();
    // {
    //     let graph = GRAPH();
    //     let g = graph.as_ref().borrow();
    //     let objects = &g.objects;
    //     println!("Graph length {:?} free {:?}", objects.len(), g.free_cnt());
    //     println!("z2={}", objects[z2.id].scalar);
    //     println!("dx={}", objects[x.id].scalar_d);
    //     println!("dy={}", objects[y.id].scalar_d);
    // }

    // *** Test 4
    let z = -x;
    let z1 = (z+(z+(z+z)));
    z1.backward();
    {
        let graph = GRAPH();
        let g = graph.as_ref().borrow();
        let objects = &g.objects;
        println!("Graph length {:?} free {:?}", objects.len(), g.free_cnt());
        println!("z2={}", objects[z1.id].scalar);
        println!("dx={}", objects[x.id].scalar_d);
        println!("dy={}", objects[y.id].scalar_d);
    }
    

    // *** Test smoothstep
    // let x1 = Dp::none_float(-0.1);
    // let z1 = x1.smoothstep(0., 1.);
    // let x2 = Dp::none_float(0.3);
    // let z2 = x2.smoothstep(0., 1.);
    // let x3 = Dp::none_float(1.1);
    // let z3 = x3.smoothstep(0., 1.);
    // z1.backward();
    // z2.backward();
    // z3.backward();
    // {
    //     let graph = GRAPH();
    //     let objects = &graph.as_ref().borrow().objects;
    //     println!("Graph length {:?}", objects.len());
    //     println!("z1={} z2={} z3={}", objects[z1.id].scalar, objects[z2.id].scalar, objects[z3.id].scalar);
    //     println!("dx1={} dx2={} dx3={}", objects[x1.id].scalar_d, objects[x2.id].scalar_d, objects[x3.id].scalar_d);
    // }


    // *** Performance test
    // for i in 0..100 {
    //     z = z + x;
    // }
    // z.backward();
    

    // *** Test node cleanup
    // {
    //     let graph = GRAPH();
    //     let g = graph.as_ref().borrow();
    //     let objects = &g.objects;
    //     println!("Graph length {:?} free {:?}", objects.len(), g.free_cnt());
    //     println!("z={}", objects[z.id].scalar);
    //     println!("dx={}", objects[x.id].scalar_d);
    //     println!("dy={}", objects[y.id].scalar_d);
    // }
    // let mut z = (x / y).sqrt() - y;
    // z.backward();
    // {
    //     let graph = GRAPH();
    //     let g = graph.as_ref().borrow();
    //     let objects = &g.objects;
    //     println!("Graph length {:?} free {:?}", objects.len(), g.free_cnt());
    //     println!("z={}", objects[z.id].scalar);
    //     println!("dx={}", objects[x.id].scalar_d);
    //     println!("dy={}", objects[y.id].scalar_d);
    // }
}


fn main() {
    task_autodiff("fractal");
    // test_autodif();
    
    return;
}