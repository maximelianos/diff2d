#![allow(dead_code)]
#![allow(warnings)]
mod loadsdf;
mod libdif;
mod libgeom;
mod hashmap;
use std::{ops::{self, Deref}, f32::consts::PI, collections::{VecDeque, HashMap}, time::Instant, rc::Rc, cell::RefCell, borrow::Borrow, fs::File};

use grid::Grid;
use image::{GrayImage, GenericImageView};
//use libdif::{DGraph, Graph, Dp, GRAPH};

use libgeom::Point;

// struct ShapeTransform {
//     scale: f32,
//     rotate: f32, // in radians
//     translate: Point
// }

// impl Default for ShapeTransform {
//     fn default() -> Self { Self{scale: 1., rotate: 0., translate: Point::new(0., 0.)} }
// }

// impl ShapeTransform {
//     fn apply(&self, p: Point) -> Point {
//         // apply same transform to translation vector
//         let translate = self.translate; //(self.translate * (1./self.scale)).rotate(-self.rotate);
//         // transform in reverse order
//         let mut q = p - translate;
//         q = q.rotate(-self.rotate);
//         q = q * (1./self.scale);
//         return q;
//     }
// }

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
    // circle
    C: Point,
    r: f32,
    th: f32,
    // tmp forward
    p_len: f32,
    sdf: f32,

    // rect
    // w: df32,
    // h: df32,
    
    dr: f32,
    dC: Point
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
    fn distance_pre(&mut self, point: Point) -> f32 {
        use ShapeType::*;
        // transform is inverse: global space -> local space -> object space
        // let p = self.local_transform.apply(self.transform.apply(point));
        let p = point;
        let sdf;
        match self.stype {
            Circle => {
                self.p_len = (p - self.C).len();
                sdf = self.p_len - self.r;
            },
            // Rectangle => {
            //     let _2 = Dp::const_f(2.0);
            //     let w2 = self.w / _2;
            //     let h2 = self.h / _2;

            //     let d1 = p.y.abs() - h2;
            //     let d2 = p.x.abs() - w2;
            //     if d1.s() < 0. && d2.s() < 0. {
            //         if d1.s() > d2.s() { return d1; } else { return d2; } // closest border
            //     } else if -w2.s() < p.x.s() && p.x.s() < w2.s() {
            //         return d1;
            //     } else if -h2.s() < p.y.s() && p.y.s() < h2.s() {
            //         return d2;
            //     } else {
            //         return Point {
            //             x: p.x.abs() - w2,
            //             y: p.y.abs() - h2
            //         }.len();
            //     }
            // },
            _ => panic!("Unimplemented shape")
        }
        self.sdf = sdf;
        return sdf;
    }

    fn w(&mut self, point: Point) -> f32 {
        let sdf = self.distance_pre(point);
        return smoothstep(0., self.th, -sdf);
    }

    fn backward(&mut self, point: Point, dldw: f32) {
        let dldsdf = dldw * smoothstep_d(0., self.th, -self.sdf) * (-1.);

        self.dr = self.dr - dldsdf;
        if self.p_len > 0.00001 { // avoid div by 0
            self.dC = self.dC + (-point + self.C) * (1. / self.p_len) * dldsdf * 10.;
        }
    }

    fn step(&mut self, lr: f32) {
        use ShapeType::*;
        match self.stype {
            Circle => {
                self.r = self.r - self.dr * lr;
                self.dr = 0.;
                self.C = self.C - self.dC * lr;
                self.dC = Point::new(0., 0.);
            },
            _ => ()
        }
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


fn circle_sdf(c: Point, r: f32, p: Point) -> f32 {
    return (p - c).len() - r;
}

pub fn smoothstep(left: f32, right: f32, x: f32) -> f32 {
    let res: f32;
    if x < left { res = 0.; } 
    else if x > right { res = 1.; }
    else {
        let y = (x - left) / (right - left);
        res = y * y * (3. - 2. * y);
    }

    return res;
}

pub fn smoothstep_d(left: f32, right: f32, x: f32) -> f32 {
    let res: f32;
    if x < left { res = 0.; }
    else if x > right { res = 0.; }
    else {
        let y = (x - left) / (right - left);
        res = 6. * y * (1. - y) / (right - left);
    }

    return res;
}


fn small(scene: SceneType, save_path: &str) {
    println!("\n=== Scene {:?}", scene);
    let mut start_time = Instant::now();

    let yellow = ShapeColor::new(255, 220, 3);
    let black = ShapeColor{r: 0., g: 0., b: 0., layer: 0};
    
    println!("Initialization took {:?}", start_time.elapsed());
    start_time = Instant::now();

    let imgx = 100;
    let imgy = imgx;

    let th = 0.2;
    let mut circ = Shape {
        stype: ShapeType::Circle,
        C: Point::new(0.3, 0.3),
        r: 0.1,
        th: th,
        ..Default::default()
    };

    let mut circ2 = Shape {
        stype: ShapeType::Circle,
        C: Point::new(-0.3, -0.1),
        r: 0.2,
        th: th,
        ..Default::default()
    };

    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);

    // ******************** FORWARD PASS

    // Iterate over the coordinates and pixels of the image
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let xf = x as f32 / imgx as f32 - 0.5;
        let yf = y as f32 / imgy as f32 - 0.5;
        //let r = (xf / imgx as f32 * 128.0) as u8;
        //let b = (yf / imgy as f32 * 128.0) as u8;

        let p: Point = Point::new(xf, yf);
        let w1 = circ.w(p);
        let w2 = circ2.w(p);
        let w;
        if circ.sdf < circ2.sdf {
            w = w1;
        } else {
            w = w2;
        }

        let mut out_color = yellow * w;
        out_color = out_color * 255.;
        *pixel = image::Rgb([
            out_color.r as u8,
            out_color.g as u8,
            out_color.b as u8]);
    }

    println!("Rendering took {:?}", start_time.elapsed());
    start_time = Instant::now();

    // Save the image as “fractal.png”, the format is deduced from the path
    imgbuf.save("reference.jpg").unwrap();

    // ******************** BACKWARD PASS

    let refimg = loadsdf::loadimage("reference.jpg");

    // let mut circ_c: Point = Point::var_f(0.2, 0.);
    // let mut circ_r = Dp::var_f(0.2);
    let mut circ = Shape {
        stype: ShapeType::Circle,
        C: Point::new(-0.2, -0.2),
        r: 0.2,
        th: th,
        ..Default::default()
    };

    let mut circ2 = Shape {
        stype: ShapeType::Circle,
        C: Point::new(0.15, 0.2),
        r: 0.3,
        th: th,
        ..Default::default()
    };


    for _it in 0..20 {
        // let mut msed = Dp::const_f(0.);
        let mut mse: f32 = 0.;
        for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
            let xf = x as f32 / imgx as f32 - 0.5;
            let yf = y as f32 / imgy as f32 - 0.5;
            //let r = (xf / imgx as f32 * 128.0) as u8;
            //let b = (yf / imgy as f32 * 128.0) as u8;

            // again forward pass
            let p: Point = Point::new(xf, yf);

            // let sdf = circ.distance(p);
            // let w = smoothstep(0., th, -sdf);
            let w1 = circ.w(p);
            let w2 = circ2.w(p);
            let w;
            if circ.sdf < circ2.sdf {
                w = w1;
            } else {
                w = w2;
            }


            // mse, derivative
            let pix: [u8; 4] = refimg.get_pixel(x, y).0;

            let rgbref: [f32; 3] = [
                pix[0] as f32 / 255.,
                pix[1] as f32 / 255.,
                pix[2] as f32 / 255.,
            ];

            let pixmse = (
                  (w * yellow.r - rgbref[0]).powf(2.)
                + (w * yellow.g - rgbref[1]).powf(2.)
                + (w * yellow.b - rgbref[2]).powf(2.)
            );
            
            let dldw = 2. * (
                  (w * yellow.r - rgbref[0]) * yellow.r
                + (w * yellow.g - rgbref[1]) * yellow.g
                + (w * yellow.b - rgbref[2]) * yellow.b
            );
            // if x%10==0 && y%10==0 {
            //     println!("{}", dwdsdf); 
            // }

            if circ.sdf < circ2.sdf {
                circ.backward(p, dldw);
            } else {
                circ2.backward(p, dldw);
            }
            
            mse += pixmse;
            // dr += -pixmse_d * smoothstep_d(0., th, -r) * (-1.) * (-1.);

            let mut out_color = yellow * w;
            out_color = out_color * 255.;
            *pixel = image::Rgb([
                out_color.r as u8,
                out_color.g as u8,
                out_color.b as u8]);
        }
        mse = mse / imgx as f32  / imgy as f32;
        println!("mse={:.3} r={:.3} dr={:.3} C={:?}", mse, circ.r, circ.dr, circ2.C);
        circ.step(0.000001);
        circ2.step(0.000001);
        
        
        let filename: String = format!("{}{:0>2}.jpg", save_path, _it);
        imgbuf.save(filename).unwrap();
    }

}






fn main() {
    let guard = pprof::ProfilerGuardBuilder::default().frequency(1000).blocklist(&["libc", "libgcc", "pthread", "vdso"]).build().unwrap();
    small(SceneType::All, "fractal");
    // test_autodif();
    if let Ok(report) = guard.report().build() {
        let file = File::create("flamegraph.svg").unwrap();
        let mut options = pprof::flamegraph::Options::default();
        options.image_width = Some(2500);
        report.flamegraph_with_options(file, &mut options).unwrap();
    };
    
    return;
}