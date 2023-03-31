#![allow(dead_code)]
#![allow(warnings)]
mod loadsdf;
mod libdif;
mod point;
mod hashmap;
use std::{ops::{self, Deref}, f32::consts::PI, collections::{VecDeque, HashMap}, time::Instant, rc::Rc, cell::RefCell, borrow::Borrow, fs::File};

use grid::Grid;
use image::{GrayImage, GenericImageView};
//use libdif::{DGraph, Graph, Dp, GRAPH};

use point::Point;

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

    dr: f32,
    drm: f32,
    drv: f32,
    dC: Point,
    dCm: Point,
    dCv: Point,

    // rect
    w: f32,
    h: f32,

    dw: f32,
    dwm: f32,
    dwv: f32,

    dh: f32,
    dhm: f32,
    dhv: f32,
    
    

    color: ShapeColor
}

#[derive(Default, Copy, Clone, Debug, PartialEq)]
struct ShapeColor {
    r: f32,
    g: f32,
    b: f32
}

impl ShapeColor {
    fn new(r: i32, g: i32, b: i32) -> ShapeColor {
        ShapeColor { r: r as f32 / 255., g: g as f32 / 255., b: b as f32 / 255. }
    }
}

impl ops::Add<ShapeColor> for ShapeColor {
    type Output = ShapeColor;

    fn add(self, _rhs: ShapeColor) -> ShapeColor {
        //println!("> Point.add(Point) was called");
        return ShapeColor {
            r: self.r + _rhs.r,
            g: self.g + _rhs.g,
            b: self.b + _rhs.b
        };
    }
}

impl ops::Mul<f32> for ShapeColor {
    type Output = ShapeColor;

    fn mul(self, rhs: f32) -> Self::Output {
        return ShapeColor {
            r: self.r * rhs,
            g: self.g * rhs,
            b: self.b * rhs,
        };
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
            Rectangle => {
                let w2 = self.w / 2.;
                let h2 = self.h / 2.;

                let y1 = p.y.abs() - h2;
                let x1 = p.x.abs() - w2;
                if y1 < 0. && x1 < 0. {
                    if y1 > x1 { sdf = y1; } else { sdf = x1; } // closest border
                } else if -w2 < p.x && p.x < w2 {
                    sdf = y1;
                } else if -h2 < p.y && p.y < h2 {
                    sdf = x1;
                } else {
                    let len = Point {
                        x: x1,
                        y: y1
                    }.len();
                    self.p_len = len;
                    sdf = len;
                }
            },
            _ => panic!("Unimplemented shape")
        }
        self.sdf = sdf;
        return sdf;
    }

    fn distance_alpha(&mut self, point: Point) -> f32 {
        let sdf = self.distance_pre(point);
        return smoothstep(0., self.th, -sdf);
    }

    fn backward(&mut self, point: Point, dldw: f32) {
        let dldsdf = dldw * smoothstep_d(0., self.th, -self.sdf) * (-1.);
        let p = point;
        unsafe {
            if PRINT_NOW {
                println!("SDF={}", self.sdf);
            }
        }
        

        use ShapeType::*;
        match self.stype {
            Circle => {
                self.dr = self.dr - dldsdf;
                if self.p_len > 0.00001 { // avoid div by 0
                    self.dC = self.dC + (-point + self.C) * (1. / self.p_len) * dldsdf * 50.;
                }
            },
            Rectangle => {
                unsafe {
                    if PRINT_NOW {
                        println!("Hello from rectangle! dsdf={}", dldsdf);
                    }
                }
                let w2 = self.w / 2.;
                let h2 = self.h / 2.;

                let y1 = p.y.abs() - h2;
                let x1 = p.x.abs() - w2;
                let dw: f32;
                let dh: f32;
                if y1 < 0. && x1 < 0. {
                    if y1 > x1 {
                        dw = 0.; dh = -0.5;
                    } else { 
                        dw = -0.5; dh = 0.;
                    } // closest border
                } else if -w2 < p.x && p.x < w2 {
                    dw = 0.; dh = -0.5;
                } else if -h2 < p.y && p.y < h2 {
                    dw = -0.5; dh = 0.;
                } else {
                    // avoid division by 0!
                    if self.p_len > 0.00001 {
                        dw = x1 * (-0.5) / self.p_len;
                        dh = y1 * (-0.5) / self.p_len;
                    } else {
                        dw = 0.;
                        dh = 0.;
                    }
                }
                unsafe {
                    if PRINT_NOW {
                        println!("dw={} dh={}", dw, dh);
                    }
                }
                self.dw += dw * dldsdf;
                self.dh += dh * dldsdf;
                unsafe {
                    if PRINT_NOW {
                        println!("self.dw={} dh={}", self.dw, self.dh);
                    }
                }
            },
            _ => ()
        }
        
    }

    fn step(&mut self, lr: f32) {
        // ** Adam params
        let beta1 = 0.9;
        let beta2 = 0.999;
        let k = 1.;
        let eps = 1e-8;

        macro_rules! adam_f{
            // match like arm for macro
               ($m:expr,$v:expr,$g:expr,$x:expr)=>{
                // input: previous momentum, v, x, current grad 
                // output: new momentum, v, x
            // macro expand to this code
                   {
            // $a and $b will be templated using the value/variable provided to macro
                        $m = $m * beta1 + $g * (1. - beta1);
                        $v = $v * beta2 + $g * $g * (1. - beta2);
                        let m_unbias = $m / (1. - beta1.powf(k));
                        let v_unbias = $v / (1. - beta2.powf(k));
                        let grad = m_unbias / ( v_unbias.powf(0.5) + eps );
                        $x = $x - grad * lr;
                   }
               }
           }

        use ShapeType::*;
        match self.stype {
            Circle => {
                // ** Adam
                
                self.drm = self.drm * beta1 + self.dr * (1. - beta1);
                self.drv = self.drv * beta2 + self.dr * self.dr * (1. - beta2);
                let m_unbias = self.drm / (1. - beta1.powf(k));
                let v_unbias = self.drv / (1. - beta2.powf(k));
                let grad = m_unbias / ( v_unbias.powf(0.5) + eps );
                self.r = self.r - grad * lr;


                // ** SGD momentum
                // self.drm = self.drm * 0.2 + self.dr * 0.8;
                // self.r = self.r - self.dr * lr;

                // ** SGD
                // self.r = self.r - self.dr * lr;
                self.dr = 0.;

                
                // ** Adam
                self.dCm = self.dCm * beta1 + self.dC * (1. - beta1);
                self.dCv = self.dCv * beta2 + self.dC * self.dC * (1. - beta2);
                let m_unbias = self.dCm * (1. / (1. - beta1.powf(k)));
                let v_unbias = self.dCv * (1. / (1. - beta2.powf(k)));
                let grad = m_unbias * ( v_unbias.powf(0.5) + eps ).powf(-1.);
                self.C = self.C - grad * lr;
                
                // ** SGD momentum
                // self.dCm = self.dCm * 0.2 + self.dC * 0.8;
                // self.C = self.C - self.dCm * lr;

                // ** SGD
                // self.C = self.C - self.dC * lr;
                self.dC = Point::new(0., 0.);
            },
            Rectangle => {
                // self.w -= self.dw * lr;
                adam_f!(self.dwm, self.dwv, self.dw, self.w);
                self.dw = 0.;

                // self.h -= self.dh * lr;
                adam_f!(self.dhm, self.dhv, self.dh, self.h);
                self.dh = 0.;
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

static mut PRINT_NOW: bool = false;

fn small(scene: SceneType, save_path: &str) {
    println!("\n=== Scene {:?}", scene);
    let mut start_time = Instant::now();

    let yellow = ShapeColor::new(255, 220, 3);
    let red = ShapeColor::new(255, 0, 0);
    let black = ShapeColor{r: 0., g: 0., b: 0.};
    
    println!("Initialization took {:?}", start_time.elapsed());
    start_time = Instant::now();

    let imgx = 100;
    let imgy = imgx;

    let th = 0.1;
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

    let mut rect = Shape {
        stype: ShapeType::Rectangle,
        w: 0.5,
        h: 0.6,
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
        let alpha1 = circ.distance_alpha(p);
        let alpha2 = rect.distance_alpha(p);
        let w1 = 0. * alpha1;
        let w2 = 1. * alpha2;
        

        let mut out_color = yellow * w1 + yellow * w2;
        out_color = out_color * 255.;
        *pixel = image::Rgb([
            out_color.r as u8,
            out_color.g as u8,
            out_color.b as u8]);
    }

    println!("Rendering took {:?}", start_time.elapsed());
    start_time = Instant::now();

    // Save the image as “fractal.png”, the format is deduced from the path
    std::fs::create_dir("anim");
    imgbuf.save("anim/reference.jpg").unwrap();

    // ******************** BACKWARD PASS

    let refimg = loadsdf::loadimage("anim/reference.jpg");

    // let mut circ = Shape {
    //     stype: ShapeType::Circle,
    //     C: Point::new(-0.2, -0.2),
    //     r: 0.1,
    //     th: th,
    //     color: yellow,
    //     ..Default::default()
    // };

    let mut rect1 = Shape {
        stype: ShapeType::Rectangle,
        w: 0.1,
        h: 0.5,
        th: th,
        color: yellow,
        ..Default::default()
    };

    let mut circ2 = Shape {
        stype: ShapeType::Circle,
        C: Point::new(0.1, 0.15),
        r: 0.3,
        th: th,
        color: yellow,
        ..Default::default()
    };


    let mut shapes: Vec<Shape> = Vec::new();
    shapes.push(rect1);
    let nshapes = shapes.len();

    // alpha = smoothstep(-sdf)
    let mut alpha: Vec<f32> = vec![0.; nshapes];
    // color = sum(weight_i * color_i)
    let mut weight1: Vec<f32> = vec![0.; nshapes];
    // compute d pixel/d smoothstep using postfix sum
    let mut dstep_cum: Vec<ShapeColor> = vec![black; nshapes];

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
            let mut out_color = black;

            // let sdf = circ.distance(p);
            // let w = smoothstep(0., th, -sdf);
            for i in 0..nshapes {
                alpha[i] = shapes[i].distance_alpha(p);
                if i == 0 {
                    weight1[i] = 1.;
                } else {
                    weight1[i] = 1. - weight1[i - 1] * alpha[i - 1];
                }
                out_color = out_color + shapes[i].color * (weight1[i] * alpha[i]);
            }


            // mse, derivative
            let pix: [u8; 4] = refimg.get_pixel(x, y).0;
            let pixref = ShapeColor::new(pix[0] as i32, pix[1] as i32, pix[2] as i32);

            let pixmse = (
                  (out_color.r - pixref.r).powf(2.)
                + (out_color.g - pixref.g).powf(2.)
                + (out_color.b - pixref.b).powf(2.)
            );

            // compute postfix sums from end to begin for d pixel/d smoothstep
            dstep_cum[nshapes-1] = shapes[nshapes-1].color;
            for i in 1..nshapes {
                let j = nshapes - i - 1;
                // had to derive this formula
                dstep_cum[j] = shapes[j].color + dstep_cum[j+1] * (-alpha[j+1]);
            }

            // pass derivative to each shape
            for i in 0..nshapes {
                // 3 components in derivative: rgb
                let dpixdw = dstep_cum[i] * weight1[i];
                let dldw = 2. * (
                    (out_color.r - pixref.r) * dpixdw.r
                  + (out_color.g - pixref.g) * dpixdw.g
                  + (out_color.b - pixref.b) * dpixdw.b
                );
                // if x == imgx / 2 && y == imgy / 2 {
                //     println!("Pixel color={:?} dpix/dw={:?}", shapes[0].color, dpixdw);
                // }
                // if x == imgx / 2 && (y == imgy / 4 || y == imgy / 4 * 3 || y == imgy/2) {
                //     unsafe {
                //         PRINT_NOW = true;
                //     }
                //     println!("> x={} y={} dldw={}", x, y, dldw);
                //     println!("dpixw.r={} out.r={} ref.r={} dldw={}", dpixdw.r, out_color.r, - pixref.r, dldw);
                // }
                shapes[i].backward(p, dldw);
                unsafe {
                    PRINT_NOW = false;
                }
            }
            
            mse += pixmse;

            
            out_color = out_color * 255.;
            *pixel = image::Rgb([
                out_color.r as u8,
                out_color.g as u8,
                out_color.b as u8]);
        }
        mse = mse / imgx as f32 / imgy as f32;
        println!("mse={:.9}", mse);
        let lr = 0.01;
        for i in 0..nshapes {
            shapes[i].step(lr);
        }
        // circ.step(lr);
        // circ2.step(lr);
        
        
        let filename: String = format!("anim/{}{:0>3}.jpg", save_path, _it);
        imgbuf.save(filename).unwrap();
    }

}






fn main() {
    // let guard = pprof::ProfilerGuardBuilder::default().frequency(1000).blocklist(&["libc", "libgcc", "pthread", "vdso"]).build().unwrap();
    small(SceneType::All, "fractal");

    // if let Ok(report) = guard.report().build() {
    //     let file = File::create("flamegraph.svg").unwrap();
    //     let mut options = pprof::flamegraph::Options::default();
    //     options.image_width = Some(2500);
    //     report.flamegraph_with_options(file, &mut options).unwrap();
    // };
    
    return;
}