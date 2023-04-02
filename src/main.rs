#![allow(dead_code)]
#![allow(warnings)]
mod loadsdf;
mod libdif;
mod point;
use std::{ops::{self, Deref}, f32::consts::PI, collections::{VecDeque, HashMap}, time::Instant, rc::Rc, cell::RefCell, borrow::Borrow, fs::File};
use nalgebra as na;

use grid::Grid;
use image::{GrayImage, GenericImageView};
use libdif::{DGraph, Graph, Dp, GRAPH};
use point::Point;

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
    dsdfw: f32,
    dw: f32,
    dwm: f32,
    dwv: f32,

    h: f32,
    dsdfh: f32,
    dh: f32,
    dhm: f32,
    dhv: f32,

    // tri
    A: Point,
    dsdfA: Point,
    dA: Point,
    dAm: Point,
    dAv: Point,

    B: Point,
    dsdfB: Point,
    dB: Point,
    dBm: Point,
    dBv: Point,

    // C: Point,
    dsdfC: Point,
    // dC: Point,
    // dCm: Point,
    // dCv: Point,

    color: vec3,
    dcol: vec3,
    dcolm: vec3,
    dcolv: vec3,
}

type vec3 = na::Vector3::<f32>;

fn vec3_sqrt(v: vec3) -> vec3 {
    vec3::new(v.x.sqrt(), v.y.sqrt(), v.z.sqrt())
}

impl Shape {
    fn build(&self) -> Self {
        let mut res = Shape { ..*self };
        use ShapeType::*;
        match res.stype {
            Triangle => {
                // order points in clockwise manner
                let a = res.B - res.A;
                let b = res.C - res.A;
                //println!("> Cross product: {}", a.cross(b));
                if a.cross(b) > 0. {
                    let t = res.B;
                    res.B = res.C;
                    res.C = t;
                }

            },
            _ => ()
        }
        return res;
    }

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

                let p = point - self.C;
                let y1 = p.y.abs() - h2;
                let x1 = p.x.abs() - w2;
                let mut dw: f32 = 0.;
                let mut dh: f32 = 0.;
                if y1 < 0. && x1 < 0. {
                    // closest border
                    if y1 > x1 { 
                        sdf = y1;
                        dw=0.; dh=-0.5;
                        let dy: f32 = if p.y > 0. {-1.} else {1.};
                        self.dsdfC = Point::new(0., dy);
                    } else { 
                        sdf = x1;
                        dw=-0.5; dh=0.;
                        let dx: f32 = if p.x > 0. {-1.} else {1.};
                        self.dsdfC = Point::new(dx, 0.);
                    } 
                } else if -w2 < p.x && p.x < w2 {
                    sdf = y1;
                    // here and forward we ditch derivative
                } else if -h2 < p.y && p.y < h2 {
                    sdf = x1;
                } else {
                    sdf = Point { x: x1, y: y1 }.len();
                }
                self.dsdfw = dw;
                self.dsdfh = dh;
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
                    sdf = ca / a.len();
                    // HERE we ditch derivative, because this branch is mul by 0 anyway
                } else if cb >= 0. && 0. <= pjb && pjb <= b.len() {
                    sdf = cb / b.len();
                } else if cc >= 0. && 0. <= pjc && pjc <= c.len() {
                    sdf = cc / c.len();
                // case B
                } else if cb > 0. && pjb < 0. || ca > 0. && !(pja < 0.)  {
                    // a-b side boundary
                    sdf = pb.len();
                    // HERE we ditch derivative, because this branch is mul by 0 anyway
                } else if cc > 0. && pjc < 0. || cb > 0. && !(pjb < 0.)  {
                    // b-c side boundary
                    sdf = pc.len();
                } else if ca > 0. && pja < 0. || cc > 0. && !(pjc < 0.)  {
                    // c-a side boundary
                    sdf = pa.len();
                } else {
                    // inside triangle. cross product is negative thus max is needed

                    sdf = f32::max(ca / a.len(), 
                        f32::max(cb / b.len(), cc / c.len()));
                    let sa = ca / a.len();
                    let sb = cb / b.len();
                    let sc = cc / c.len();
                    if sa > sb && sa > sc {
                        let len = a.len();
                        if len > 0.0001 {
                            self.dsdfA = Point::new(-pa.y+a.y, pa.x-a.x) / len + a * (ca / len.powf(3.));
                            self.dsdfB = Point::new(pa.y, -pa.x) / len - a * (ca / len.powf(3.));
                        }
                    } else if sb > sc {
                        let len = b.len();
                        if len > 0.0001 {
                            self.dsdfB = Point::new(-pb.y+b.y, pb.x-b.x) / len + b * (cb / len.powf(3.));
                            self.dsdfC = Point::new(pb.y, -pb.x) / len - b * (cb / len.powf(3.));
                        }
                    } else {
                        let len = c.len();
                        if len > 0.0001 {
                            self.dsdfC = Point::new(-pc.y+c.y, pc.x-c.x) / len + c * (cc / len.powf(3.));
                            self.dsdfA = Point::new(pc.y, -pc.x) / len - c * (cc / len.powf(3.));
                        }
                    }
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

    fn backward(&mut self, point: Point, dldw: f32, drgb: vec3) {
        let dldsdf = dldw * smoothstep_d(0., self.th, -self.sdf) * (-1.);
        let p = point;
        unsafe {
            if PRINT_NOW {
                println!("SDF={}", self.sdf);
            }
        }

        self.dcol += drgb;

        use ShapeType::*;
        match self.stype {
            Circle => {
                self.dr = self.dr - dldsdf;
                if self.p_len > 0.00001 { // avoid div by 0
                    self.dC = self.dC + (-point + self.C) * (1. / self.p_len) * dldsdf * 50.;
                }
            },
            Rectangle => {
                self.dw += self.dsdfw * dldsdf;
                self.dh += self.dsdfh * dldsdf;
                self.dC = self.dC + self.dsdfC * dldsdf;
                self.dsdfw = 0.;
                self.dsdfh = 0.;
                self.dsdfC = Point::new(0., 0.);
            },
            Triangle => {
                self.dA = self.dA + self.dsdfA * dldsdf;
                self.dB = self.dB + self.dsdfB * dldsdf;
                self.dC = self.dC + self.dsdfC * dldsdf;
                self.dsdfA = Point::new(0., 0.);
                self.dsdfB = Point::new(0., 0.);
                self.dsdfC = Point::new(0., 0.);
            },
            _ => panic!("Unimplemented backward")
        }
        
    }

    fn step(&mut self, lr: f32) {
        // ** Adam params
        let beta1 = 0.9;
        let beta2 = 0.999;
        let k = 1.;
        let eps = 1e-8;
        let vec3_eps = vec3::new(eps, eps, eps);

        macro_rules! adam{
            // match like arm for macro
               ($m:expr,$v:expr,$g:expr,$x:expr)=>{
                // input: prev momentum, prev v, cur grad, prev x, lr
                // output: new momentum, new v, new wx
                // macro expand to this code
                {
                    // $a and $b will be templated using the value/variable provided to macro
                    $m = $m * beta1 + $g * (1. - beta1);
                    $v = $v * beta2 + $g * $g * (1. - beta2);
                    let m_unbias = $m / (1. - beta1.powf(k));
                    let v_unbias = $v / (1. - beta2.powf(k));
                    let grad = m_unbias / (v_unbias.powf(0.5) + eps);
                    $x = $x - grad * lr;
                }
            }
        }

        macro_rules! adam_vec{
            // match like arm for macro
               ($m:expr,$v:expr,$g:expr,$x:expr)=>{
                // input: prev momentum, prev v, cur grad, prev x, lr
                // output: new momentum, new v, new wx
                // macro expand to this code
                {
                    // $a and $b will be templated using the value/variable provided to macro
                    $m = $m * beta1 + $g * (1. - beta1);
                    $v = $v * beta2 + $g.component_mul(&$g) * (1. - beta2);
                    let m_unbias = $m * (1. / (1. - beta1.powf(k)));
                    let v_unbias = $v * (1. / (1. - beta2.powf(k)));
                    let grad = m_unbias.component_div(&(vec3_sqrt(v_unbias) + vec3_eps));
                    $x = $x - grad * lr;
                }
            }
        }

        adam_vec!(self.dcolm, self.dcolv, self.dcol, self.color);
        // self.color = self.color - self.dcol * lr;
        self.dcol = vec3::new(0., 0., 0.);

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
                adam!(self.dwm, self.dwv, self.dw, self.w);
                self.dw = 0.;

                // self.h -= self.dh * lr;
                adam!(self.dhm, self.dhv, self.dh, self.h);
                self.dh = 0.;

                adam!(self.dCm, self.dCv, self.dC, self.C);
                self.dC = Point::new(0., 0.);
            },
            Triangle => {
                adam!(self.dAm, self.dAv, self.dA, self.A);
                adam!(self.dBm, self.dBv, self.dB, self.B);
                adam!(self.dCm, self.dCv, self.dC, self.C);
                // self.A = self.A - self.dA * lr;
                // self.B = self.B - self.dB * lr;
                // self.C = self.C - self.dC * lr;
                self.dA = Point::new(0., 0.);
                self.dB = Point::new(0., 0.);
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
//     colors: Vec<vec3>, // for each face
//     tc: Vec<Point>,
//     texture: GrayImage,
// }

// impl TriangleMesh {
//     fn render(&self, point: Point, color: &mut vec3) -> bool {
//         *color = vec3 { r: 1., g: 1., b: 1., layer: 0 };
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

    let yellow = vec3::new(255./255., 220./255., 3./255.);
    let red = vec3::new(255./255., 0., 0.);
    let black = vec3::new(0., 0., 0.);
    
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
        C: Point::new(0., 0.),
        w: 0.5,
        h: 0.6,
        th: th,
        ..Default::default()
    };

    let mut tri = Shape {
        stype: ShapeType::Triangle,
        A: Point::new(0.2, 0.2),
        B: Point::new(-0.2, 0.2),
        C: Point::new(-0.2, -0.2),
        th: th,
        ..Default::default()
    }.build();



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
        let alpha2 = tri.distance_alpha(p);
        let w1 = 0. * alpha1;
        let w2 = 1. * alpha2;
        

        let mut out_color = yellow * w1 + yellow * w2;
        out_color = out_color * 255.;
        *pixel = image::Rgb([
            out_color.x as u8,
            out_color.y as u8,
            out_color.z as u8]);
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
        C: Point::new(0.3, 0.15),
        w: 0.1,
        h: 0.2,
        th: th,
        color: black,
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

    let mut tri = Shape {
        stype: ShapeType::Triangle,
        A: Point::new(0.4, 0.5),
        B: Point::new(-0.4, 0.2),
        C: Point::new(-0.3, -0.5),
        th: th,
        color: red,
        ..Default::default()
    }.build();


    let mut shapes: Vec<Shape> = Vec::new();
    shapes.push(tri);
    let nshapes = shapes.len();

    // alpha = smoothstep(-sdf)
    let mut smstep: Vec<f32> = vec![0.; nshapes];
    // color = sum(weight_i * color_i)
    let mut weight1: Vec<f32> = vec![0.; nshapes];
    // compute d pixel/d smoothstep using postfix sum
    let mut dstep_cum: Vec<vec3> = vec![black; nshapes];

    for _it in 0..80 {
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
                smstep[i] = shapes[i].distance_alpha(p);
                if i == 0 {
                    weight1[i] = 1.;
                } else {
                    weight1[i] = 1. - weight1[i - 1] * smstep[i - 1];
                }
                out_color = out_color + shapes[i].color * (weight1[i] * smstep[i]);
            }


            // mse, derivative
            let pix: [u8; 4] = refimg.get_pixel(x, y).0;
            let pixref = vec3::new(pix[0] as f32 / 255., pix[1] as f32 / 255., pix[2] as f32 / 255.);
            let pixdif = out_color - pixref;

            let pixmse = pixdif.dot(&pixdif); // sum of squares

            // compute postfix sums from end to begin for d pixel/d smoothstep
            dstep_cum[nshapes-1] = shapes[nshapes-1].color;
            for i in 1..nshapes {
                let j = nshapes - i - 1;
                // had to derive this formula
                dstep_cum[j] = shapes[j].color + dstep_cum[j+1] * (-smstep[j+1]);
            }

            // pass derivative to each shape
            for i in 0..nshapes {
                // 3 components in derivative: rgb
                let dpixdw: vec3 = dstep_cum[i] * weight1[i];
                let dldw = pixdif.dot(&dpixdw) * 2.; // coors are factored out of rgb channels, sum is ok
                
                let dpixdrgb = weight1[i] * smstep[i];
                let drgb = pixdif * dpixdrgb;

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
                shapes[i].backward(p, dldw, drgb);
                unsafe {
                    PRINT_NOW = false;
                }
            }
            
            mse += pixmse;

            
            out_color = out_color * 255.;
            *pixel = image::Rgb([
                out_color.x as u8,
                out_color.y as u8,
                out_color.z as u8]);
        }
        mse = mse / imgx as f32 / imgy as f32;
        println!("mse={:.9}", mse);
        let mut lr = 0.01;
        if _it > 40 {
            lr *= 15.;
        }
        for i in 0..nshapes {
            shapes[i].step(lr);
        }
        
        
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
