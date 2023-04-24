#![allow(dead_code)]
#![allow(warnings)]
mod loadsdf;
mod libdif;
mod point;
use std::{ops::{self, Deref}, f32::consts::PI, collections::{VecDeque, HashMap}, time::Instant, rc::Rc, cell::RefCell, borrow::Borrow, fs::File};
use rand;
use nalgebra as na;
use ndarray::{Array3, Array};

use image::{GrayImage, GenericImageView, RgbImage};
use libdif::{DGraph, Graph, Dp, GRAPH};
use point::Point;

#[derive(Default, Copy, Clone, Debug, PartialEq)]
enum ShapeType {
    #[default]
    Circle,
    
    Rectangle,
    Triangle,
    LineSegment,
    Bitmap,
    Star5,
    Horseshoe,
    Moon,
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
    dsdfr: f32,
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

    // image sdf
    // C - corner
    // w, h - size
    test_bitmap: bool,
    bitmap: Array3<f32>,
    dbitmap: Array3<f32>,
    dbitmapm: Array3<f32>,
    dbitmapv: Array3<f32>,
    field_scale: f32, // d *= field_scale

    // Star5
    rs: f32,
    dsdfrs: f32,
    drs: f32,
    drsm: f32,
    drsv: f32,
    
    rb: f32,
    dsdfrb: f32,
    drb: f32,
    drbm: f32,
    drbv: f32,

    no_smooth: bool,
}

type vec3 = na::Vector3::<f32>;

fn vec3_sqrt(v: vec3) -> vec3 {
    vec3::new(v.x.sqrt(), v.y.sqrt(), v.z.sqrt())
}

impl Shape {
    fn build_tri(mut self) -> Self {
        use ShapeType::*;
        match self.stype {
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
            _ => ()
        }
        return self;
    }

    fn build_bitmap(mut self, bitmap_h: usize, bitmap_w: usize) -> Self {
        // without this initialization derivative will be 0
        self.bitmap = Array3::<f32>::zeros((bitmap_h as usize, bitmap_w as usize, 1)) + 0.5 + self.th / 2.;

        self.dbitmap = Array3::<f32>::zeros((bitmap_h as usize, bitmap_w as usize, 1));
        self.dbitmapm = Array3::<f32>::zeros((bitmap_h as usize, bitmap_w as usize, 1));
        self.dbitmapv = Array3::<f32>::zeros((bitmap_h as usize, bitmap_w as usize, 1));
        return self;
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
            Bitmap => {
                let (bitmap_h, bitmap_w, _c) = self.bitmap.dim();
                let q = (p - self.C) * Point::new(bitmap_w as f32 / self.w, bitmap_h as f32 / self.h);

                let scale: f32 = if self.field_scale != 0. { self.field_scale } else { 1. };
                let dist_large: f32 = 1000.;

                let x = q.x.floor() as i32;
                let y = q.y.floor() as i32;

                let kx = q.x - x as f32;
                let ky = q.y - y as f32;

                let k1 = (1.-kx)*(1.-ky);
                let k2 = (1.-kx)*(ky);
                let k3 = (kx)*(ky);
                let k4 = (kx)*(1.-ky);

                let mut get_sdf = |px: i32, py: i32, k: f32| -> f32 {
                    if 0 <= px && px < bitmap_w as i32 && 0 <= py && py < bitmap_h as i32 {
                        0.5 - self.bitmap[[py as usize, px as usize, 0]]
                    } else {
                        dist_large
                    }
                };
                
                let v1 = get_sdf(x, y, k1);
                let v2 = get_sdf(x, y+1, k2);
                let v3 = get_sdf(x+1, y+1, k3);
                let v4 = get_sdf(x+1, y, k4);
                
                sdf = (k1*v1 + k2*v2 + k3*v3 + k4*v4)*scale; // TODO: this is not true distance, please think about it
            },
            Star5 => {
                let p = p - self.C;
                let p = Point::new(p.x, -p.y);
                let alpha = PI / 5.;
                let p = p.abs_x().rotate(alpha*2.)
                    .abs_x().rotate(alpha*2.)
                    .abs_x().rotate(alpha)
                    .abs_x();
                let vec_small = Point::new(0., self.rs);
                let vec_big = Point::new(0., self.rb).rotate(-alpha);
                let sb = vec_big - vec_small;

                let sp = p - vec_small;
                let perp = sp.cross(sb) / sb.len();

                // drs
                let cross = sp.cross(sb);
                let len = sb.len();
                let dcross = sb.x - sp.x;
                let dlen = -sb.y * (1. / len);
                self.dsdfrs = (cross * dlen - dcross * len) / (len*len);

                // drb
                let dcross = -sp.y * (-alpha.sin()) + sp.x * alpha.cos();
                let dlen = (sb.x*(-alpha.sin()) + sb.y*alpha.cos()) * (1. / len);
                self.dsdfrb = (cross * dlen - dcross * len) / (len*len);

                sdf = -perp;
            },
            Horseshoe => {
                let p = (point - self.C).abs_x();
                let a = (self.rs + self.rb) / 2.;
                let a2 = (self.rb - self.rs) / 2.;
                let alpha = PI / 10.;
                let w = Point::new(alpha.cos(), alpha.sin());
                
                let p1 = w * a;
                let b = 0.1;
                
                if w.cross(p) < 0. {
                    let pc = point - self.C;
                    let sdf1 = p.len() - a;
                    sdf = sdf1.abs() - (a - self.rs);
                    if p.len() > 0.0001 {
                        self.dsdfC = p / p.len() * sdf1.signum() * Point::new(pc.x.signum(), 1.) * (-1.);
                    }
                    self.dsdfrs = sdf1.signum() * (-0.5) + 0.5;
                    self.dsdfrb = sdf1.signum() * (-0.5) - 0.5;
                } else {
                    let p2 = p1 + w.normal() * b;
                    let p1p2 = p2 - p1;
                    let p1p = p - p1;
                    let perp = p1p.cross(p1p2).abs() / p1p2.len() - a2;
                    let p2p = p - p2;
                    let perp2 = p1.cross(p2p).abs() / p1.len() - a2;
                    if w.cross(p2p) < 0. {
                        sdf = perp;
                    } else {
                        sdf = perp.max(perp2);
                    }
                }
            },
            Moon => {
                let s1 = (point - self.C).len() - self.r;
                let s2 = (point - self.B).len() - self.rb;
                if -s1 < s2 {
                    sdf = s1;
                    let p = point - self.C;
                    if p.len() > 0.0001 {
                        self.dsdfC = -p / p.len();
                    }
                    self.dsdfr = -1.;
                } else {
                    sdf = -s2;
                    let p = point - self.B;
                    if p.len() > 0.0001 {
                        self.dsdfB = p / p.len();
                    }
                    self.dsdfrb = 1.;
                };
            },
            _ => panic!("Unimplemented shape")
        }
        self.sdf = sdf;
        return sdf;
    }

    fn distance_alpha(&mut self, point: Point) -> f32 {
        let sdf = self.distance_pre(point);
        if self.no_smooth { return if sdf < 0. {1.} else {0.}; }
        return smoothstep(0., self.th, -sdf);
    }

    fn backward(&mut self, point: Point, dldw: f32, drgb: vec3) {
        let dldsdf = dldw * smoothstep_d(0., self.th, -self.sdf) * (-1.);
        let p = point;

        self.dcol += drgb;

        use ShapeType::*;
        // Extra: SDF edge sampling
        if self.no_smooth {
            if self.sdf > 0. || self.sdf < -1. / 128. {return};
            match self.stype {
                Circle => {
                    self.dr += dldw;
                    let normal = p - self.C; // insignificat error
                    self.dC += normal * dldw;
                },
                Rectangle => {
                    let w2 = self.w / 2.;
                    let h2 = self.h / 2.;

                    let p = point - self.C;
                    let y1 = p.y.abs() - h2;
                    let x1 = p.x.abs() - w2;
                    if y1 > x1 {
                        self.dh += dldw; 
                        let normal = if p.y > 0. {1.} else {-1.};
                        self.dC += Point::new(0., normal) * dldw;
                    } else {
                        self.dw += dldw;
                        let normal = if p.x > 0. {1.} else {-1.};
                        self.dC += Point::new(normal, 0.) * dldw;
                    }
                },
                _ => panic!("Unimplemented backward")
            }
            return;
        }

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
            Bitmap => {
                let (bitmap_h, bitmap_w, _c) = self.bitmap.dim();
                let q = (p - self.C) * Point::new(bitmap_w as f32 / self.w, bitmap_h as f32 / self.h);

                let scale: f32 = if self.field_scale != 0. { self.field_scale } else { 1. };
                let dist_large: f32 = 1000.;

                let x = q.x.floor() as i32;
                let y = q.y.floor() as i32;

                let kx = q.x - x as f32;
                let ky = q.y - y as f32;

                let k1 = (1.-kx)*(1.-ky);
                let k2 = (1.-kx)*(ky);
                let k3 = (kx)*(ky);
                let k4 = (kx)*(1.-ky);

                let mut get_sdf = |px: i32, py: i32, k: f32| -> f32 {
                    if 0 <= px && px < bitmap_w as i32 && 0 <= py && py < bitmap_h as i32 {
                        self.dbitmap[[py as usize, px as usize, 0]] += k * dldsdf;
                        0.
                    } else {
                        dist_large
                    }
                };
                
                let v1 = get_sdf(x, y, k1);
                let v2 = get_sdf(x, y+1, k2);
                let v3 = get_sdf(x+1, y+1, k3);
                let v4 = get_sdf(x+1, y, k4);
            },
            Star5 => {
                self.drs += self.dsdfrs * dldsdf;
                self.drb += self.dsdfrb * dldsdf;
                self.dsdfrs = 0.;
                self.dsdfrb = 0.;

                let pc = self.C - point;
                if pc.len() > 0.00001 {
                    self.dC = self.dC + pc / pc.len() * dldsdf;
                }
            },
            Horseshoe => {
                self.dC = self.dC + self.dsdfC * dldsdf;
                self.dsdfC = Point::new(0., 0.);

                self.drs += self.dsdfrs * dldsdf;
                self.drb += self.dsdfrb * dldsdf;
                self.dsdfrs = 0.;
                self.dsdfrb = 0.;
            },
            Moon => {
                self.dC = self.dC + self.dsdfC * dldsdf;
                self.dsdfC = Point::new(0., 0.);
                self.dB = self.dB + self.dsdfB * dldsdf;
                self.dsdfB = Point::new(0., 0.);

                self.dr += self.dsdfr * dldsdf;
                self.drb += self.dsdfrb * dldsdf;
                self.dsdfr = 0.;
                self.dsdfrb = 0.;
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
                // output: new momentum, new v, new x
                // macro expand to this code
                {
                    // $a and $b will be templated using the value/variable provided to macro
                    $m = $m * beta1 + $g * (1. - beta1);
                    $v = $v * beta2 + $g * $g * (1. - beta2);
                    let m_unbias = $m * (1. / (1. - beta1.powf(k)));
                    let v_unbias = $v * (1. / (1. - beta2.powf(k)));
                    let grad = m_unbias / (v_unbias.powf(0.5) + eps);
                    $x = $x - grad * lr;
                }
            }
        }

        macro_rules! adam_vec{
            // match like arm for macro
               ($m:expr,$v:expr,$g:expr,$x:expr,$lr:expr)=>{
                // input: prev momentum, prev v, cur grad, prev x, lr
                // output: new momentum, new v, new x
                // macro expand to this code
                {
                    // $a and $b will be templated using the value/variable provided to macro
                    $m = $m * beta1 + $g * (1. - beta1);
                    $v = $v * beta2 + $g.component_mul(&$g) * (1. - beta2);
                    let m_unbias = $m * (1. / (1. - beta1.powf(k)));
                    let v_unbias = $v * (1. / (1. - beta2.powf(k)));
                    let grad = m_unbias.component_div( &(vec3_sqrt(v_unbias) + vec3_eps) );
                    $x = $x - grad * $lr;
                }
            }
        }

        macro_rules! adam_ndarray{
            // match like arm for macro
               ($m:expr,$v:expr,$g:expr,$x:expr)=>{
                // input: prev momentum, prev v, cur grad, prev x, lr
                // output: new momentum, new v, new x
                // macro expand to this code
                {
                    // $a and $b will be templated using the value/variable provided to macro
                    $m = &$m * beta1 + &$g * (1. - beta1);
                    $v = &$v * beta2 + &$g * &$g * (1. - beta2);
                    let m_unbias = &$m * (1. / (1. - beta1.powf(k)));
                    let mut v_unbias = &$v * (1. / (1. - beta2.powf(k)));
                    let grad = &m_unbias / &(v_unbias.mapv(f32::sqrt) + eps);
                    $x = &$x - &grad * lr;
                }
            }
        }

        

        adam_vec!(self.dcolm, self.dcolv, self.dcol, self.color, lr*20.);
        // self.color = self.color - self.dcol * (lr*0.0001);
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
            Bitmap => {
                // self.bitmap += &(&self.dbitmap * lr); // + here because sdf_pix = 0.5 - bitmap[pix]
                self.dbitmap *= -1.;
                adam_ndarray!(self.dbitmapm, self.dbitmapv, self.dbitmap, self.bitmap);
                self.dbitmap.fill(0.);
            },
            Star5 => {
                adam!(self.drsm, self.drsv, self.drs, self.rs);
                adam!(self.drbm, self.drbv, self.drb, self.rb);
                adam!(self.dCm, self.dCv, self.dC, self.C);
                // self.rs -= self.drs * lr*0.0001;
                // self.rb -= self.drb * lr*0.001;
                // self.C = self.C - self.dC * lr*0.001;
                self.drs = 0.;
                self.drb = 0.;
                self.dC = Point::new(0., 0.);
            },
            Horseshoe => {
                adam!(self.dCm, self.dCv, self.dC, self.C);
                adam!(self.drsm, self.drsv, self.drs, self.rs);
                adam!(self.drbm, self.drbv, self.drb, self.rb);
                self.dC = Point::new(0., 0.);
                self.drs = 0.;
                self.drb = 0.;
            },
            Moon => {
                adam!(self.dCm, self.dCv, self.dC, self.C);
                adam!(self.dBm, self.dBv, self.dB, self.B);
                self.dC = Point::new(0., 0.);
                self.dB = Point::new(0., 0.);

                adam!(self.drbm, self.drbv, self.drb, self.rb);
                adam!(self.drm, self.drv, self.dr, self.r);
                self.dr = 0.;
                self.drb = 0.;
            },
            _ => ()
        }
    }
}

#[derive(Default)]
struct TriangleMesh {
    verticies: Vec<Point>,
    dverticies: Vec<Point>,
    dverticiesm: Vec<Point>,
    dverticiesv: Vec<Point>,

    indices: Vec<u32>, // triangles, must be counter-clockwise (in standard non-flipped XY plane!!!)
    colors: Vec<vec3>, // for each vertex
    is_textured: Vec<bool>, // each face is colored or textured
    tc: Vec<Point>, // texture coors for each vertex
    input_texture: RgbImage, // type int [0, 255]

    dcolors: Vec<vec3>,
    dcolorsm: Vec<vec3>,
    dcolorsv: Vec<vec3>,

    texture: Array3<f32>, // type f32 [0, 1]
    dtexture: Array3<f32>,
    dtexturem: Array3<f32>,
    dtexturev: Array3<f32>,
}

impl TriangleMesh {
    fn build(mut self, zero_init: bool, bitmap_h: usize, bitmap_w: usize) -> Self {
        self.dverticies.resize(self.verticies.len(), Point::default());
        self.dverticiesm.resize(self.verticies.len(), Point::default());
        self.dverticiesv.resize(self.verticies.len(), Point::default());

        self.dcolors.resize(self.verticies.len(), vec3::default());
        self.dcolorsm.resize(self.verticies.len(), vec3::default());
        self.dcolorsv.resize(self.verticies.len(), vec3::default());

        // init texture
        let h: usize;
        let w: usize;
        if !zero_init {
            let dim = self.input_texture.dimensions();
            w = dim.0 as usize;
            h = dim.1 as usize;
        } else {
            h = bitmap_h;
            w = bitmap_w;
        }
        // self.texture = Array3::<f32>::zeros((h, w, 3));
        self.texture = Array3::<f32>::from_elem((h, w, 3), 0.1);
        self.dtexture = Array3::<f32>::zeros((h, w, 3));
        self.dtexturem = Array3::<f32>::zeros((h, w, 3));
        self.dtexturev = Array3::<f32>::zeros((h, w, 3));
        if !zero_init {
            for i in 0..h {
                for j in 0..w {
                    for c in 0..3 {
                        self.texture[[i, j, c]] = self.input_texture.get_pixel(j as u32, i as u32).0[c] as f32 / 255.;
                    }
                }
            }
        }
        return self;
    }

    fn render(&mut self, point: Point, color: &mut vec3, backward: bool, dldcolor: vec3) -> bool {
        let p = point;
        *color = vec3::new(1., 1., 1.);
        for i in 0..self.indices.len() / 3 {
            let j = i * 3;
            // check triangle
            let ia = self.indices[j] as usize;
            let ib = self.indices[j+1] as usize;
            let ic = self.indices[j+2] as usize;
            let A = self.verticies[ia];
            let B = self.verticies[ib];
            let C = self.verticies[ic];
            

            let a = B - A;
            let b = C - B;
            let c = A - C;
            let pa = p - A;
            let pb = p - B;
            let pc = p - C;
            // on which side of triangle is our point
            let ca = a.cross(pa);
            let cb = b.cross(pb);
            let cc = c.cross(pc);

            // inside
            if ca >= 0. && cb >= 0. && cc >= 0. {
                let perpa = cb / b.len();
                let perpb = cc / c.len();
                let perpc = ca / a.len();
                let s = perpa + perpb + perpc;
                let xa = perpa / s;
                let xb = perpb / s;
                let xc = perpc / s;

                // interpolate texture attribute
                let mut sample_texture = |p: Point| -> vec3 {
                    let (bitmap_h, bitmap_w, _c): (usize, usize, usize) = self.texture.dim();
                    let q = p * Point::new(bitmap_w as f32, bitmap_h as f32);

                    let x = q.x.floor() as i32;
                    let y = q.y.floor() as i32;

                    let kx = q.x - x as f32;
                    let ky = q.y - y as f32;

                    let k1 = (1.-kx)*(1.-ky);
                    let k2 = (1.-kx)*(ky);
                    let k3 = (kx)*(ky);
                    let k4 = (kx)*(1.-ky);

                    let mut get_pixel = |px: i32, py: i32, k: f32| -> vec3 {
                        if 0 <= px && px < bitmap_w as i32 && 0 <= py && py < bitmap_h as i32 {
                            if backward {
                                self.dtexture[[py as usize, px as usize, 0]] += k * dldcolor.x; // TODO probably there's a better method in ndarray
                                self.dtexture[[py as usize, px as usize, 1]] += k * dldcolor.y;
                                self.dtexture[[py as usize, px as usize, 2]] += k * dldcolor.z;
                            }
                            vec3::new(
                            self.texture[[py as usize, px as usize, 0]],
                            self.texture[[py as usize, px as usize, 1]],
                            self.texture[[py as usize, px as usize, 2]]
                            )
                            // 0.5 - self.bitmap[[py as usize, px as usize, 0]]
                        } else {
                            vec3::default()
                        }
                    };
                    
                    let v1 = get_pixel(x, y, k1);
                    let v2 = get_pixel(x, y+1, k2);
                    let v3 = get_pixel(x+1, y+1, k3);
                    let v4 = get_pixel(x+1, y, k4);
                    
                    k1*v1 + k2*v2 + k3*v3 + k4*v4
                };

                // *** Derivative of points
                if backward {
                    let na = a.normal() / a.len();
                    let nb = b.normal() / b.len();
                    let nc = c.normal() / c.len();
                    
                    let dldcolor_sum = dldcolor.dot(&vec3::new(1., 1., 1.));
                    unsafe {
                        if PRINT_NOW {
                            println!("dldcolor_sum={}", dldcolor_sum);
                            PRINT_NOW = false;
                        }
                    }
                    let th = 2. / 256 as f32;
                    if perpa < th { // TODO which value to use here?
                        let t = b.proj(pb) / b.len();
                        self.dverticies[ib] += nb * (1. - t) * (-1.) * dldcolor_sum; // normal points inside
                        self.dverticies[ic] += nb * t * (-1.) * dldcolor_sum;
                    }
                    if perpb < th {
                        let t = c.proj(pc) / c.len();
                        self.dverticies[ic] += nc * (1. - t) * (-1.) * dldcolor_sum;
                        self.dverticies[ia] += nc * t * (-1.) * dldcolor_sum;
                    }
                    if perpc < th {
                        let t = a.proj(pa) / a.len();
                        self.dverticies[ia] += na * (1. - t) * (-1.) * dldcolor_sum;
                        self.dverticies[ib] += na * t * (-1.) * dldcolor_sum;
                    }
                }

                if self.is_textured[i] {
                    *color = sample_texture(
                        self.tc[ia] * xa
                        + self.tc[ib] * xb
                        + self.tc[ic] * xc
                    );
                } else {
                    self.dcolors[ia] += dldcolor * xa;
                    self.dcolors[ib] += dldcolor * xb;
                    self.dcolors[ic] += dldcolor * xc;
                    *color = (
                        self.colors[ia] * xa
                        + self.colors[ib] * xb
                        + self.colors[ic] * xc
                    );
                }
                return true;
            }
        }
        return false;
    }

    fn step(&mut self, lr: f32) {
        let beta1 = 0.9;
        let beta2 = 0.999;
        let k = 1.;
        let eps = 1e-8;
        let vec3_eps = vec3::new(eps, eps, eps);

        macro_rules! adam{
            // match like arm for macro
               ($m:expr,$v:expr,$g:expr,$x:expr,$lr:expr)=>{
                // input: prev momentum, prev v, cur grad, prev x, lr
                // output: new momentum, new v, new x
                // macro expand to this code
                {
                    // $a and $b will be templated using the value/variable provided to macro
                    $m = $m * beta1 + $g * (1. - beta1);
                    $v = $v * beta2 + $g * $g * (1. - beta2);
                    let m_unbias = $m * (1. / (1. - beta1.powf(k)));
                    let v_unbias = $v * (1. / (1. - beta2.powf(k)));
                    let grad = m_unbias / (v_unbias.powf(0.5) + eps);
                    $x = $x - grad * $lr;
                }
            }
        }

        macro_rules! adam_vec{
            // match like arm for macro
               ($m:expr,$v:expr,$g:expr,$x:expr,$lr:expr)=>{
                // input: prev momentum, prev v, cur grad, prev x, lr
                // output: new momentum, new v, new x
                // macro expand to this code
                {
                    // $a and $b will be templated using the value/variable provided to macro
                    $m = $m * beta1 + $g * (1. - beta1);
                    $v = $v * beta2 + $g.component_mul(&$g) * (1. - beta2);
                    let m_unbias = $m * (1. / (1. - beta1.powf(k)));
                    let v_unbias = $v * (1. / (1. - beta2.powf(k)));
                    let grad = m_unbias.component_div( &(vec3_sqrt(v_unbias) + vec3_eps) );
                    $x = $x - grad * $lr;
                }
            }
        }

        macro_rules! adam_ndarray{
            // match like arm for macro
               ($m:expr,$v:expr,$g:expr,$x:expr,$lr:expr)=>{
                // input: prev momentum, prev v, cur grad, prev x, lr
                // output: new momentum, new v, new x
                // macro expand to this code
                {
                    // $a and $b will be templated using the value/variable provided to macro
                    $m = &$m * beta1 + &$g * (1. - beta1);
                    $v = &$v * beta2 + &$g * &$g * (1. - beta2);
                    let m_unbias = &$m * (1. / (1. - beta1.powf(k)));
                    let mut v_unbias = &$v * (1. / (1. - beta2.powf(k)));
                    let grad = &m_unbias / &(v_unbias.mapv(f32::sqrt) + eps);
                    $x = &$x - &grad * $lr;
                }
            }
        }

        adam_ndarray!(self.dtexturem, self.dtexturev, self.dtexture, self.texture, lr);
        // self.texture = &self.texture - &self.dtexture * lr;
        self.dtexture.fill(0.); // TODO use this in other places

        for i in 0..self.verticies.len() {
            adam_vec!(self.dcolorsm[i], self.dcolorsv[i], self.dcolors[i], self.colors[i], lr*0.5);
            // self.colors[i] = self.colors[i] - self.dcolors[i] * lr * 0.001;
            self.dcolors[i] = vec3::default();
            
            adam!(self.dverticiesm[i], self.dverticiesv[i], self.dverticies[i], self.verticies[i], lr);
            // self.verticies[i] = self.verticies[i] - self.dverticies[i] * lr * (1.);
            self.dverticies[i] = Point::default();
        }
    }
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

pub fn task_sdf(save_path: &str) {
    println!("\n=== Task: 2 circles, 1 box ");
    let mut start_time = Instant::now();

    let yellow = vec3::new(255./255., 220./255., 3./255.);
    let red = vec3::new(255./255., 0., 0.);
    let black = vec3::new(0., 0., 0.);

    let imgx = 256;
    let imgy = imgx;
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);

    let th = 0.025;

    std::fs::create_dir("anim");

    // ******************** BACKWARD PASS

    let refimg = loadsdf::loadimage("resources/02_reference.png");

    let mut circ = Shape {
        stype: ShapeType::Circle,
        C: Point::new(0.2, 0.2),
        r: 0.1,
        th: th,
        color: vec3::new(0.543, 0.2232, 0.42),
        ..Default::default()
    };

    let mut circ2 = Shape {
        stype: ShapeType::Circle,
        C: Point::new(-0.27, 0.22),
        r: 0.12,
        th: th,
        color: vec3::new(0.1, 0.6, 1.),
        ..Default::default()
    };

    let mut rect = Shape {
        stype: ShapeType::Rectangle,
        C: Point::new(-0.19, -0.22),
        w: 0.25,
        h: 0.25,
        th: th,
        color: vec3::new(0.7, 0.7, 0.),
        ..Default::default()
    };

    let mut shapes: Vec<Shape> = Vec::new();
    shapes.push(circ);
    shapes.push(circ2);
    shapes.push(rect);
    let nshapes = shapes.len();

    
    let mut smstep: Vec<f32> = vec![0.; nshapes]; // smstep = smoothstep(-sdf)
    let mut weight1: Vec<f32> = vec![0.; nshapes]; // weight1[i] * smstep[i] * color[i]
    // compute d pixel/d smoothstep using postfix sum
    let mut dstep_cum: Vec<vec3> = vec![black; nshapes];

    println!("Initialization took {:?}", start_time.elapsed());
    start_time = Instant::now();

    for _it in 0..60 {
        let mut mse: f32 = 0.;
        for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
            let xf = x as f32 / imgx as f32 - 0.5;
            let yf = y as f32 / imgy as f32 - 0.5;
            //let r = (xf / imgx as f32 * 128.0) as u8;
            //let b = (yf / imgy as f32 * 128.0) as u8;

            // again forward pass
            let p: Point = Point::new(xf, yf);
            let mut out_color = black;

            for i in 0..nshapes {
                smstep[i] = shapes[i].distance_alpha(p);
                if i == 0 {
                    weight1[i] = 1.;
                } else {
                    weight1[i] = 1. - weight1[i - 1] * smstep[i - 1];
                }
                out_color = out_color + shapes[i].color * (weight1[i] * smstep[i]);
            }


            // *** Compute derivatives, MSE
            let pix = refimg.get_pixel(x, y).0;
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

                shapes[i].backward(p, dldw, drgb);
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
        for i in 0..nshapes {
            shapes[i].step(lr);
        }
        
        let filename: String = format!("anim/{}{:0>3}.jpg", save_path, _it);
        imgbuf.save(filename).unwrap();
    }

    println!("Rendering took {:?}", start_time.elapsed());
    start_time = Instant::now();
}

pub fn task_bitmap(save_path: &str) {
    println!("\n=== Scene");
    let mut start_time = Instant::now();

    let yellow = vec3::new(255./255., 220./255., 3./255.);
    let red = vec3::new(255./255., 0., 0.);
    let black = vec3::new(0., 0., 0.);

    let imgx = 1024;
    let imgy = imgx;
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);

    let th = 0.025;

    std::fs::create_dir("anim");

    // ******************** BACKWARD PASS

    let refimg = loadsdf::loadimage("resources/ACDC_logo.jpg");

    let mut bmp1 = Shape {
        stype: ShapeType::Bitmap,
        C: Point::new(-0.5, -0.5),
        w: 1.,
        h: 1.,
        field_scale: 1.,
        test_bitmap: true,
        th: th,
        color: red,
        ..Default::default()
    }.build_bitmap(64, 64);


    let mut shapes: Vec<Shape> = Vec::new();
    shapes.push(bmp1);
    let nshapes = shapes.len();

    
    let mut smstep: Vec<f32> = vec![0.; nshapes]; // smstep = smoothstep(-sdf)
    let mut weight1: Vec<f32> = vec![0.; nshapes]; // weight1[i] * smstep[i] * color[i]
    // compute d pixel/d smoothstep using postfix sum
    let mut dstep_cum: Vec<vec3> = vec![black; nshapes];

    let iterations = 100;
    for _it in 0..iterations {
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


            // *** Compute derivatives, MSE
            let pix = refimg.get_pixel(x, y).0;
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

            // *** pass derivative to each shape
            for i in 0..nshapes {
                // 3 components in derivative: rgb
                let dpixdw: vec3 = dstep_cum[i] * weight1[i];
                let dldw = pixdif.dot(&dpixdw) * 2.; // coors are factored out of rgb channels, sum is ok
                
                let dpixdrgb = weight1[i] * smstep[i];
                let drgb = pixdif * dpixdrgb;

                shapes[i].backward(p, dldw, drgb);
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
        let mut lr = 0.02;
        for i in 0..nshapes {
            shapes[i].step(lr);
        }
        
        
        if _it == iterations - 1 {
            let (dh, dw, _dc) = shapes[0].bitmap.dim();
            let mut img: RgbImage = image::ImageBuffer::new(dw as u32, dh as u32);
            for (x, y, pixel) in img.enumerate_pixels_mut() {
                let v = (shapes[0].bitmap[[y as usize, x as usize, 0]] - 0.5) * 255.;
                *pixel = image::Rgb([
                    v as u8,
                    v as u8,
                    v as u8]);
            }
            img.save("anim/sdf.png");  

            let filename: String = format!("anim/{}{:0>3}.jpg", save_path, _it);
            imgbuf.save(filename).unwrap();
        }
        
    }
}

pub fn task_edge_sampling(save_path: &str) {
    println!("\n=== Scene");
    let mut start_time = Instant::now();

    let yellow = vec3::new(255./255., 220./255., 3./255.);
    let red = vec3::new(255./255., 0., 0.);
    let red2 = vec3::new(128./255., 0., 0.);
    let black = vec3::new(0., 0., 0.);
    let blue = vec3::new(0., 0., 1.);

    let imgx = 256;
    let imgy = imgx;
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);
    
    let th = 0.025;
    std::fs::create_dir("anim");

    let mut mesh: TriangleMesh = TriangleMesh { 
        // verticies: vec![Point::new(-0.25, -0.25), Point::new(0.25, -0.25), Point::new(-0.25, 0.25)], 
        // indices: vec![0, 1, 2], 
        // colors: vec![red, red, red],
        // is_textured: vec![false],
        // tc: vec![Point::new(0., 0.), Point::new(1., 0.), Point::new(1., 1.)],

        verticies: vec![Point::new(-0.5, -0.5), Point::new(-0., -0.5), Point::new(-0.5, 0.),
            Point::new(-0., -0.5), Point::new(0.5, -0.5), Point::new(-0., -0.),
            Point::new(-0.5, -0.), Point::new(-0., 0.), Point::new(-0.5, 0.5)], 
        indices: vec![0, 1, 2,
            3, 4, 5,
            6, 7, 8], 
        colors: vec![red, red, red, red, red, red, red, red, red,],
        is_textured: vec![false, false, false],
        tc: vec![Point::new(0., 0.), Point::new(1., 0.), Point::new(0., 1.),
            Point::new(0., 0.), Point::new(1., 0.), Point::new(0., 1.),
            Point::new(0., 0.), Point::new(1., 0.), Point::new(0., 1.)],

        input_texture: loadsdf::loadimage("resources/brick.jpg").to_rgb8(),
        ..Default::default()
    }.build(false, 0, 0);


    {
        let mut shapes: Vec<Shape> = Vec::new();
        let nshapes = shapes.len();

        let mut smstep: Vec<f32> = vec![0.; nshapes]; // smstep = smoothstep(-sdf)
        let mut weight1: Vec<f32> = vec![0.; nshapes]; // weight1[i] * smstep[i] * color[i]

        // Create a new ImgBuf with width: imgx and height: imgy
        
        println!("Initialization took {:?}", start_time.elapsed());
        start_time = Instant::now();

        // ******************** FORWARD PASS
        

        // Iterate over the coordinates and pixels of the image
        for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
            let xf = x as f32 / imgx as f32 - 0.5;
            let yf = y as f32 / imgy as f32 - 0.5;
            //let r = (xf / imgx as f32 * 128.0) as u8;
            //let b = (yf / imgy as f32 * 128.0) as u8;

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

            // *** Mesh
            
            let mut mesh_color = black;
            if mesh.render(p, &mut mesh_color, true, vec3::default()) {
                out_color = mesh_color;
            }

            
            out_color = out_color * 255.;
            *pixel = image::Rgb([
                out_color.x as u8,
                out_color.y as u8,
                out_color.z as u8]);
        }

        println!("Rendering took {:?}", start_time.elapsed());
        start_time = Instant::now();

        // Save the image as “fractal.png”, the format is deduced from the path
        // imgbuf.save("anim/reference.png");

        // return;
    }

    // ******************** BACKWARD PASS

    let refimg = loadsdf::loadimage("resources/01_reference.png");

    let mut mesh: TriangleMesh = TriangleMesh { 
        // verticies: vec![Point::new(-0.2, -0.2), Point::new(0.25, -0.25), Point::new(-0.25, 0.25)], 
        // indices: vec![0, 1, 2], 
        // colors: vec![red2, red2, red2],
        // is_textured: vec![false],
        // tc: vec![Point::new(0., 0.), Point::new(1., 0.), Point::new(1., 1.)],

        verticies: vec![Point::new(-0.43, -0.45), Point::new(-0.04, -0.453), Point::new(-0.44, 0.01),
            Point::new(-0.05, -0.5), Point::new(0.47, -0.41), Point::new(-0.09, -0.05),
            Point::new(-0.5, -0.05), Point::new(-0.05, 0.01), Point::new(-0.5, 0.36)], 
        indices: vec![0, 1, 2,
            3, 4, 5,
            6, 7, 8], 
        colors: vec![
            red2, red2, red2,
            red2, red2, red2,
            red2, red2, red2,],
        is_textured: vec![true, true, false],
        tc: vec![Point::new(0., 0.), Point::new(1., 0.), Point::new(0., 1.),
            Point::new(0., 0.), Point::new(1., 0.), Point::new(0., 1.),
            Point::new(0., 0.), Point::new(1., 0.), Point::new(0., 1.)],

        input_texture: loadsdf::loadimage("resources/brick.jpg").to_rgb8(),
        ..Default::default()
    }.build(true, 1024, 1024);

    let mut circ = Shape {
        stype: ShapeType::Circle,
        C: Point::new(0.2, 0.2),
        r: 0.1,
        th: th,
        color: vec3::new(0.543, 0.2232, 0.42),
        ..Default::default()
    };

    {
        let mut shapes: Vec<Shape> = Vec::new();
        shapes.push(circ);
        let nshapes = shapes.len();

        let mut smstep: Vec<f32> = vec![0.; nshapes]; // smstep = smoothstep(-sdf)
        let mut weight1: Vec<f32> = vec![0.; nshapes]; // weight1[i] * smstep[i] * color[i]
        // compute d pixel/d smoothstep using postfix sum
        let mut dstep_cum: Vec<vec3> = vec![black; nshapes];

        // Create a new ImgBuf with width: imgx and height: imgy
        
        println!("Initialization took {:?}", start_time.elapsed());
        start_time = Instant::now();

        for _it in 0..60 {
            let mut mse = 0.;
            // unsafe {
            //     PRINT_NOW = true;
            // }

            // Iterate over the coordinates and pixels of the image
            for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
                let xf = x as f32 / imgx as f32 - 0.5;
                let yf = y as f32 / imgy as f32 - 0.5;
                //let r = (xf / imgx as f32 * 128.0) as u8;
                //let b = (yf / imgy as f32 * 128.0) as u8;

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

                // *** Mesh
                // if x==imgx/3 && y==imgy/4 {
                //     unsafe {
                //         PRINT_NOW = true;
                //     }
                // }
                // if x==imgx/2 && y==imgy/2 {
                //     unsafe {
                //         PRINT_NOW = true;
                //     }
                // }

                let mut mesh_color = black;
                let nsamples = 4;
                for sample in 0..nsamples {
                    let x_add = (sample % 2) as f32;
                    let y_add = (sample / 2) as f32;
                    let p_add = Point::new(
                        (0.25 + x_add * 0.5) / imgx as f32,
                        (0.25 + y_add * 0.5) / imgy as f32);
                    let p_sample = p + p_add;
                    if mesh.render(p_sample, &mut mesh_color, false, vec3::default()) {
                        out_color = out_color + mesh_color / nsamples as f32;
                    }
                }
                // if mesh.render(p, &mut mesh_color, false, vec3::default()) {
                //     out_color = mesh_color;
                // }
                

                // *** Compute derivatives, MSE
                let pix = refimg.get_pixel(x, y).0;
                let pixref: vec3 = vec3::new(pix[0] as f32 / 255., pix[1] as f32 / 255., pix[2] as f32 / 255.);
                let pixdif: vec3 = out_color - pixref;

                let pixmse = pixdif.dot(&pixdif); // sum of squares
                mse += pixmse;


                // *** pass derivative to each shape

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

                    shapes[i].backward(p, dldw, drgb);
                }
                
    
                let dldi: vec3 = pixdif * 2.;
                mesh.render(p, &mut mesh_color, true, dldi);
                unsafe {
                    if PRINT_NOW {
                        println!("dldi={:.4}", dldi);
                    }
                }
                
                unsafe {
                    PRINT_NOW = false;
                }

                
                out_color = out_color * 255.;
                *pixel = image::Rgb([
                    out_color.x as u8,
                    out_color.y as u8,
                    out_color.z as u8]);
            }
            mse = mse / imgx as f32 / imgy as f32;
            println!("MSE={:.3}", mse);

            let lr = 0.01;
            mesh.step(lr);
            for i in 0..nshapes {
                shapes[i].step(lr);
            }

            let filename: String = format!("anim/{}{:0>3}.jpg", save_path, _it);
            imgbuf.save(filename);
        }

        println!("Rendering took {:?}", start_time.elapsed());
        start_time = Instant::now();

        return;
    }

}

pub fn task_sdf_sampling(save_path: &str) {
    println!("\n=== Task SDF sampling ");
    let mut start_time = Instant::now();

    let yellow = vec3::new(255./255., 220./255., 3./255.);
    let red = vec3::new(255./255., 0., 0.);
    let black = vec3::new(0., 0., 0.);

    let imgx = 256;
    let imgy = imgx;
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);

    let th = 0.025;

    std::fs::create_dir("anim");

    // ******************** BACKWARD PASS

    let refimg = loadsdf::loadimage("resources/02_reference.png");

    // let mut circ = Shape {
    //     stype: ShapeType::Circle,
    //     C: Point::new(-0.2, 0.1),
    //     r: 0.2,
    //     th: th,
    //     color: red,
    //     no_smooth: true,
    //     ..Default::default()
    // };

    // let mut rect = Shape {
    //     stype: ShapeType::Rectangle,
    //     C: Point::new(0., 0.),
    //     w: 0.1,
    //     h: 0.1,
    //     th: th,
    //     color: yellow,
    //     no_smooth: true,
    //     ..Default::default()
    // };

    let mut circ = Shape {
        stype: ShapeType::Circle,
        C: Point::new(0.2, 0.2),
        r: 0.1,
        th: th,
        color: vec3::new(0.543, 0.2232, 0.42),
        no_smooth: true,
        ..Default::default()
    };

    let mut circ2 = Shape {
        stype: ShapeType::Circle,
        C: Point::new(-0.27, 0.22),
        r: 0.12,
        th: th,
        color: vec3::new(0.1, 0.6, 1.),
        no_smooth: true,
        ..Default::default()
    };

    let mut rect = Shape {
        stype: ShapeType::Rectangle,
        C: Point::new(-0.19, -0.22),
        w: 0.25,
        h: 0.25,
        th: th,
        color: vec3::new(0.7, 0.7, 0.),
        no_smooth: true,
        ..Default::default()
    };

    let mut shapes: Vec<Shape> = Vec::new();
    shapes.push(circ);
    shapes.push(circ2);
    shapes.push(rect);
    let nshapes = shapes.len();

    
    let mut smstep: Vec<f32> = vec![0.; nshapes]; // smstep = smoothstep(-sdf)
    let mut weight1: Vec<f32> = vec![0.; nshapes]; // weight1[i] * smstep[i] * color[i]
    // compute d pixel/d smoothstep using postfix sum
    let mut dstep_cum: Vec<vec3> = vec![black; nshapes];

    println!("Initialization took {:?}", start_time.elapsed());
    start_time = Instant::now();

    for _it in 0..150 {
        let mut mse: f32 = 0.;
        for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
            let xf = x as f32 / imgx as f32 - 0.5;
            let yf = y as f32 / imgy as f32 - 0.5;
            //let r = (xf / imgx as f32 * 128.0) as u8;
            //let b = (yf / imgy as f32 * 128.0) as u8;

            // again forward pass
            let p: Point = Point::new(xf, yf);
            let mut out_color = black;

            for i in 0..nshapes {
                let mut sdf = 0.;
                let nsamples = 5;
                for sample in 0..nsamples {
                    let x_add = (sample % 2) as f32;
                    let y_add = (sample / 2) as f32;
                    let p_add = Point::new(
                        rand::random::<f32>() / imgx as f32,
                        rand::random::<f32>() / imgy as f32);
                        // (0.25 - x_add * 0.5) / imgx as f32,
                        // (0.25 - y_add * 0.5) / imgy as f32);
                    let p_sample = p + p_add;
                    sdf += shapes[i].distance_alpha(p_sample);
                }

                smstep[i] = sdf / nsamples as f32;
                if i == 0 {
                    weight1[i] = 1.;
                } else {
                    weight1[i] = 1. - weight1[i - 1] * smstep[i - 1];
                }
                out_color = out_color + shapes[i].color * (weight1[i] * smstep[i]);
            }


            // *** Compute derivatives, MSE
            let pix = refimg.get_pixel(x, y).0;
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

                shapes[i].backward(p, dldw, drgb);
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
        let mut lr = 0.002;
        if _it > 30 {
            lr *= 2.;
        }
        for i in 0..nshapes {
            shapes[i].step(lr);
        }
        
        let filename: String = format!("anim/{}{:0>3}.jpg", save_path, _it);
        imgbuf.save(filename).unwrap();
    }

    println!("Rendering took {:?}", start_time.elapsed());
    start_time = Instant::now();
}

pub fn task_new_sdf(save_path: &str) {
    println!("\n=== Task SDF sampling ");
    let mut start_time = Instant::now();

    let yellow = vec3::new(255./255., 220./255., 3./255.);
    let red = vec3::new(255./255., 0., 0.);
    let black = vec3::new(0., 0., 0.);
    let sky = vec3::new(142., 202., 230.)/255.;
    let azul = vec3::new(0., 102., 230.)/255.;
    let fur = vec3::new(230., 230., 230.)/255.;

    let imgx = 256;
    let imgy = imgx;
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);

    let th = 0.025;

    std::fs::create_dir("anim");

    // ******************** BACKWARD PASS

    // let refimg = loadsdf::loadimage("resources/moon.jpg");
    let refimg = loadsdf::loadimage("resources/moon_stars.jpg");

    // let s: f32 = 0.5 / 7.;
    // let horse = Shape {
    //     stype: ShapeType::Horseshoe,
    //     C: Point::new(0., -2.5*s),
    //     rs: 2.5*s,
    //     rb: 3.5*s,
    //     th: th,
    //     color: sky,
    //     ..Default::default()
    // };

    // let moon = Shape {
    //     stype: ShapeType::Moon,
    //     C: Point::new(0., 2.5*s),
    //     r: 4.5*s,
    //     B: Point::new(0., 0.),
    //     rb: 3.5*s,
    //     th: th,
    //     color: yellow,
    //     ..Default::default()
    // };

    // let star1 = Shape {
    //     stype: ShapeType::Star5,
    //     C: Point::new(0., -1.5*s),
    //     rs: 0.08,
    //     rb: 0.15,
    //     th: th,
    //     color: azul,
    //     ..Default::default()
    // };

    let s: f32 = 0.5 / 7.;
    let horse = Shape {
        stype: ShapeType::Horseshoe,
        C: Point::new(0., -1.5*s),
        rs: 1.5*s,
        rb: 4.5*s,
        th: th,
        color: fur,
        ..Default::default()
    };

    let moon = Shape {
        stype: ShapeType::Moon,
        C: Point::new(0.3, 0.*s),
        r: 4.5*s,
        B: Point::new(0., 0.5*s),
        rb: 1.5*s,
        th: th,
        color: azul,
        ..Default::default()
    };

    let star1 = Shape {
        stype: ShapeType::Star5,
        C: Point::new(-0.1, -0.2),
        rs: 0.18,
        rb: 0.35,
        th: th,
        color: red,
        ..Default::default()
    };

    let mut shapes: Vec<Shape> = Vec::new();
    shapes.push(star1);
    shapes.push(horse);
    shapes.push(moon);
    
    
    
    let nshapes = shapes.len();

    
    let mut smstep: Vec<f32> = vec![0.; nshapes]; // smstep = smoothstep(-sdf)
    let mut weight1: Vec<f32> = vec![0.; nshapes]; // weight1[i] * smstep[i] * color[i]
    // compute d pixel/d smoothstep using postfix sum
    let mut dstep_cum: Vec<vec3> = vec![black; nshapes];

    println!("Initialization took {:?}", start_time.elapsed());
    start_time = Instant::now();

    for _it in 0..120 {
        let mut mse: f32 = 0.;
        for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
            let xf = x as f32 / imgx as f32 - 0.5;
            let yf = y as f32 / imgy as f32 - 0.5;
            //let r = (xf / imgx as f32 * 128.0) as u8;
            //let b = (yf / imgy as f32 * 128.0) as u8;

            // again forward pass
            let p: Point = Point::new(xf, yf);
            let mut out_color = black;

            for i in 0..nshapes {
                let mut sdf = 0.;
                let nsamples = 5;
                for sample in 0..nsamples {
                    let x_add = (sample % 2) as f32;
                    let y_add = (sample / 2) as f32;
                    let p_add = Point::new(
                        rand::random::<f32>() / imgx as f32,
                        rand::random::<f32>() / imgy as f32);
                        // (0.25 - x_add * 0.5) / imgx as f32,
                        // (0.25 - y_add * 0.5) / imgy as f32);
                    let p_sample = p + p_add;
                    sdf += shapes[i].distance_alpha(p_sample);
                }

                smstep[i] = sdf / nsamples as f32;
                if i == 0 {
                    weight1[i] = 1.;
                } else {
                    weight1[i] = 1. - weight1[i - 1] * smstep[i - 1];
                }
                out_color = out_color + shapes[i].color * (weight1[i] * smstep[i]);
            }


            // *** Compute derivatives, MSE
            let pix = refimg.get_pixel(x, y).0;
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

                shapes[i].backward(p, dldw, drgb);
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
        let mut lr = 0.002;
        for i in 0..nshapes {
            shapes[i].step(lr);
        }
        
        let filename: String = format!("anim/{}{:0>3}.jpg", save_path, _it);
        imgbuf.save(filename).unwrap();
    }

    println!("Rendering took {:?}", start_time.elapsed());
    start_time = Instant::now();
}


fn main() {
    // let guard = pprof::ProfilerGuardBuilder::default().frequency(1000).blocklist(&["libc", "libgcc", "pthread", "vdso"]).build().unwrap();
    
    // task1("fractal");
    // task_bitmap("fractal");
    // task_edge_sampling("fractal");
    // task_sdf_sampling("fractal");
    task_new_sdf("fractal");

    // if let Ok(report) = guard.report().build() {
    //     let file = File::create("flamegraph.svg").unwrap();
    //     let mut options = pprof::flamegraph::Options::default();
    //     options.image_width = Some(2500);
    //     report.flamegraph_with_options(file, &mut options).unwrap();
    // };

    
    return;
}
