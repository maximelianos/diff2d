use std::borrow::Borrow;
use std::{ops, collections::VecDeque};
use std::cell::{RefCell, RefMut, Ref, Cell};
use std::rc::Rc; // immutable! https://doc.rust-lang.ru/stable/rust-by-example/std/rc.html
use std::collections::{BTreeSet, LinkedList};


#[derive(Default, Debug, Clone, Copy)]
enum Dtype {
    #[default]
    SCALAR,
    MATRIX
}

#[derive(Default, Debug, Clone, Copy, PartialEq)]
enum OpType {
    #[default]
    Const,
    Var,
    ScPlus,
    ScSub,
    ScMul,
    ScDiv,
    ScNeg,
    ScSqrt,
    ScSquare,
    ScSmoothstep,
    ScAbs
}

#[derive(Default, Debug, Clone, Copy)]
pub struct D {
    dtype: Dtype,
    pub scalar: f32,
    pub scalar_d: f32,
    id: usize,
    order: usize,
    optype: OpType,
    inp: [usize; 2],

    smleft: f32, // smoothstep
    smright: f32,
}

#[derive(Default, Debug, Clone, Copy)]
pub struct Dp {
    //pub Rc<RefCell<D>> // borrow here breaks everything
    pub id: usize,
    pub print: bool
    // g: Graph
} // Rc is for multiple references (but immutable), RefCell is for mutability

impl Dp {
    pub fn new_f(scalar: f32, is_const: bool) -> Dp {
        // New variable without associated operation
        let graph: Graph = GRAPH();
        let mut g: RefMut<DGraph> = graph.as_ref().borrow_mut();
        let node = D {
            dtype: Dtype::SCALAR,
            scalar: scalar,
            optype: if is_const { OpType::Const } else { OpType::Var },
            ..Default::default()
        };

        return g.push(node);
        
        // println!("New variable, len={}", g.objects.len());
    }

    pub fn var_f(scalar: f32) -> Dp {
        return Dp::new_f(scalar, false);
    }

    pub fn const_f(scalar: f32) -> Dp {
        return Dp::new_f(scalar, true);
    }

    fn op_float(lhs: &Dp, rhs: &Dp, op: OpType) -> Dp {
        // Create result of operation
        let graph: Graph = GRAPH();
        let mut g: RefMut<DGraph> = graph.as_ref().borrow_mut();

        let a = &g.objects[lhs.id];
        let b = &g.objects[rhs.id];
        let res: f32;
        use OpType::*;
        match op {
            ScPlus => { res = a.scalar + b.scalar; },
            ScSub => { res = a.scalar - b.scalar; },
            ScMul => { res = a.scalar * b.scalar; },
            ScDiv => { res = a.scalar / b.scalar; },
            ScNeg => { res = -a.scalar; },
            ScSqrt => { res = a.scalar.powf(0.5); },
            ScSquare => { res = a.scalar.powf(2.); },
            ScAbs => { res = a.scalar.abs() },
            _ => {
                panic!("Unimplemented op");
            }
        }

        let c = D {
            dtype: Dtype::SCALAR,
            scalar: res,
            optype: op,
            inp: [a.id, b.id],
            ..Default::default()
        };
        
        return g.push(c);
        // println!("New result, len={}", g.objects.len());
    }

    pub fn obj(&self) -> D {
        let graph: Graph = GRAPH();
        let mut g: RefMut<DGraph> = graph.as_ref().borrow_mut();
        return g.objects[self.id];
    }

    pub fn s(&self) -> f32 {
        let graph: Graph = GRAPH();
        let mut g: RefMut<DGraph> = graph.as_ref().borrow_mut();
        return g.objects[self.id].scalar;
    }

    pub fn set(&self, scalar: f32, scalar_d: f32) {
        let graph: Graph = GRAPH();
        let mut g: RefMut<DGraph> = graph.as_ref().borrow_mut();
        g.objects[self.id].scalar = scalar;
        g.objects[self.id].scalar_d = scalar_d;
    }


    pub fn backward(&self) {
        // dL / dx_i
        let graph: Graph = GRAPH();
        let mut g: RefMut<DGraph> = graph.as_ref().borrow_mut();
        let cur_id: usize = self.id;
        // let objects: &mut Vec<D> = &mut g.objects;
        
        g.objects[cur_id].scalar_d = 1.;
        
        // "almost" BFS traverse the differentiation graph
        let objects_len = g.objects.len();
        g.visited.clear();
        g.visited.resize(objects_len, false);
        g.visited[cur_id] = true;
        // let mut deque: VecDeque<usize> = VecDeque::new();
        // deque.push_back(cur_id);

        // *** always select object with maximum order num
        let cur_order = g.objects[cur_id].order;
        g.bfs_order.insert( (cur_order, cur_id) );
        let mut _it = 0;
        while !g.bfs_order.is_empty() {
            _it += 1;
            // *** current node
            let (_, cur_id) = g.bfs_order.pop_last().unwrap(); //deque.pop_front().unwrap();
            let a = g.objects[cur_id];

            let i1 = a.inp[0];
            let i2 = a.inp[1];
            let b1 = g.objects[i1];
            let b2 = g.objects[i2];
            use OpType::*;
            match a.optype {
                Const => (),
                Var => (),
                ScPlus => {
                    // println!("backward plus {} {}", i1, i2);
                    // *** previous nodes
                    // same indexes are fine here!
                    g.objects[i1].scalar_d += a.scalar_d;
                    g.objects[i2].scalar_d += a.scalar_d;
                },
                ScSub => {
                    g.objects[i1].scalar_d += a.scalar_d;
                    g.objects[i2].scalar_d += -a.scalar_d;
                },
                ScMul => {
                    // println!("backward mul {} {}", i1, i2);
                    g.objects[i1].scalar_d += a.scalar_d * b2.scalar;
                    g.objects[i2].scalar_d += a.scalar_d * b1.scalar;
                },
                ScSqrt => {
                    // println!("backward sqrt {}", i1);
                    g.objects[i1].scalar_d += a.scalar_d / (2. * b1.scalar.sqrt());
                },
                ScSquare => {
                    g.objects[i1].scalar_d += a.scalar_d * 2. * b1.scalar;
                },
                ScNeg => {
                    g.objects[i1].scalar_d += -a.scalar_d;
                },
                ScDiv => {
                    g.objects[i1].scalar_d += a.scalar_d / b2.scalar;
                    g.objects[i2].scalar_d += a.scalar_d * (-b1.scalar) / (b2.scalar*b2.scalar);
                },
                ScSmoothstep => {
                    let x = b1.scalar;
                    let d_b1: f32;
                    if x < a.smleft { d_b1 = 0. }
                    else if x > a.smright { d_b1 = 0. }
                    else {
                        let y = (x - a.smleft) / (a.smright - a.smleft);
                        d_b1 = 6. * y * (1. - y) / (a.smright - a.smleft);
                    }
                    
                    g.objects[i1].scalar_d += a.scalar_d * d_b1;
                },
                ScAbs => {
                    g.objects[i1].scalar_d += if b1.scalar > 0. { a.scalar_d } else { -a.scalar_d };
                },
                _ => {
                    panic!("Unimplemented operation {:?}", a.optype);
                }
            }
            match a.optype {
                Const => {
                    g.free_node(a.id);
                },
                Var => (),
                _ => { 
                    // Visited is not duplicated when i1 == i2 :)
                    // Order is VERY important, wrong order leads to wrong deriv
                    if !g.visited[i1] { 
                        let i1_order = g.objects[i1].order;
                        g.bfs_order.insert( (i1_order, i1) );
                        g.visited[i1] = true;
                    }
                    if !g.visited[i2] { 
                        let i2_order = g.objects[i2].order;
                        g.bfs_order.insert( (i2_order, i2) );
                        g.visited[i2] = true;
                    }
                    
                    g.free_node(a.id);
                }
            }
        }
        if self.print {
            println!("Nodes visited in backward {}", _it);
        }
    }

    pub fn sqrt(&self) -> Dp {
        // println!("> Dp.sqrt() was called");
        return Dp::op_float(&self, &self, OpType::ScSqrt);
    }

    pub fn square(&self) -> Dp {
        // println!("> Dp.sqrt() was called");
        return Dp::op_float(&self, &self, OpType::ScSquare);
    }

    pub fn smoothstep(&self, left: f32, right: f32) -> Dp {
        let graph: Graph = GRAPH();
        let mut g: RefMut<DGraph> = graph.as_ref().borrow_mut();

        let a = &g.objects[self.id];

        let x: f32 = a.scalar;
        let res: f32;
        if x < left {
            res = 0.;
        } else if x > right {
            res = 1.;
        } else {
            let y = (x - left) / (right - left);
            res = y * y * (3. - 2. * y);
        }

        let c = D {
            dtype: Dtype::SCALAR,
            scalar: res,
            optype: OpType::ScSmoothstep,
            inp: [a.id, a.id],
            smleft: left,
            smright: right,
            ..Default::default()
        };
        
        return g.push(c);
    }

    pub fn abs(&self) -> Dp {
        return Dp::op_float(&self, &self, OpType::ScAbs);
    }
}


// https://doc.rust-lang.org/core/ops/
impl ops::Add<Dp> for Dp {
    type Output = Dp;

    fn add(self, _rhs: Dp) -> Dp {
        // println!("> Dp.add(Dp) was called");
        return Dp::op_float(&self, &_rhs, OpType::ScPlus);
    }
}

impl ops::Sub<Dp> for Dp {
    type Output = Dp;

    fn sub(self, _rhs: Dp) -> Dp {
        return Dp::op_float(&self, &_rhs, OpType::ScSub);
    }
}

impl ops::Mul<Dp> for Dp {
    type Output = Dp;

    fn mul(self, _rhs: Dp) -> Dp {
        return Dp::op_float(&self, &_rhs, OpType::ScMul);
    }
}

impl ops::Neg for Dp {
    type Output = Dp;

    fn neg(self) -> Dp {
        return Dp::op_float(&self, &self, OpType::ScNeg);
    }
}

impl ops::Div<Dp> for Dp {
    type Output = Dp;

    fn div(self, _rhs: Dp) -> Dp {
        return Dp::op_float(&self, &_rhs, OpType::ScDiv);
    }
}


#[derive(Default, Debug)]
pub struct DGraph {
    pub objects: Vec<D>,
    free_nodes: LinkedList<usize>,
    is_free: Vec<bool>,
    tmp_nodes: LinkedList<usize>,

    // --> temporary for backward
    visited: Vec<bool>, 
    bfs_order: BTreeSet<(usize, usize)>,
    // <--

    order: usize,
    pub eval: bool
}

impl DGraph {
    pub fn free_cnt(&self) -> usize {
        return self.free_nodes.len();
    }

    fn push(&mut self, mut node: D) -> Dp {
        node.order = self.order;
        self.order += 1;

        let id;
        if self.free_nodes.is_empty() {
            // push to vec
            id = self.objects.len();
            self.is_free.push(false);
            node.id = id;
            self.objects.push(node);
        } else {
            // assign to some vec elem
            id = self.free_nodes.pop_back().unwrap();
            self.is_free[id] = false;
            node.id = id;
            self.objects[id] = node;
        }

        if node.optype != OpType::Var { self.tmp_nodes.push_back(id); }
        // println!("Create node id {}", id);
        return Dp { id: id, ..Default::default() };
    }

    fn free_node(&mut self, id: usize) {
        // println!("Free node id {}", node.id);
        if !self.is_free[id] { // the vector always keeps right info
            self.free_nodes.push_back(id);
            self.is_free[id] = true;
        }
    }

    pub fn clean_graph(&mut self) {
        while !self.tmp_nodes.is_empty() {
            let id = self.tmp_nodes.pop_back().unwrap();
            self.free_node(id);
        }
    }
}

pub type Graph = Rc<RefCell<DGraph>>;

thread_local!(static GLOBAL_DATA: Graph = Default::default());

pub fn GRAPH() -> Graph {
    let mut ptr: Graph = Default::default();
    GLOBAL_DATA.with(|data| {
        ptr = data.clone();
    });
    return ptr;
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
    // let z = -x;
    // let z1 = (z+(z+(z+z)));
    // z1.backward();
    // {
    //     let graph = GRAPH();
    //     let g = graph.as_ref().borrow();
    //     let objects = &g.objects;
    //     println!("Graph length {:?} free {:?}", objects.len(), g.free_cnt());
    //     println!("z2={}", objects[z1.id].scalar);
    //     println!("dx={}", objects[x.id].scalar_d);
    //     println!("dy={}", objects[y.id].scalar_d);
    // }
    

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





// ************************ Point with differentiable coordinates

#[derive(Default, Debug, Clone, Copy)]
struct Point {
    x: Dp,
    y: Dp
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

    fn len(&self) -> Dp {
        // sqrt here causes numerical instability, clip
        let r = self.x * self.x + self.y * self.y;
        if r.obj().scalar < 0.001 { return Dp::const_f(0.); }
        else { return r.sqrt(); }
    }

    fn dot(&self, b: &Point) -> Dp {
        return self.x * b.x + self.y * b.y;
    }

    fn cross(&self, b: &Point) -> Dp {
        return self.x * b.y - b.x * self.y;
    }

    fn proj(&self, b: &Point) -> Dp {
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