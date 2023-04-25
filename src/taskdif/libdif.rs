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

#[derive(Default, Debug, Clone, Copy)]
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
        while !g.bfs_order.is_empty() {
            // *** current node
            let (_, cur_id) = g.bfs_order.pop_last().unwrap(); //deque.pop_front().unwrap();
            let a = g.objects[cur_id];

            let i1 = a.inp[0];
            let i2 = a.inp[1];
            let b1 = g.objects[i1];
            let b2 = g.objects[i2];
            use OpType::*;
            if !g.freeze { // unimportant
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
                    _ => {
                        panic!("Unimplemented operation {:?}", a.optype);
                    }
                }
            }
            match a.optype {
                Const => {
                    g.free_node(&a);
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
                    
                    g.free_node(&a);
                }
            }
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

    // --> temporary for backward
    visited: Vec<bool>, 
    bfs_order: BTreeSet<(usize, usize)>,
    // <--

    order: usize,
    pub freeze: bool // only free nodes during backward
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
        // println!("Create node id {}", id);
        return Dp { id: id };
    }

    fn free_node(&mut self, node: &D) {
        // println!("Free node id {}", node.id);
        if !self.is_free[node.id] { // the vector always keeps right info
            self.free_nodes.push_back(node.id);
            self.is_free[node.id] = true;
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
