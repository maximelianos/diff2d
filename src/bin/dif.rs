use std::{ops, collections::VecDeque};
use std::cell::{RefCell, RefMut, Ref};
use std::rc::Rc; // immutable! https://doc.rust-lang.ru/stable/rust-by-example/std/rc.html

use grid::Grid;

#[derive(Default, Debug)]
enum Dtype {
    #[default]
    SCALAR,
    MATRIX
}

#[derive(Default, Debug)]
enum OpType {
    #[default]
    NONE,
    SC_PLUS,
    SC_MINUS,
    SC_MUL,
    SC_SQRT
}

#[derive(Default, Debug)]
struct D {
    dtype: Dtype,
    scalar: f32,
    //matrix: Grid<(f32, f32)>,
    scalar_d: f32,
    //matrix_d: Grid<(f32, f32)>,
    id: usize,
    optype: OpType,
    inp: [usize; 2],

    g: Graph
}

#[derive(Default, Debug, Clone)]
struct Dp(Rc<RefCell<D>>); // Rc is for multiple references (but immutable), RefCell is for mutability

impl Dp {
    fn new_float(scalar: f32, graph: &Graph) -> Dp {
        // New variable without associated operation
        let mut g: RefMut<DGraph> = graph.as_ref().borrow_mut();
        let dp = Dp(Rc::new(RefCell::new( D {
            dtype: Dtype::SCALAR,
            scalar: scalar,
            id: g.objects.len() as usize,
            optype: OpType::NONE,
            g: graph.clone(),
            ..Default::default()
        } )));
        
        g.objects.push(dp.clone());
        // println!("New variable, len={}", g.objects.len());
        return dp;
    }

    fn float(base: D, graph: &Graph) -> Dp {
        // Create result of operation
        let mut g: RefMut<DGraph> = graph.as_ref().borrow_mut();
        let dp = Dp(Rc::new(RefCell::new(D {
            id: g.objects.len() as usize,
            g: graph.clone(),
            ..base
        } )));
        
        g.objects.push(dp.clone());
        // println!("New result, len={}", g.objects.len());
        return dp;
    }



    fn backward(&self) {
        // dL / dx_i
        let gref: Graph;
        let cur_id: usize;
        {
            let mut a = self.0.as_ref().borrow_mut();
            a.scalar_d = 1.;
            gref = a.g.clone();
            cur_id = a.id;
            // now a is free
        }
        let g: RefMut<DGraph> = gref.as_ref().borrow_mut();

        // BFS traverse the differentiation graph
        let mut visited: Vec<bool> = vec![false; g.objects.len()];
        let mut deque: VecDeque<usize> = VecDeque::new();
        deque.push_back(cur_id);
        while !deque.is_empty() {
            let cur_id: usize = deque.pop_front().unwrap();
            let a = g.objects[cur_id].0.as_ref().borrow();
            if visited[a.id] {
                continue;
            }
            visited[a.id] = true;

            use OpType::*;
            let i1 = a.inp[0];
            let i2 = a.inp[1];
            match a.optype {
                NONE => (),
                SC_PLUS => {
                    // println!("backward plus {} {}", i1, i2);
                    let mut b1 = g.objects[i1].0.as_ref().borrow_mut();
                    b1.scalar_d += a.scalar_d;
                    let mut b2 = g.objects[i2].0.as_ref().borrow_mut();
                    b2.scalar_d += a.scalar_d;
                    deque.push_back(b1.id);
                    deque.push_back(b2.id);

                },
                SC_MUL => {
                    // println!("backward mul {} {}", i1, i2);
                    let mut b1 = g.objects[i1].0.as_ref().borrow_mut();
                    let mut b2 = g.objects[i2].0.as_ref().borrow_mut();
                    b1.scalar_d += a.scalar_d * b2.scalar;
                    b2.scalar_d += a.scalar_d * b1.scalar;
                    deque.push_back(b1.id);
                    deque.push_back(b2.id);
                },
                SC_SQRT => {
                    // println!("backward sqrt {}", i1);
                    let mut b1 = g.objects[i1].0.as_ref().borrow_mut();
                    b1.scalar_d += a.scalar_d / (2. * b1.scalar.sqrt());
                    deque.push_back(b1.id);
                }
                _ => ()
            }
        }
    }


}

impl ops::Add<Dp> for Dp {
    type Output = Dp;

    fn add(self, _rhs: Dp) -> Dp {
        // println!("> Dp.add(Dp) was called");
        let a: Ref<D> = self.0.as_ref().borrow();
        let b: Ref<D> = _rhs.0.as_ref().borrow();

        let c = D {
            dtype: Dtype::SCALAR,
            scalar: a.scalar + b.scalar,
            optype: OpType::SC_PLUS,
            inp: [a.id, b.id],
            ..Default::default()
        };
        let cp = Dp::float(c, &a.g);

        return cp;
    }
}

impl ops::Mul<&Dp> for &Dp {
    type Output = Dp;

    fn mul(self, _rhs: &Dp) -> Dp {
        // println!("> Dp.mul(Dp) was called");
        let a: Ref<D> = self.0.as_ref().borrow();
        let b: Ref<D> = _rhs.0.as_ref().borrow();

        let c = D {
            dtype: Dtype::SCALAR,
            scalar: a.scalar * b.scalar,
            optype: OpType::SC_MUL,
            inp: [a.id, b.id],
            ..Default::default()
        };
        let cp = Dp::float(c, &a.g);

        return cp;
    }
}

impl Dp {
    fn sqrt(&self) -> Dp {
        // println!("> Dp.sqrt() was called");
        let a: Ref<D> = self.0.as_ref().borrow();

        let c = D {
            dtype: Dtype::SCALAR,
            scalar: a.scalar.sqrt(),
            optype: OpType::SC_SQRT,
            inp: [a.id, 0],
            ..Default::default()
        };
        let cp = Dp::float(c, &a.g);

        return cp;
    }
}


#[derive(Default, Debug)]
struct DGraph {
    objects: Vec<Dp>
}

type Graph = Rc<RefCell<DGraph>>;



// #[thread_local]
// static mut GRAPH: DGraph = Default::default();
//thread_local!(static GRAPH: RefCell<DGraph> = RefCell::new(Default::default()));


fn main() {
    // Statements here are executed when the compiled binary is called

    // Print text to the console
    //let n: i8 = 257; // rust checks range of number and compilation fails :)
    //println!("Hello World! x = {:>5} {}", n); // println checks number of arguments :)

    // println!("{:?}", GRAPH.borrow_mut());
    // let GRAPH: RefCell<DGraph> = RefCell::new(Default::default());
    // println!("{}", GRAPH);
    let graph: Graph = Rc::new(RefCell::new(DGraph::default()));

    let x: Dp = Dp::new_float(2., &graph);
    let y: Dp = Dp::new_float(3., &graph);
    // let w: Dp = Dp::new_float(1., graph.clone());

    let mut z = &x * &y;
    for i in 0..1000000 {
        z = z.clone() + x.clone();
        //z = z.sqrt();
    }
    

    z.backward();

    use std::io::{stdin,stdout,Write};
    let mut s=String::new();
    print!("Please enter some text: ");
    let _=stdout().flush();
    stdin().read_line(&mut s).expect("Did not enter a correct string");
    if let Some('\n')=s.chars().next_back() {
        s.pop();
    }
    if let Some('\r')=s.chars().next_back() {
        s.pop();
    }
    println!("You typed: {}",s);


    println!("Graph length {:?}", graph.as_ref().borrow().objects.len());
    // println!("z1={} dz1={}", z1.0.as_ref().borrow().scalar, z1.0.as_ref().borrow().scalar_d);

    // println!("z={} dz={}", z.0.as_ref().borrow().scalar, z.0.as_ref().borrow().scalar_d);

    
    // println!("dz={}", z.0.as_ref().borrow().scalar_d);

    println!("dx={}", x.0.as_ref().borrow().scalar_d);
    println!("dy={}", y.0.as_ref().borrow().scalar_d);
    
    let msg = "Hello World!";
    println!("{}", msg);
}