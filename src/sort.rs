fn get_sorted(v: &[f32]) -> Vec<usize> {
    let mut idx: Vec<usize> = vec![0; v.len()];
    for i in 0..v.len() {
        idx[i] = i;
    }
    radsort::sort_by_key(&mut idx, |i| v[*i]);
    return idx;
}

fn main() {
    let mut v: Vec<f32> = vec![3., 2., 1., 4., -0.5, 4.];
    let idx = get_sorted(&v);
    for i in 0..v.len() {
        print!("{} ", v[idx[i]]);
    }
}