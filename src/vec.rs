fn main() {
    let mut arr = Vec::new();
    arr.push(1);
    arr.push(2);
    arr.push(30);
    println!("arr[0]={}", arr[0]);
    println!("arr[3]={}", arr[3]); // bound checking!
}