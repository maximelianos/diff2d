use std::collections::HashMap;

fn main() {
    // *** How to use HashMap
    let mut books = HashMap::new();
    books.insert(1, "Voyna i Mir".to_string());
    books.insert(10, "Crime and Punishment".to_string());
    println!("{}", books.contains_key(&10));
    println!("{}", books.get(&1).unwrap());
    for (key, val) in books.iter() {
        println!("key: {key} val: {val}");
    }
    println!("{}", books[&2]);
}