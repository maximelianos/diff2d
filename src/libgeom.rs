use std::ops;

#[derive(Default, Copy, Clone, Debug, PartialEq)]
pub struct Point {
    pub x: f32,
    pub y: f32
}

impl Point {
    pub fn new(x: f32, y: f32) -> Point {
        return Point{x: x, y: y};
    }

    pub fn len(&self) -> f32 {
        return (self.x * self.x + self.y * self.y).sqrt();
    }

    pub fn dot(&self, b: Point) -> f32 {
        return self.x * b.x + self.y * b.y;
    }

    pub fn cross(&self, b: Point) -> f32 {
        return self.x * b.y - b.x * self.y;
    }

    pub fn proj(&self, b: Point) -> f32 {
        return self.dot(b) / self.len();
    }

    pub fn rotate(&self, alpha: f32) -> Point {
        return Point::new(
            self.x * alpha.cos() - self.y * alpha.sin(),
            self.x * alpha.sin() + self.y * alpha.cos());
    }
}

impl ops::Add<Point> for Point {
    type Output = Point;

    fn add(self, _rhs: Point) -> Point {
        //println!("> Point.add(Point) was called");
        return Point {x: self.x + _rhs.x, y: self.y + _rhs.y};
    }
}

impl ops::Mul<f32> for Point {
    type Output = Point;

    fn mul(self, _rhs: f32) -> Point {
        return Point {x: self.x * _rhs, y: self.y * _rhs};
    }
}

impl ops::Neg for Point {
    type Output = Self;

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
