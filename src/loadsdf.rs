use image::{GrayImage, DynamicImage, GenericImageView};

pub fn loadsdf(path: &str) -> GrayImage {
    // Use the open function to load an image from a Path.
    // `open` returns a `DynamicImage` on success.
    // Image color type: Luma8, range 0-255
    let img = image::open(path).unwrap();
    let grayscale = img.to_luma8();

    // *** The dimensions method returns the images width and height.
    // println!("dimensions {:?}", img.dimensions());
    // let (w, h) = img.dimensions();

    // *** The color method returns the image's `ColorType`.
    // println!("{:?}", img.color());

    // *** Probe pixel values
    // for (x, y, pixel) in grayscale.enumerate_pixels() {
    //     if (x == 0 || x == w-1) && (y == 0 || y == h-1) {
    //         let a = pixel.0;
    //         println!("pixel coors x={x} y={y} pixel={a:?}");
    //     }
    // }

    // return ImageBuffer, to get pixels use smth like imgbuf.get_pixel(x, y)
    return grayscale;
}

pub fn loadimage(path: &str) -> DynamicImage {
    // Use the open function to load an image from a Path.
    // `open` returns a `DynamicImage` on success.
    let img = image::open(path).unwrap();

    // The dimensions method returns the images width and height.
    println!("dimensions {:?}", img.dimensions());

    // The color method returns the image's `ColorType`.
    println!("{:?}", img.color());

    return img;
}