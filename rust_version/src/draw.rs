use crate::palettes;
use image::{
    ImageBuffer, 
    RgbImage, DynamicImage,
};

pub fn plot_2d(map: crate::GridValues, width: usize, height: usize, flnm: &str) {
    let mut img: RgbImage = ImageBuffer::new(map.nx as u32, map.ny as u32);
    let mut idx: usize;

    for jx in 0..map.nx {
        for jy in 0..map.ny {
            idx = jx * map.ny + jy;
            let u = map.uxy[idx];
            let col = palette_rev(u, map.umin, map.umax);
            img.put_pixel(jx as u32, jy as u32, image::Rgb(col));
        }
    }

    let mut img2 = DynamicImage::ImageRgb8(img);
    img2 = img2.resize(width as u32, height as u32, image::imageops::FilterType::Nearest);

    img2.save(&flnm).unwrap();
}

fn _palette(f: f64, f_min: f64, f_max: f64) -> [u8; 3] {
    let z = ((f - f_min) / (f_max - f_min) * 99.0) as usize;

    palettes::GRAY_PAL[z]
}

fn palette_rev(f: f64, f_min: f64, f_max: f64) -> [u8; 3] {
    let z = ((f - f_min) / (f_max - f_min) * 99.0) as usize;

    palettes::GRAY_PAL[99-z]
}
