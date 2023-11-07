use crate::palettes;
use image::{
    ImageBuffer, 
    RgbImage,
};

pub fn plot_2d(map: crate::GridValues, flnm: &str) {
    let mut img: RgbImage = ImageBuffer::new(map.nx as u32, map.ny as u32);

    for jx in 0..map.nx {
        for jy in 0..map.ny {
            let u = map.uxy[jx * map.ny + jy];
            let col = palette(u, map.umin, map.umax);
            img.put_pixel(jx as u32, jy as u32, image::Rgb(col));
        }
    }

    img.save(&flnm).unwrap();
}

fn palette(f: f64, f_min: f64, f_max: f64) -> [u8; 3] {
    let z = ((f - f_min) / (f_max - f_min) * 99.0) as usize;

    palettes::RDBU_PAL[z]
}
