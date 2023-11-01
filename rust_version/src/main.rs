// Program to simulate movement of classical electrons through a 2D array of nanostructures (quantum dots)
// More comments forthcoming
use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;
use rand::Rng;

//Reciprocal electron mass in nm^2/fs^2/eV
const RME: f64 = 0.176;

// Input parameters
const ME: f64 = 0.067;
const TAU: f64 = 4000.0;
const U0: f64 = -0.2;
const A: f64 = 20.0;
const B: f64 = 18.0;
const L: f64 = 1500.0;
const NN: usize = 150;
const DT: f64 = 1.0;
const EE: f64 = 0.0001;
const EB: f64 = 0.0;
const NE: usize = 100;
const NC1: usize = 10;

// Auxillary parameters
const NC: usize = NC1 * NC1;
const TF: f64 = 6.0 * TAU;
const NT: usize = (TF / DT) as usize;
const CM: f64 = RME / ME;
const CTAU: f64 = 1.0 / TAU;
const CA: f64 = 0.5 / A / A;
const CB: f64 = 0.5 / B / B;
const CWX: f64 = 2.0 * CM * CA * U0;
const CWY: f64 = 2.0 * CM * CB * U0;
const DC: f64 = L / (NC1 as f64);
const W: f64 = 3.0 * A;

fn main() {
    // measure time
    let now = Instant::now();

    //make directory
    let dir_name = "results";
    fs::create_dir_all(dir_name).expect("Error creating directory");

    // random number generator
    let mut rng = rand::thread_rng();

    // coordinates of nanostructures, cell numbers of nanostructures
    let mut xn: [f64; NN] = [0.0; NN];
    let mut yn: [f64; NN] = [0.0; NN];
    let mut cn: [usize; NN] = [1; NN];
    
    // initiate random coordinates, but not too close to the boundary
    for jn in 1..(NN+1) {
        let rng_x: f64 = rng.gen();
        let rng_y: f64 = rng.gen();
        xn[jn-1] = W + rng_x * (L - 2.0 * W);
        yn[jn-1] = W + rng_y * (L - 2.0 * W);
    }

    //file for nanostructures data
    let fl_name = "nanostructures.dat";
        let file_path: PathBuf = [dir_name, fl_name].iter().collect();
        let mut my_file = fs::File::create(file_path).expect("Error creating file");

    // define cell number for each nanostructure, write data to file
    for jn in 1..(NN+1) {
        cn[jn-1] = cell_numbers(xn[jn-1], yn[jn-1])[0];
        writeln!(my_file, "{} {} {}", xn[jn-1], yn[jn-1], cn[jn-1]).expect("Error writing to file");
    }

    //file for cells data
    let fl_name = "cells.dat";
        let file_path: PathBuf = [dir_name, fl_name].iter().collect();
        let mut my_file = fs::File::create(file_path).expect("Error creating file");

    let mut nano_in_cell: [usize; NC] = [0; NC];
    let mut nano_list: [[(f64, f64); NN]; NC] = [[(0.0, 0.0); NN]; NC];

    for jn in 1..(NN+1) {
        let jc = cn[jn-1];
        let nc = nano_in_cell[jc-1];
        nano_in_cell[jc-1] = nc + 1;
        nano_list[jc-1][nc].0 = xn[jn-1];
        nano_list[jc-1][nc].1 = yn[jn-1];
    }

    for jc in 1..(NC+1) {
        writeln!(my_file, "{}", nano_in_cell[jc-1]).expect("Error writing to file");
    }

    // Start calculations
    for je in 1..(NE+1) {
        //file for trajectory data for a single electron
        let fl_name = format!("Trajectory_{}.dat", je);
            let file_path: PathBuf = [dir_name, &fl_name].iter().collect();
            let mut my_file = fs::File::create(file_path).expect("Error creating file");

        // initial data
        let rnx: f64 = rng.gen();
        let rny: f64 = rng.gen();
        let rnvx: f64 = rng.gen();
        let rnvy: f64 = rng.gen();
        let mut x0: f64= W + (L - 2.0 * W) * rnx;
        let mut y0: f64 = W + (L - 2.0 * W) * rny;
        let mut vx0: f64 = 0.01 * rnvx;
        let mut vy0: f64 = 0.01 * rnvy;
        
        for jt in 1..NT {
            let t: f64 = DT * (jt as f64);
            writeln!(my_file, "{} {} {} {} {}", t, x0, y0, vx0, vy0).expect("Error writing to file");

            // Runge-Kutta method start
            let rk_vx1 = fx(x0, y0, vx0, vy0, &nano_in_cell, &nano_list);
            let rk_vy1 = fy(x0, y0, vx0, vy0, &nano_in_cell, &nano_list);
            let rk_x1 = vx0;
            let rk_y1 = vy0;

            let x1 = x0 + rk_x1 * DT * 0.5;
            let y1 = y0 + rk_y1 * DT * 0.5;
            let vx1 = vx0 + rk_vx1 * DT * 0.5;
            let vy1 = vy0 + rk_vy1 * DT * 0.5;

            let rk_vx2 = fx(x1, y1, vx1, vy1, &nano_in_cell, &nano_list);
            let rk_vy2 = fy(x1, y1, vx1, vy1, &nano_in_cell, &nano_list);
            let rk_x2 = vx1;
            let rk_y2 = vy1;

            let x1 = x0 + rk_x2 * DT * 0.5;
            let y1 = y0 + rk_y2 * DT * 0.5;
            let vx1 = vx0 + rk_vx2 * DT * 0.5;
            let vy1 = vy0 + rk_vy2 * DT * 0.5;

            let rk_vx3 = fx(x1, y1, vx1, vy1, &nano_in_cell, &nano_list);
            let rk_vy3 = fy(x1, y1, vx1, vy1, &nano_in_cell, &nano_list);
            let rk_x3 = vx1;
            let rk_y3 = vy1;

            let x1 = x0 + rk_x3 * DT;
            let y1 = y0 + rk_y3 * DT;
            let vx1 = vx0 + rk_vx3 * DT;
            let vy1 = vy0 + rk_vy3 * DT;

            let rk_vx4 = fx(x1, y1, vx1, vy1, &nano_in_cell, &nano_list);
            let rk_vy4 = fy(x1, y1, vx1, vy1, &nano_in_cell, &nano_list);
            let rk_x4 = vx1;
            let rk_y4 = vy1;

            x0 = x0 + (rk_x1 + 2.0 * rk_x2 + 2.0 * rk_x3 + rk_x4) * DT / 6.0;
            y0 = y0 + (rk_y1 + 2.0 * rk_y2 + 2.0 * rk_y3 + rk_y4) * DT / 6.0;
            vx0 = vx0 + (rk_vx1 + 2.0 * rk_vx2 + 2.0 * rk_vx3 + rk_vx4) * DT / 6.0;
            vy0 = vy0 + (rk_vy1 + 2.0 * rk_vy2 + 2.0 * rk_vy3 + rk_vy4) * DT / 6.0;

            // Periodic boundary conditions
            if x0 > L {x0 = x0 - L};
            if x0 < 0.0 {x0 = x0 + L};
            if y0 > L {y0 = y0 - L};
            if y0 < 0.0 {y0 = y0 + L};
        }

            
    }


    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
}

// function to define a singe nanostructure
fn u1(x: f64, y: f64, xn: f64, yn: f64) -> f64 {
    if (x - xn).abs() < (5.0 * A) && (y - yn).abs() < (5.0 * B) {
        (- CA * (x - xn) * (x - xn) - CB * (y - yn) * (y - yn)).exp()
    } else { 
        0.0
    }
}
// accelerations due to a single nanostructure
fn wx1(x: f64, y: f64, xn: f64, yn: f64) -> f64 {
    (x - xn) * u1(x, y, xn, yn)
}
fn wy1(x: f64, y: f64, xn: f64, yn: f64) -> f64 {
    (y - yn) * u1(x, y, xn, yn)
}

// function to check cell number for any pair of coordinates, as well as its 8 neighbors
fn cell_numbers(x: f64, y: f64) -> [usize; 9] {
    let mut cx = (x / DC).floor() as usize + 1;
    let mut cy = (y / DC).floor() as usize + 1;
    if cx > NC1 {cx = NC1};
    if cy > NC1 {cy = NC1};
    if cx < 1 {cx = 1};
    if cy < 1 {cy = 1};
    let cyd = if cy == 1 {NC1} else {cy - 1};
    let cyu = if cy == NC1 {1} else {cy + 1};
    let cxl = if cx == 1 {NC1} else {cx - 1};
    let cxr = if cx == NC1 {1} else {cx + 1};
    [NC1 * (cy - 1) + cx, NC1 * (cy - 1) + cxl, NC1 * (cy - 1) + cxr, NC1 * (cyu - 1) + cx, NC1 * (cyd - 1) + cx, 
    NC1 * (cyu - 1) + cxl, NC1 * (cyu - 1) + cxr, NC1 * (cyd - 1) + cxl, NC1 * (cyd - 1) + cxr]
}

// general functions for accelerations
fn wx(x: f64, y: f64, nano_in_cell: &[usize; NC], nano_list: &[[(f64, f64); NN]; NC]) -> f64 {
    let nc = cell_numbers(x, y);
    let mut wres = 0.0;
    for jc in 0..9 {
        let sc = nano_in_cell[nc[jc]-1];
        if sc > 0 {
            for j in 1..sc {
                let n_l = nano_list[nc[jc]-1][j-1];
                wres = wres + wx1(x, y, n_l.0, n_l.1);
            }
        }
    }
    wres
}

fn wy(x: f64, y: f64, nano_in_cell: &[usize; NC], nano_list: &[[(f64, f64); NN]; NC]) -> f64 {
    let nc = cell_numbers(x, y);
    let mut wres = 0.0;
    for jc in 0..9 {
        let sc = nano_in_cell[nc[jc]-1];
        if sc > 0 {
            for j in 1..sc {
                let n_l = nano_list[nc[jc]-1][j-1];
                wres = wres + wy1(x, y, n_l.0, n_l.1);
            }
        }
    }
    wres
}

// Equations right hand sides
fn fx(x: f64, y: f64, vx: f64, vy: f64, nano_in_cell: &[usize; NC], nano_list: &[[(f64, f64); NN]; NC]) -> f64 {
    CWX * wx(x, y, nano_in_cell, nano_list) + CM * EE + CM * EB * vy - CTAU * vx
}
fn fy(x: f64, y: f64, vx: f64, vy: f64, nano_in_cell: &[usize; NC], nano_list: &[[(f64, f64); NN]; NC]) -> f64 {
    CWY * wy(x, y, nano_in_cell, nano_list) - CM * EB * vx - CTAU * vy
}