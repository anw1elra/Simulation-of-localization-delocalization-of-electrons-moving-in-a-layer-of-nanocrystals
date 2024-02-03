use rand::Rng;
use std::fs;
use std::io::Write;
use std::io;
use std::path::PathBuf;
use std::time::Instant;
use std::env;

mod draw;
mod palettes;
mod plot;
mod file;

//Reciprocal electron mass in nm^2/fs^2/eV
const RME: f64 = 0.176;

// // Input parameters
const ME: f64 = 0.067; // electron effective mass
const TAU: f64 = 4000.0; // scattering time, fs
const U0: f64 = -0.2; // depth of QD potential well, eV
const A: f64 = 50.0; // quantum dot size x, nm
const B: f64 = 40.0; // quantum dot size y, nm
const L: f64 = 2000.0; // size of calculation area, nm
const NP: usize = 2000; // number of grid points for plotting U(x,y) map
const WIDTH: usize = 2000; // U(x,y) png map resolution
const ELECTRIC_MAX: f64 = 5.0 * 1e-4; // maximum electric field, eV/nm
const VMAX: f64 = 5.0;
const IV_POINTS: usize = 50; // number of points on IV curve
const MAGNETIC: f64 = 0.0; // magnetic field induction,
const DT: f64 = 1.0; // delta_time, fs
const NE: usize = 100; // number of electrons
const NEST: usize = 10; // number of electron trajectories to save
const NC1: usize = 10; // number of cells (to sort electtons and QD)

// Auxillary parameters
const NC: usize = NC1 * NC1;
const TF: f64 = 5.0 * TAU;
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

    // input number of quantum dots
    let mut input_line = String::new();
    io::stdin().read_line(&mut input_line).expect("Failed to read line");
    let nqd: usize = input_line.trim().parse().expect("Input not an integer");

    // Path::new("/etc/hosts").exists()
    let cur_dir = env::current_dir().unwrap();
    let ancestor_dir = &cur_dir.ancestors().nth(1).unwrap().join("data");

    if !ancestor_dir.exists() {
        fs::create_dir_all(ancestor_dir).expect("Error creating directory");
    }

    // create dir
    let dir_name = &(ancestor_dir.join(&format!("NQD={}/results", nqd)).into_os_string().into_string().unwrap());
    println!("{}", dir_name);
    fs::create_dir_all(dir_name).expect("Error creating directory");

    // random number generator
    let mut rng = rand::thread_rng();

    // coordinates of nanostructures, cell numbers of nanostructures
    let mut xn = vec![0.0; nqd]; // : [f64; NQD] = [0.0; NQD];
    let mut yn = vec![0.0; nqd];//: [f64; NQD] = [0.0; NQD];
    let mut cn = vec![1; nqd];//: [usize; NQD] = [1; NQD];

    // initiate random coordinates, but not too close to the boundary
    for jn in 1..(nqd + 1) {
        let rng_x: f64 = rng.gen();
        let rng_y: f64 = rng.gen();
        xn[jn - 1] = W + rng_x * (L - 2.0 * W);
        yn[jn - 1] = W + rng_y * (L - 2.0 * W);
    }

    //file for nanostructures data
    let fl_name = "nanostructures.dat";
    let file_path: PathBuf = [dir_name, fl_name].iter().collect();
    let mut my_file_1 = fs::File::create(file_path).expect("Error creating file");

    // define cell number for each nanostructure, write data to file
    for jn in 1..(nqd + 1) {
        cn[jn - 1] = cell_numbers(xn[jn - 1], yn[jn - 1])[0];
        writeln!(my_file_1, "{} {} {}", xn[jn - 1], yn[jn - 1], cn[jn - 1])
            .expect("Error writing to file");
    }

    //file for cells data
    let fl_name = "cells.dat";
    let file_path: PathBuf = [dir_name, fl_name].iter().collect();
    let mut my_file_2 = fs::File::create(file_path).expect("Error creating file");

    let mut nano_in_cell: [usize; NC] = [0; NC];
    let mut nano_list = vec![vec![(0.0, 0.0); nqd]; NC];

    for jn in 1..(nqd + 1) {
        let jc = cn[jn - 1];
        let nc = nano_in_cell[jc - 1];
        nano_in_cell[jc - 1] = nc + 1;
        nano_list[jc - 1][nc].0 = xn[jn - 1];
        nano_list[jc - 1][nc].1 = yn[jn - 1];
    }

    for jc in 1..(NC + 1) {
        writeln!(my_file_2, "{}", nano_in_cell[jc - 1]).expect("Error writing to file");
    }

    //PLOT POTENTIAL MAP U(x,y)
    let mut uxy_map: GridValues = GridValues::new_square(L, NP);
    uxy_map.apply_func_u(&u,&nano_in_cell,&nano_list);
    draw::plot_2d(uxy_map, WIDTH, WIDTH, "potential_map.png");

    //panic!("Just needed to plot the potential");

    //file for IV curve data
    let fl_name = "IV_curve.dat";
    let file_path: PathBuf = [dir_name,fl_name].iter().collect();
    let mut my_file_3 = fs::File::create(file_path).expect("Error creating file");

    // Start calculations
    // Electric field loop
    for j_field in 1..=IV_POINTS {
        let ee = j_field as f64 * ELECTRIC_MAX / IV_POINTS as f64;
        let mut v_array = [0.0; NE];
        let mut v_av_time = 0.0;
        let mut en_array = [0.0; NE];
        let mut en_av_time = 0.0;

        //file for IV curve data
        let fl_name = format!("Velocity_distribution_E-{}.dat", j_field);
        let file_path: PathBuf = [dir_name, &fl_name].iter().collect();
        let mut my_file_4 = fs::File::create(file_path).expect("Error creating file");

        // PARTICLE LOOP
        for je in 0..NE {
            //file for trajectory data for a single electron

            let mut fl_name = format!("Empty file.dat");
            // ONLY RECORD 5 TRAJECTORIES
            if je < NEST {
                fl_name = format!("Trajectory_e{}_E-{}.dat", je + 1, j_field)
            }
            let file_path: PathBuf = [dir_name, &fl_name].iter().collect();
            let mut my_file = fs::File::create(file_path).expect("Error creating file");

            // initial data
            let rnx: f64 = rng.gen();
            let rny: f64 = rng.gen();
            let mut x0: f64 = W + (L - 2.0 * W) * rnx;
            let mut y0: f64 = W + (L - 2.0 * W) * rny;
            let mut vx0: f64 = 0.0;
            let mut vy0: f64 = 0.0;
            // total energy (kinetic + potential)
            let mut energy: f64 =
                vx0 * vx0 + vy0 * vy0 + 2.0 * CM * U0 * u(x0, y0, &nano_in_cell, &nano_list)
                    - 2.0 * CM * ee * x0;

            // TIME LOOP
            for _jt in 0..NT {
                if je < 10 {
                    writeln!(my_file, "{} {}", x0, y0)
                        .expect("Error writing to file");
                }
                // calculating average velocity for a single electron over time
                v_av_time += vx0;
                en_av_time += energy;

                // Runge-Kutta method start
                let rk_vx1 = fx(x0, y0, vx0, vy0, ee, &nano_in_cell, &nano_list);
                let rk_vy1 = fy(x0, y0, vx0, vy0, &nano_in_cell, &nano_list);
                let rk_x1 = vx0;
                let rk_y1 = vy0;

                let x1 = x0 + rk_x1 * DT * 0.5;
                let y1 = y0 + rk_y1 * DT * 0.5;
                let vx1 = vx0 + rk_vx1 * DT * 0.5;
                let vy1 = vy0 + rk_vy1 * DT * 0.5;

                let rk_vx2 = fx(x1, y1, vx1, vy1, ee, &nano_in_cell, &nano_list);
                let rk_vy2 = fy(x1, y1, vx1, vy1, &nano_in_cell, &nano_list);
                let rk_x2 = vx1;
                let rk_y2 = vy1;

                let x1 = x0 + rk_x2 * DT * 0.5;
                let y1 = y0 + rk_y2 * DT * 0.5;
                let vx1 = vx0 + rk_vx2 * DT * 0.5;
                let vy1 = vy0 + rk_vy2 * DT * 0.5;

                let rk_vx3 = fx(x1, y1, vx1, vy1, ee, &nano_in_cell, &nano_list);
                let rk_vy3 = fy(x1, y1, vx1, vy1, &nano_in_cell, &nano_list);
                let rk_x3 = vx1;
                let rk_y3 = vy1;

                let x1 = x0 + rk_x3 * DT;
                let y1 = y0 + rk_y3 * DT;
                let vx1 = vx0 + rk_vx3 * DT;
                let vy1 = vy0 + rk_vy3 * DT;

                let rk_vx4 = fx(x1, y1, vx1, vy1, ee, &nano_in_cell, &nano_list);
                let rk_vy4 = fy(x1, y1, vx1, vy1, &nano_in_cell, &nano_list);
                let rk_x4 = vx1;
                let rk_y4 = vy1;

                x0 = x0 + (rk_x1 + 2.0 * rk_x2 + 2.0 * rk_x3 + rk_x4) * DT / 6.0;
                y0 = y0 + (rk_y1 + 2.0 * rk_y2 + 2.0 * rk_y3 + rk_y4) * DT / 6.0;
                vx0 = vx0 + (rk_vx1 + 2.0 * rk_vx2 + 2.0 * rk_vx3 + rk_vx4) * DT / 6.0;
                vy0 = vy0 + (rk_vy1 + 2.0 * rk_vy2 + 2.0 * rk_vy3 + rk_vy4) * DT / 6.0;

                energy =
                    vx0 * vx0 + vy0 * vy0 + 2.0 * CM * U0 * u(x0, y0, &nano_in_cell, &nano_list)
                        - 2.0 * CM * ee * x0;

                // Periodic boundary conditions
                if x0 > L {
                    x0 = x0 - L
                };
                if x0 < 0.0 {
                    x0 = x0 + L
                };
                if y0 > L {
                    y0 = y0 - L
                };
                if y0 < 0.0 {
                    y0 = y0 + L
                };
            }

            // dividing by number of time steps
            v_av_time = v_av_time / NT as f64;
            en_av_time = en_av_time / NT as f64;
            // recording the average velocity in an array
            v_array[je] = v_av_time;
            en_array[je] = en_av_time;
        }

        // write data for velocity and energy distributions
        for je in 0..NE {
            writeln!(my_file_4, "{} {}", v_array[je], en_array[je]).expect("Error writing to file");
        }

        // write IV data to file, I is velocity, V is ee
        writeln!(my_file_3, "{} {} {}", ee * 1e4, aver(v_array), aver(en_array))
            .expect("Error writing to file");
    }

    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);

    plot::plot_data(nqd, ancestor_dir).ok();
}

// function to define a singe nanostructure
fn u1(x: f64, y: f64, xn: f64, yn: f64) -> f64 {
    if (x - xn).abs() < (5.0 * A) && (y - yn).abs() < (5.0 * B) {
        (-CA * (x - xn) * (x - xn) - CB * (y - yn) * (y - yn)).exp()
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
    if cx > NC1 {
        cx = NC1
    };
    if cy > NC1 {
        cy = NC1
    };
    if cx < 1 {
        cx = 1
    };
    if cy < 1 {
        cy = 1
    };
    let cyd = if cy == 1 { NC1 } else { cy - 1 };
    let cyu = if cy == NC1 { 1 } else { cy + 1 };
    let cxl = if cx == 1 { NC1 } else { cx - 1 };
    let cxr = if cx == NC1 { 1 } else { cx + 1 };
    [
        NC1 * (cy - 1) + cx,
        NC1 * (cy - 1) + cxl,
        NC1 * (cy - 1) + cxr,
        NC1 * (cyu - 1) + cx,
        NC1 * (cyd - 1) + cx,
        NC1 * (cyu - 1) + cxl,
        NC1 * (cyu - 1) + cxr,
        NC1 * (cyd - 1) + cxl,
        NC1 * (cyd - 1) + cxr,
    ]
}

// general function for potential energy
fn u(x: f64, y: f64, nano_in_cell: &[usize; NC], nano_list: &[Vec<(f64, f64)>]) -> f64 {
    let nc = cell_numbers(x, y);
    let mut ures = 0.0;
    for jc in 0..9 {
        let sc = nano_in_cell[nc[jc] - 1];
        if sc > 0 {
            for j in 1..sc {
                let n_l = nano_list[nc[jc] - 1][j - 1];
                ures = ures + u1(x, y, n_l.0, n_l.1);
            }
        }
    }
    ures
}

// general functions for accelerations
fn wx(x: f64, y: f64, nano_in_cell: &[usize; NC], nano_list: &[Vec<(f64, f64)>]) -> f64 {
    let nc = cell_numbers(x, y);
    let mut wres = 0.0;
    for jc in 0..9 {
        let sc = nano_in_cell[nc[jc] - 1];
        if sc > 0 {
            for j in 1..sc {
                let n_l = nano_list[nc[jc] - 1][j - 1];
                wres = wres + wx1(x, y, n_l.0, n_l.1);
            }
        }
    }
    wres
}

fn wy(x: f64, y: f64, nano_in_cell: &[usize; NC], nano_list: &[Vec<(f64, f64)>]) -> f64 {
    let nc = cell_numbers(x, y);
    let mut wres = 0.0;
    for jc in 0..9 {
        let sc = nano_in_cell[nc[jc] - 1];
        if sc > 0 {
            for j in 1..sc {
                let n_l = nano_list[nc[jc] - 1][j - 1];
                wres = wres + wy1(x, y, n_l.0, n_l.1);
            }
        }
    }
    wres
}

// Equations right hand sides
fn fx(
    x: f64,
    y: f64,
    vx: f64,
    vy: f64,
    ee: f64,
    nano_in_cell: &[usize; NC],
    nano_list: &[Vec<(f64, f64)>],
) -> f64 {
    CWX * wx(x, y, nano_in_cell, nano_list) + CM * ee + CM * MAGNETIC * vy - CTAU * vx
}
fn fy(
    x: f64,
    y: f64,
    vx: f64,
    vy: f64,
    nano_in_cell: &[usize; NC],
    nano_list: &[Vec<(f64, f64)>],
) -> f64 {
    CWY * wy(x, y, nano_in_cell, nano_list) - CM * MAGNETIC * vx - CTAU * vy
}

// Average velocity find
fn aver(v_array: [f64; NE]) -> f64 {
    let mut v_av = 0.0;
    for je in 0..NE {
        v_av += v_array[je]
    }
    v_av / NE as f64
}

pub struct GridValues {
    uxy: Vec<f64>,
    umax: f64,
    umin: f64,
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
}

impl GridValues {
    fn new_square(l: f64, n: usize) -> GridValues {
        GridValues {
            uxy: Vec::new(),
            umax: 0.0,
            umin: 0.0,
            nx: n,
            ny: n,
            dx: l / n as f64,
            dy: l / n as f64,
        }
    }

    fn apply_func_u(&mut self,
        u: &dyn Fn(f64, f64, &[usize; NC], &[Vec<(f64, f64)>]) -> f64,    
        nano_in_cell: &[usize; NC], 
        nano_list: &[Vec<(f64, f64)>],
    ) {
        for jx in 0..self.nx {
            for jy in 0..self.ny {
                let x = self.dx * jx as f64;
                let y = self.dy * jy as f64;
                self.uxy.push(u(x, y, nano_in_cell, nano_list))
            }
        }

        self.umax = self.uxy.iter().cloned().fold(0./0., f64::max);
        self.umin = self.uxy.iter().cloned().fold(0./0., f64::min);
    }

    fn _apply_func(&mut self,
        f: &dyn Fn(f64, f64) -> f64,    
    ) {
        for jx in 0..self.nx {
            for jy in 0..self.ny {
                let x = self.dx * jx as f64;
                let y = self.dy * jy as f64;
                self.uxy.push(f(x, y))
            }
        }

        self.umax = self.uxy.iter().cloned().fold(0./0., f64::max);
        self.umin = self.uxy.iter().cloned().fold(0./0., f64::min);
    }
}

fn _plot_chek(x:  f64, y: f64) -> f64 {
    (-2.0*(x-0.5*L)*(x-0.5*L)/L/L-2.0*(y-0.5*L)*(y-0.5*L)/L/L).exp()
}