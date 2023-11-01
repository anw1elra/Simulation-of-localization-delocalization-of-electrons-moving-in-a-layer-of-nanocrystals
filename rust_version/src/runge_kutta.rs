use std::vec::Vec;

pub struct MyVec {
    vec: Vec<f64>,
    n: usize,
}

impl MyVec {
    pub fn _new_emp(n: usize) -> MyVec {
        MyVec {
            vec: vec![0.0; n],
            n,
        }
    }

    pub fn new_from(a: &Vec<f64>) -> MyVec {
        MyVec {
            vec: a.to_vec(),
            n: a.len(),
        }
    }

    pub fn add(&self, a: &MyVec) -> MyVec {
        if self.n != a.n {
            panic!("Vector addition panic: lengths not equal")
        }
        let mut c = vec![0.0; self.n];
        for j in 0..self.n {
            c[j] = self.vec[j] + a.vec[j]
        }
        MyVec::new_from(&c)
    }

    pub fn scale(&self, s: f64) -> MyVec {
        let mut c = vec![0.0; self.n];
        for j in 0..self.n {
            c[j] = self.vec[j] * s
        }
        MyVec::new_from(&c)
    }
}

pub fn runge_kutta(
    vars: &MyVec,
    pars: &Vec<f64>,
    rhs: &dyn Fn(&MyVec, &Vec<f64>) -> MyVec,
    dt: f64,
) -> MyVec {
    let rk_1 = rhs(vars, pars);
    let rk_2 = rhs(&vars.add(&rk_1.scale(dt / 2.0)), pars);
    let rk_3 = rhs(&vars.add(&rk_2.scale(dt / 2.0)), pars);
    let rk_4 = rhs(&vars.add(&rk_3.scale(dt)), pars);

    let vars_new = vars
        .add(&rk_1.scale(dt / 6.0))
        .add(&rk_2.scale(dt / 3.0))
        .add(&rk_3.scale(dt / 3.0))
        .add(&rk_4.scale(dt / 6.0));
    vars_new
}