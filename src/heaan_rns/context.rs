use crate::math::prime::is_prime;

const Q0_BIT_SIZE: u64 = 61;

pub struct Context {
    n: u64, // ring dimension
    max_level: usize,
    num_special_modulus: usize, // usually max_level + 1
    p: u64,
    p_vec: Vec<u64>, // p_0, p_1, ..., p_{k-1}, where k is the number of special modulus
    q_vec: Vec<u64>, // q_0, q_1, ..., q_{L-1}, where L is maximal level
}

impl Context {
    pub fn new(n: u64, l: usize, k: usize, p: u64) -> Self {
        let mut res = Self {
            n,
            max_level: l,
            num_special_modulus: k,
            p,
            p_vec: Vec::with_capacity(k),
            q_vec: Vec::with_capacity(l)
        };
        res.generate_primes();
        res
    }

    fn generate_primes(&mut self) {
        let m = 2 * self.n;
        let mut bnd = 1u64;
        let mut cnt = 1u64;

        loop {
            let q0 = 1u64 << Q0_BIT_SIZE + bnd * m + 1;
            if is_prime(q0) {
                self.q_vec[0] = q0;
                break;
            }
            bnd += 1;
        }

        bnd = 1;
        while cnt < self.max_level as u64 {
            let p1 = self.p + bnd * m + 1;
            if is_prime(p1) {
                self.q_vec[cnt as usize] = p1;
                cnt += 1;
            }
            if cnt == self.max_level as u64 {
                break;
            }
            let p2 = self.p - bnd * m + 1;
            if is_prime(p2) {
                self.q_vec[cnt as usize] = p2;
                cnt += 1;
            }
            bnd += 1;
        }

        if self.p < 1024 * m * ((bnd as f64).log2().ceil() as u64) {
            panic!("ERROR: too small number of precision\nTry to use larger p or smaller depth");
        }

        // generate special primes
        cnt = 0;
        while cnt < self.num_special_modulus as u64 {
            let p1 = self.p + bnd * m + 1;
            if is_prime(p1) {
                self.p_vec[cnt as usize] = p1;
                cnt += 1;
            }
            if cnt == self.num_special_modulus as u64 {
                break;
            }
            let p2 = self.p - bnd * m + 1;
            if is_prime(p2) {
                self.p_vec[cnt as usize] = p2;
                cnt += 1;
            }
            bnd += 1;
        }
    }
}
