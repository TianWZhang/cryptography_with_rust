use num_complex::Complex64;
use rand::Rng;
use std::f64::consts::PI;

use crate::math::{
    fft::{fft_special, fft_special_inv},
    finite_field::{add_mod, inv_mod, mul_mod, neg_mod, sub_mod},
    ntt::{inv_ntt_radix2_u64, ntt_radix2_u64},
    prime::{get_primitive_root_of_unity, is_prime},
};

use super::approx_modulus_switch::fast_basis_conversion;

const Q0_BIT_SIZE: u32 = 61;

pub struct Context {
    pub(crate) n: u64, // ring dimension
    pub(crate) max_level: usize,
    pub(crate) num_special_modulus: usize, // usually max_level + 1
    pub(crate) p: u64,
    p_vec: Vec<u64>, // p_0, p_1, ..., p_{k-1}, where k is the number of special modulus
    pub(crate) q_vec: Vec<u64>, // q_0, q_1, ..., q_{L-1}, where L is maximal level
    p_roots: Vec<u64>,
    q_roots: Vec<u64>,
    p_roots_inv: Vec<u64>,
    q_roots_inv: Vec<u64>,
    p_roots_pows: Vec<Vec<u64>>,
    q_roots_pows: Vec<Vec<u64>>,
    p_roots_pows_inv: Vec<Vec<u64>>,
    q_roots_pows_inv: Vec<Vec<u64>>,
    p_hat_mod_p: Vec<u64>,
    q_hat_mod_q: Vec<Vec<u64>>,
    p_hat_inv_mod_p: Vec<u64>,
    q_hat_inv_mod_q: Vec<Vec<u64>>,
    q_hat_mod_p: Vec<Vec<Vec<u64>>>,
    p_hat_mod_q: Vec<Vec<u64>>,
    p_mod_q: Vec<u64>,
    p_inv_mod_q: Vec<u64>,
    q_mod_p: Vec<Vec<u64>>,
    q_inv_mod_p: Vec<Vec<u64>>,
    q_inv_mod_q: Vec<Vec<u64>>,
    rot_group: Vec<usize>,
    ksi_pows: Vec<Complex64>,
    sigma: f64,
}

impl Context {
    pub fn new(n: u64, l: usize, k: usize, p: u64, sigma: f64) -> Self {
        let mut res = Self {
            n,
            max_level: l,
            num_special_modulus: k,
            p,
            p_vec: Vec::with_capacity(k),
            q_vec: Vec::with_capacity(l),
            p_roots: Vec::with_capacity(k),
            q_roots: Vec::with_capacity(l),
            p_roots_inv: Vec::with_capacity(k),
            q_roots_inv: Vec::with_capacity(l),
            p_roots_pows: vec![vec![0; n as usize]; k],
            q_roots_pows: vec![vec![0; n as usize]; l],
            p_roots_pows_inv: vec![vec![0; n as usize]; k],
            q_roots_pows_inv: vec![vec![0; n as usize]; l],
            p_hat_mod_p: vec![1; k],
            q_hat_mod_q: Vec::with_capacity(l),
            p_hat_inv_mod_p: vec![1; k],
            q_hat_inv_mod_q: Vec::with_capacity(l),
            q_hat_mod_p: Vec::with_capacity(l),
            p_hat_mod_q: Vec::with_capacity(k),
            p_mod_q: vec![1; l],
            p_inv_mod_q: vec![1; l],
            q_mod_p: vec![vec![1; k]; l],
            q_inv_mod_p: vec![vec![1; k]; l],
            q_inv_mod_q: vec![vec![0; l]; l],
            rot_group: Vec::with_capacity(n as usize / 2),
            ksi_pows: Vec::with_capacity(2 * (n as usize) + 1),
            sigma,
        };
        res.generate_primes();
        res.generate_primitive_m_th_roots();
        res.generate_rot_group();
        res.generate_ksi_powers();
        res
    }

    fn generate_primes(&mut self) {
        println!("N = {}", self.n);
        println!("M = {}", 2 * self.n);
        let m = 2 * self.n;
        let mut bnd = 1u64;
        let mut cnt = 1u64;

        loop {
            let q0 = 1u64.wrapping_shl(Q0_BIT_SIZE) + bnd * m + 1;
            if is_prime(q0) {
                self.q_vec.push(q0);
                break;
            }
            bnd += 1;
        }

        bnd = 1;
        while cnt < self.max_level as u64 {
            let p1 = self.p + bnd * m + 1;
            if is_prime(p1) {
                self.q_vec.push(p1);
                cnt += 1;
            }
            if cnt == self.max_level as u64 {
                bnd += 1;
                break;
            }
            let p2 = self.p - bnd * m + 1;
            if is_prime(p2) {
                self.q_vec.push(p2);
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
                self.p_vec.push(p1);
                cnt += 1;
            }
            if cnt == self.num_special_modulus as u64 {
                break;
            }
            let p2 = self.p - bnd * m + 1;
            if is_prime(p2) {
                self.p_vec.push(p2);
                cnt += 1;
            }
            bnd += 1;
        }

        println!("q_vec:");
        for q in &self.q_vec {
            println!("{}", q);
        }
        println!("p_vec:");
        for p in &self.p_vec {
            println!("{}", p);
        }

        // generate q_hat_mod_q
        // for l in 0..L, for i in 0..=l,
        // q_hat_mod_q[l][i] = \prod_{j != i} q_vec[j] mod q_vec[i]
        for l in 0..self.max_level {
            self.q_hat_mod_q.push(vec![1; l + 1]);
            for i in 0..=l {
                for j in 0..i {
                    self.q_hat_mod_q[l][i] =
                        mul_mod(self.q_hat_mod_q[l][i], self.q_vec[j], self.q_vec[i]);
                }
                for j in i + 1..=l {
                    self.q_hat_mod_q[l][i] =
                        mul_mod(self.q_hat_mod_q[l][i], self.q_vec[j], self.q_vec[i]);
                }
            }

            self.q_hat_inv_mod_q.push(vec![1; l + 1]);
            for i in 0..=l {
                self.q_hat_inv_mod_q[l][i] = inv_mod(self.q_hat_mod_q[l][i], self.q_vec[i]);
            }
        }

        for k in 0..self.num_special_modulus {
            for i in 0..k {
                self.p_hat_mod_p[k] = mul_mod(self.p_hat_mod_p[k], self.p_vec[i], self.p_vec[k]);
            }
            for i in k + 1..self.num_special_modulus {
                self.p_hat_mod_p[k] = mul_mod(self.p_hat_mod_p[k], self.p_vec[i], self.p_vec[k]);
            }

            self.p_hat_inv_mod_p[k] = inv_mod(self.p_hat_mod_p[k], self.p_vec[k]);
        }

        // generate q_hat_mod_p
        // for l in 0..L, for i in 0..=l, for k in 0..K
        // q_hat_mod_p[l][i][k] = \prod_{j != i, 0 <= j <= l} q_vec[j] mod p_vec[k]
        for l in 0..self.max_level {
            self.q_hat_mod_p
                .push(vec![vec![1; self.num_special_modulus]; l + 1]);
            for i in 0..=l {
                for k in 0..self.num_special_modulus {
                    for j in 0..i {
                        self.q_hat_mod_p[l][i][k] =
                            mul_mod(self.q_hat_mod_p[l][i][k], self.q_vec[j], self.p_vec[k]);
                    }
                    for j in i + 1..=l {
                        self.q_hat_mod_p[l][i][k] =
                            mul_mod(self.q_hat_mod_p[l][i][k], self.q_vec[j], self.p_vec[k]);
                    }
                }
            }
        }

        // generate p_hat_mod_q
        // for k in 0..K, for i in 0..L,
        // p_hat_mod_q[k][i] = \prod_{s != k} p_vec[s] mod q_vec[i]
        for k in 0..self.num_special_modulus {
            self.p_hat_mod_q.push(vec![1; self.max_level]);
            for i in 0..self.max_level {
                for s in 0..k {
                    self.p_hat_mod_q[k][i] =
                        mul_mod(self.p_hat_mod_q[k][i], self.p_vec[s], self.q_vec[i]);
                }
                for s in k + 1..self.num_special_modulus {
                    self.p_hat_mod_q[k][i] =
                        mul_mod(self.p_hat_mod_q[k][i], self.p_vec[s], self.q_vec[i]);
                }
            }
        }

        // generate PModq, the product of p_vec mod q_vec[i] and PInvModq
        for l in 0..self.max_level {
            for k in 0..self.num_special_modulus {
                self.p_mod_q[l] = mul_mod(self.p_mod_q[l], self.p_vec[k], self.q_vec[l]);
            }
            self.p_inv_mod_q[l] = inv_mod(self.p_mod_q[l], self.q_vec[l]);
        }

        // generate QModp, QModp[l][k] = \prod_{0 <= j <= l} q_vec[j] mod p_vec[k]
        for l in 0..self.max_level {
            for k in 0..self.num_special_modulus {
                for j in 0..=l {
                    self.q_mod_p[l][k] = mul_mod(self.q_mod_p[l][k], self.q_vec[j], self.p_vec[k]);
                }
                self.q_inv_mod_p[l][k] = inv_mod(self.q_mod_p[l][k], self.p_vec[k]);
            }
        }

        // generate QInvModq, q_inv_mod_q[i][j] = q_i^{-1} mod q_j
        for i in 0..self.max_level {
            for j in 0..i {
                self.q_inv_mod_q[i][j] = inv_mod(self.q_vec[i], self.q_vec[j]);
            }
            for j in (i + 1)..self.max_level {
                self.q_inv_mod_q[i][j] = inv_mod(self.q_vec[i], self.q_vec[j]);
            }
        }
    }

    fn generate_primitive_m_th_roots(&mut self) {
        let n = self.n;
        for i in 0..self.max_level {
            self.q_roots
                .push(get_primitive_root_of_unity(n * 2, self.q_vec[i]));
            self.q_roots_inv
                .push(inv_mod(self.q_roots[i], self.q_vec[i]));
            let mut power = 1;
            let mut power_inv = 1;
            for j in 0..(n as usize) {
                let jprime = ((j as u32).reverse_bits() >> (32 - (n as f64).log2() as u32)) as usize;
                self.q_roots_pows[i][jprime] = power;
                self.q_roots_pows_inv[i][jprime] = power_inv;
                if j < (n - 1) as usize {
                    power = mul_mod(power, self.q_roots[i], self.q_vec[i]);
                    power_inv = mul_mod(power_inv, self.q_roots_inv[i], self.q_vec[i]);
                }
            }
        }

        for i in 0..self.num_special_modulus {
            self.p_roots
                .push(get_primitive_root_of_unity(n * 2, self.p_vec[i]));
            self.p_roots_inv
                .push(inv_mod(self.p_roots[i], self.p_vec[i]));
            let mut power = 1;
            let mut power_inv = 1;
            for j in 0..(n as usize) {
                let jprime = ((j as u32).reverse_bits() >> (32 - (n as f64).log2() as u32)) as usize;
                self.p_roots_pows[i][jprime] = power;
                self.p_roots_pows_inv[i][jprime] = power_inv;
                if j < (n - 1) as usize {
                    power = mul_mod(power, self.p_roots[i], self.p_vec[i]);
                    power_inv = mul_mod(power_inv, self.p_roots_inv[i], self.p_vec[i]);
                }
            }
        }
    }

    fn generate_rot_group(&mut self) {
        let mut five_pow = 1;
        for _ in 0..(self.n as usize) / 2 {
            self.rot_group.push(five_pow);
            five_pow = (5 * five_pow) % (2 * self.n as usize);
        }
    }

    fn generate_ksi_powers(&mut self) {
        let m = 2 * self.n as usize;
        for i in 0..m {
            let angle = 2.0 * PI * i as f64 / m as f64;
            self.ksi_pows.push(Complex64::new(angle.cos(), angle.sin()));
        }
        self.ksi_pows.push(self.ksi_pows[0]);
    }

    pub fn encode(&self, v: &Vec<Complex64>, l: usize) -> Vec<u64> {
        let n = self.n as usize;
        let mut uvals = v.clone();
        let slots = v.len();
        let gap = n / (2 * slots);

        // encode slots complex numbers `v` into `uvals`
        fft_special_inv(&mut uvals, &self.ksi_pows, &self.rot_group);
        for i in 0..slots {
            uvals[i] *= self.p as f64;
        }

        let mut res = vec![0; l * n];
        for i in 0..l {
            let mut idx = 0;
            let mut jdx = n / 2;
            for j in 0..slots {
                let mir = uvals[j].re as i64;
                let mii = uvals[j].im as i64;
                res[i * n + idx] = if mir >= 0 {
                    mir as u64
                } else {
                    (self.q_vec[i] as i64 + mir) as u64
                };
                res[i * n + jdx] = if mii >= 0 {
                    mii as u64
                } else {
                    (self.q_vec[i] as i64 + mii) as u64
                };
                idx += gap;
                jdx += gap;
            }

            // `res[i * n..(i + 1) * n]` stores the components mod qi
            ntt_radix2_u64(
                &mut res[i * n..(i + 1) * n],
                &self.q_roots_pows[i],
                self.q_vec[i],
            );
        }
        res
    }

    /// The absolute value of Re(v) and Im(v) should be less than one.
    pub fn encode_single(&self, v: Complex64, l: usize) -> Vec<u64> {
        let n = self.n as usize;
        let vr = v.re * self.p as f64;
        let vi = v.im * self.p as f64;

        let mut res = vec![0; l * n];
        for i in 0..l {
            res[i * n] = if vr >= 0.0 {
                vr as u64
            } else {
                (self.q_vec[i] as f64 + vr) as u64
            };
            res[i * n + n / 2] = if vi >= 0.0 {
                vi as u64
            } else {
                (self.q_vec[i] as f64 + vi) as u64
            };

            ntt_radix2_u64(
                &mut res[i * n..(i + 1) * n],
                &self.q_roots_pows[i],
                self.q_vec[i],
            );
        }
        res
    }

    pub fn decode(&self, v: &Vec<u64>, slots: usize) -> Vec<Complex64> {
        let n = self.n as usize;
        let mut uvals = v.clone();
        let gap = n / (2 * slots);
        inv_ntt_radix2_u64(&mut uvals[..n], &self.q_roots_pows_inv[0], self.q_vec[0]);

        let mut res = Vec::with_capacity(slots);
        let mut idx = 0;
        let mut jdx = n / 2;
        for _ in 0..slots {
            let mir = if uvals[idx] <= self.q_vec[0] / 2 {
                uvals[idx] as f64 / self.p as f64
            } else {
                (uvals[idx] as f64 - self.q_vec[0] as f64) / self.p as f64
            };
            let mii = if uvals[jdx] <= self.q_vec[0] / 2 {
                uvals[jdx] as f64 / self.p as f64
            } else {
                (uvals[jdx] as f64 - self.q_vec[0] as f64) / self.p as f64
            };
            res.push(Complex64::new(mir, mii));
            idx += gap;
            jdx += gap;
        }
        fft_special(&mut res, &self.ksi_pows, &self.rot_group);
        res
    }

    pub fn decode_single(&self, v: &Vec<u64>) -> Complex64 {
        let n = self.n as usize;
        let mut uvals = v.clone();
        inv_ntt_radix2_u64(&mut uvals[..n], &self.q_roots_pows_inv[0], self.q_vec[0]);

        let vr = if uvals[0] <= self.q_vec[0] / 2 {
            uvals[0] as f64 / self.p as f64
        } else {
            (uvals[0] as f64 - self.q_vec[0] as f64) / self.p as f64
        };
        let vi = if uvals[n / 2] <= self.q_vec[0] / 2 {
            uvals[n / 2] as f64 / self.p as f64
        } else {
            (uvals[n / 2] as f64 - self.q_vec[0] as f64) / self.p as f64
        };
        Complex64::new(vr, vi)
    }

    /// Output a vector of length N * (l + k), each element is 0 with prob. 0.5, is 1 with prob. 0.25, is -1
    /// with prob. 0.25
    pub fn sample_zo(&self, l: usize, k: usize) -> Vec<u64> {
        let n = self.n as usize;
        let mut res = vec![0; n * (l + k)];
        let mut rng = rand::thread_rng();
        for i in 0..n {
            let zo = if rng.gen_bool(0.5) {
                0
            } else {
                if rng.gen_bool(0.5) {
                    1
                } else {
                    -1
                }
            };
            for j in 0..l {
                res[j * n + i] = if zo >= 0 {
                    zo as u64
                } else {
                    self.q_vec[j] - 1
                };
            }
            for j in 0..k {
                res[(j + l) * n + i] = if zo >= 0 {
                    zo as u64
                } else {
                    self.p_vec[j] - 1
                };
            }
        }
        res
    }

    pub fn sample_guass(&self, l: usize, k: usize) -> Vec<u64> {
        let n = self.n as usize;
        let mut res = vec![0; n * (l + k)];
        let mut rng = rand::thread_rng();
        const BIGNUM: u64 = 0xfffffff;

        for i in (0..n).step_by(2) {
            let r1 = (1.0 + rng.gen_range(0..BIGNUM) as f64) / (BIGNUM as f64 + 1.0);
            let r2 = (1.0 + rng.gen_range(0..BIGNUM) as f64) / (BIGNUM as f64 + 1.0);
            let theta = 2.0 * PI * r1;
            let rr = (-2.0 * r2.ln()).sqrt() * self.sigma;

            let g1 = (rr * theta.cos() + 0.5).floor() as i64;
            let g2 = (rr * theta.sin() + 0.5).floor() as i64;

            for j in 0..l {
                res[j * n + i] = if g1 >= 0 {
                    g1 as u64
                } else {
                    (self.q_vec[j] as i64 + g1) as u64
                };
                res[j * n + i + 1] = if g2 >= 0 {
                    g2 as u64
                } else {
                    (self.q_vec[j] as i64 + g2) as u64
                };
            }
            for j in 0..k {
                res[(j + l) * n + i] = if g1 >= 0 {
                    g1 as u64
                } else {
                    (self.p_vec[j] as i64 + g1) as u64
                };
                res[(j + l) * n + i + 1] = if g2 >= 0 {
                    g2 as u64
                } else {
                    (self.p_vec[j] as i64 + g2) as u64
                };
            }
        }
        res
    }

    pub fn sample_uniform(&self, l: usize, k: usize) -> Vec<u64> {
        let n = self.n as usize;
        let mut res = vec![0; n * (l + k)];
        let mut rng = rand::thread_rng();

        for i in 0..n {
            for j in 0..l {
                res[j * n + i] = rng.gen_range(0..self.q_vec[j]);
            }
            for j in 0..k {
                res[(j + l) * n + i] = rng.gen_range(0..self.p_vec[j]);
            }
        }
        res
    }

    pub fn sample_hwt(&self, l: usize, k: usize) -> Vec<u64> {
        let n = self.n as usize;
        let mut res = vec![0; n * (l + k)];
        let mut rng = rand::thread_rng();

        let h = 64;
        let mut idx = 0;
        while idx < h {
            let i = rng.gen_range(0..n);
            if res[i] == 0 {
                let hwt = if rng.gen_bool(0.5) { 1 } else { -1 };
                for j in 0..l {
                    res[j * n + i] = if hwt >= 0 { 1 } else { self.q_vec[j] - 1 };
                }
                for j in 0..k {
                    res[(j + l) * n + i] = if hwt >= 0 { 1 } else { self.p_vec[j] - 1 };
                }
                idx += 1;
            }
        }
        res
    }

    pub fn ntt(&self, x: &mut [u64], l: usize, k: usize) {
        let n = self.n as usize;
        for i in 0..l {
            ntt_radix2_u64(
                &mut x[i * n..(i + 1) * n],
                &self.q_roots_pows[i],
                self.q_vec[i],
            );
        }
        for i in 0..k {
            ntt_radix2_u64(
                &mut x[(i + l) * n..(i + l + 1) * n],
                &self.p_roots_pows[i],
                self.p_vec[i],
            );
        }
    }

    pub fn inv_ntt(&self, x: &mut [u64], l: usize, k: usize) {
        let n = self.n as usize;
        for i in 0..l {
            inv_ntt_radix2_u64(
                &mut x[i * n..(i + 1) * n],
                &self.q_roots_pows_inv[i],
                self.q_vec[i],
            );
        }
        for i in 0..k {
            inv_ntt_radix2_u64(
                &mut x[(i + l) * n..(i + l + 1) * n],
                &self.p_roots_pows_inv[i],
                self.p_vec[i],
            );
        }
    }

    pub fn mul(&self, a: &[u64], b: &[u64], l: usize, k: usize) -> Vec<u64> {
        let n = self.n as usize;
        let mut res = vec![0; n * (l + k)];
        for i in 0..l {
            for j in 0..n {
                res[i * n + j] = mul_mod(a[i * n + j], b[i * n + j], self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                res[(i + l) * n + j] =
                    mul_mod(a[(i + l) * n + j], b[(i + l) * n + j], self.p_vec[i]);
            }
        }
        res
    }

    pub fn mul_key(&self, a: &[u64], key: &[u64], l: usize) -> Vec<u64> {
        let n = self.n as usize;
        let k = self.num_special_modulus;
        let mut res = vec![0; n * (l + k)];
        for i in 0..l {
            for j in 0..n {
                res[i * n + j] = mul_mod(a[i * n + j], key[i * n + j], self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                res[(i + l) * n + j] = mul_mod(
                    a[(i + l) * n + j],
                    key[(i + self.max_level) * n + j],
                    self.p_vec[i],
                );
            }
        }
        res
    }

    pub fn mul_by_bigp(&self, a: &mut [u64], l: usize) {
        let ring_dim = self.n as usize;
        self.inv_ntt(a, l, 0);
        for i in 0..l {
            for n in 0..ring_dim {
                a[i * ring_dim + n] = mul_mod(a[i * ring_dim + n], self.p_mod_q[i], self.q_vec[i]);
            }
        }

        let k = a.len() / ring_dim - l;
        for i in 0..k {
            for n in 0..ring_dim {
                a[(i + l) * ring_dim + n] = 0;
            }
        }
        self.ntt(a, l, 0);
    }

    pub fn mul_inplace(&self, a: &mut [u64], b: &[u64], l: usize, k: usize) {
        let n = self.n as usize;
        for i in 0..l {
            for j in 0..n {
                a[i * n + j] = mul_mod(a[i * n + j], b[i * n + j], self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                a[(i + l) * n + j] = mul_mod(a[(i + l) * n + j], b[(i + l) * n + j], self.p_vec[i]);
            }
        }
    }

    pub fn sub(&self, a: &[u64], b: &[u64], l: usize, k: usize) -> Vec<u64> {
        let n = self.n as usize;
        let mut res = vec![0; n * (l + k)];
        for i in 0..l {
            for j in 0..n {
                res[i * n + j] = sub_mod(a[i * n + j], b[i * n + j], self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                res[(i + l) * n + j] =
                    sub_mod(a[(i + l) * n + j], b[(i + l) * n + j], self.p_vec[i]);
            }
        }
        res
    }

    pub fn add(&self, a: &[u64], b: &[u64], l: usize, k: usize) -> Vec<u64> {
        let n = self.n as usize;
        let mut res = vec![0; n * (l + k)];
        for i in 0..l {
            for j in 0..n {
                res[i * n + j] = add_mod(a[i * n + j], b[i * n + j], self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                res[(i + l) * n + j] =
                    add_mod(a[(i + l) * n + j], b[(i + l) * n + j], self.p_vec[i]);
            }
        }
        res
    }

    pub fn add_const(&self, a: &[u64], c: u64, l: usize, k: usize) -> Vec<u64> {
        let n = self.n as usize;
        let mut res = vec![0; n * (l + k)];
        for i in 0..l {
            for j in 0..n {
                res[i * n + j] = add_mod(a[i * n + j], c, self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                res[(i + l) * n + j] = add_mod(a[(i + l) * n + j], c, self.p_vec[i]);
            }
        }
        res
    }

    pub fn add_const_inplace(&self, a: &mut [u64], c: u64, l: usize, k: usize) {
        let n = self.n as usize;
        for i in 0..l {
            for j in 0..n {
                a[i * n + j] = add_mod(a[i * n + j], c, self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                a[(i + l) * n + j] = add_mod(a[(i + l) * n + j], c, self.p_vec[i]);
            }
        }
    }

    pub fn sub_const(&self, a: &[u64], c: u64, l: usize, k: usize) -> Vec<u64> {
        let n = self.n as usize;
        let mut res = vec![0; n * (l + k)];
        for i in 0..l {
            for j in 0..n {
                res[i * n + j] = sub_mod(a[i * n + j], c, self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                res[(i + l) * n + j] = sub_mod(a[(i + l) * n + j], c, self.p_vec[i]);
            }
        }
        res
    }

    pub fn sub_const_inplace(&self, a: &mut [u64], c: u64, l: usize, k: usize) {
        let n = self.n as usize;
        for i in 0..l {
            for j in 0..n {
                a[i * n + j] = sub_mod(a[i * n + j], c, self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                a[(i + l) * n + j] = sub_mod(a[(i + l) * n + j], c, self.p_vec[i]);
            }
        }
    }

    pub fn mul_const(&self, a: &[u64], c: u64, l: usize, k: usize) -> Vec<u64> {
        let n = self.n as usize;
        let mut res = vec![0; n * (l + k)];
        for i in 0..l {
            for j in 0..n {
                res[i * n + j] = mul_mod(a[i * n + j], c, self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                res[(i + l) * n + j] = mul_mod(a[(i + l) * n + j], c, self.p_vec[i]);
            }
        }
        res
    }

    pub fn mul_const_inplace(&self, a: &mut [u64], c: u64, l: usize, k: usize) {
        let n = self.n as usize;
        for i in 0..l {
            for j in 0..n {
                a[i * n + j] = mul_mod(a[i * n + j], c, self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                a[(i + l) * n + j] = mul_mod(a[(i + l) * n + j], c, self.p_vec[i]);
            }
        }
    }

    pub fn sub_inplace(&self, a: &mut [u64], b: &[u64], l: usize, k: usize) {
        let n = self.n as usize;
        for i in 0..l {
            for j in 0..n {
                a[i * n + j] = sub_mod(a[i * n + j], b[i * n + j], self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                a[(i + l) * n + j] = sub_mod(a[(i + l) * n + j], b[(i + l) * n + j], self.p_vec[i]);
            }
        }
    }

    pub fn add_inplace(&self, a: &mut [u64], b: &[u64], l: usize, k: usize) {
        let n = self.n as usize;
        for i in 0..l {
            for j in 0..n {
                a[i * n + j] = add_mod(a[i * n + j], b[i * n + j], self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                a[(i + l) * n + j] = add_mod(a[(i + l) * n + j], b[(i + l) * n + j], self.p_vec[i]);
            }
        }
    }

    pub fn negate(&self, a: &[u64], l: usize, k: usize) -> Vec<u64> {
        let n = self.n as usize;
        let mut res = vec![0; n * (l + k)];
        for i in 0..l {
            for j in 0..n {
                res[i * n + j] = neg_mod(a[i * n + j], self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                res[(i + l) * n + j] = neg_mod(a[(i + l) * n + j], self.p_vec[i]);
            }
        }
        res
    }

    pub fn negate_inplace(&self, a: &mut [u64], l: usize, k: usize) {
        let n = self.n as usize;
        for i in 0..l {
            for j in 0..n {
                a[i * n + j] = neg_mod(a[i * n + j], self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                a[(i + l) * n + j] = neg_mod(a[(i + l) * n + j], self.p_vec[i]);
            }
        }
    }

    /// b = a - b
    pub fn sub2_inplace(&self, a: &[u64], b: &mut [u64], l: usize, k: usize) {
        let n = self.n as usize;
        for i in 0..l {
            for j in 0..n {
                b[i * n + j] = sub_mod(a[i * n + j], b[i * n + j], self.q_vec[i]);
            }
        }
        for i in 0..k {
            for j in 0..n {
                b[(i + l) * n + j] = sub_mod(a[(i + l) * n + j], b[(i + l) * n + j], self.p_vec[i]);
            }
        }
    }

    pub fn conjugate(&self, a: &[u64], l: usize) -> Vec<u64> {
        let n = self.n as usize;
        let mut res = vec![0; n * l];
        for i in 0..l {
            for j in 0..n {
                res[i * n + j] = a[n - 1 - j + i * n];
            }
        }
        res
    }

    pub fn conjugate_inplace(&self, a: &mut [u64], l: usize) {
        let ring_dim = self.n as usize;
        for i in 0..l {
            for n in 0..ring_dim {
                a.swap(n + ring_dim * i, ring_dim - 1 - n + i * ring_dim);
            }
        }
    }

    /// Conv_{C_l->B}([a]_{C_l}) = [a + Q_l * e]_B, the input and output are of
    /// coefficient form.
    fn fast_basis_conversion_c2b(&self, a: &[u64], l: usize) -> Vec<u64> {
        // the length of a is l * ring_dim
        let ring_dim = self.n as usize;
        fast_basis_conversion(
            a,
            &self.q_vec[..l],
            &self.p_vec,
            ring_dim,
            &self.q_hat_inv_mod_q[l - 1],
            &self.q_hat_mod_p[l - 1],
        )
    }

    /// Conv_{B->C_l}([a]_B) = [a + P * e]_{C_l}, the input and output are of
    /// coefficient form.
    fn fast_basis_conversion_b2c(&self, a: &[u64], l: usize) -> Vec<u64> {
        // the length of a is num_special_modulus * ring_dim
        let ring_dim = self.n as usize;
        // the length of res is l * ring_dim
        fast_basis_conversion(
            a,
            &self.p_vec,
            &self.q_vec[..l],
            ring_dim,
            &self.p_hat_inv_mod_p,
            &self.p_hat_mod_q,
        )
    }

    /// The length of `a` is l * ring_dim. The input and output are of NTT form.
    pub fn mod_up(&self, a: &Vec<u64>, l: usize) -> Vec<u64> {
        let mut res = vec![];

        let mut a = a.clone();
        // the first l * ring_dim elements of res are of NTT form
        res.extend(a.iter());
        self.inv_ntt(&mut a, l, 0);

        let ring_dim = self.n as usize;
        res.extend(self.fast_basis_conversion_c2b(&a, l));
        // only need to convert the last K * ring_dim elements of res to NTT form
        self.ntt(&mut res[l * ring_dim..], 0, self.num_special_modulus);

        // the length of res is (l + K) * ring_dim
        res
    }

    pub fn mod_down(&self, a: &[u64], l: usize, dl: usize) -> Vec<u64> {
        let mut res = vec![];
        let ring_dim = self.n as usize;
        res.extend(&a[..(l - dl) * ring_dim]);
        res
    }

    pub fn approx_modulus_reduction(&self, b_tilde: &Vec<u64>, l: usize) -> Vec<u64> {
        let ring_dim = self.n as usize;
        let mut b_tilde_coeffi = b_tilde.clone();
        self.inv_ntt(&mut b_tilde_coeffi, l, self.num_special_modulus);
        let a_tilde = self.fast_basis_conversion_b2c(&b_tilde_coeffi[l * ring_dim..], l);

        let mut res = vec![];
        for i in 0..l {
            res.extend((0..ring_dim).map(|n| {
                let tmp = sub_mod(
                    b_tilde_coeffi[n + i * ring_dim],
                    a_tilde[n + i * ring_dim],
                    self.q_vec[i],
                );
                mul_mod(tmp, self.p_inv_mod_q[i], self.q_vec[i])
            }));
        }
        self.ntt(&mut res, l, 0);
        res
    }

    pub fn rescale(&self, a: &Vec<u64>, l: usize) -> Vec<u64> {
        let ring_dim = self.n as usize;
        let mut res = vec![0; (l - 1) * ring_dim];
        let mut a = a.clone();
        inv_ntt_radix2_u64(
            &mut a[(l - 1) * ring_dim..l * ring_dim],
            &self.q_roots_pows_inv[l - 1],
            self.q_vec[l - 1],
        );

        for i in 0..(l - 1) {
            for n in 0..ring_dim {
                res[i * ring_dim + n] = a[(l - 1) * ring_dim + n] % self.q_vec[i];
            }
            ntt_radix2_u64(
                &mut res[i * ring_dim..(i + 1) * ring_dim],
                &self.q_roots_pows[i],
                self.q_vec[i],
            );
            for n in 0..ring_dim {
                res[i * ring_dim + n] =
                    sub_mod(a[i * ring_dim + n], res[i * ring_dim + n], self.q_vec[i]);
                res[i * ring_dim + n] = mul_mod(
                    res[i * ring_dim + n],
                    self.q_inv_mod_q[l - 1][i],
                    self.q_vec[i],
                );
            }
        }
        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::heaan_rns::{approx_modulus_switch::crt_reconstruct, utils::{equal_up_to_epsilon, gen_random_complex_vector}};
    use num_bigint::{BigUint, RandBigInt};

    #[test]
    fn test_encode_decode_single() {
        // N: 16384
        // q0: 2305843009214414849
        // p: 36028797018963968
        let context = Context::new(1 << 14, 1, 1, 1 << 55, 3.2);
        let mut v = Complex64::new(0.47, 0.97);
        let mut encoding = context.encode_single(v, 1);
        let mut v_decoded = context.decode_single(&encoding);
        assert!(equal_up_to_epsilon(&[v], &[v_decoded], 0.000000000001));

        v = Complex64::new(0.31, -0.24);
        encoding = context.encode_single(v, 1);
        v_decoded = context.decode_single(&encoding);
        assert!(equal_up_to_epsilon(&[v], &[v_decoded], 0.000000000001));

        v = Complex64::new(-0.86, 0.12);
        encoding = context.encode_single(v, 1);
        v_decoded = context.decode_single(&encoding);
        assert!(equal_up_to_epsilon(&[v], &[v_decoded], 0.000000000001));
    }

    #[test]
    fn test_encode_decode_batch() {
        let l = 6;
        let k = l;
        let context = Context::new(1 << 15, l, k, 1 << 55, 3.2);
        let slots = 8;
        
        let v = gen_random_complex_vector(slots);
        let encoding = context.encode(&v, l);
        let v_decoded = context.decode(&encoding, slots);
        assert!(equal_up_to_epsilon(&v, &v_decoded, 0.001));
    }

    #[test]
    fn test_rns_ntt() {
        let mut rng = rand::thread_rng();
        let l = 6;
        let k = 6;
        let n = 1 << 15;
        let context = Context::new(n as u64, l, k, 1 << 55, 3.2);

        let length = n * (l + k);
        let mut x: Vec<u64> = (0..length).map(|_| rng.gen_range(0..1000)).collect();
        let x_original = x.clone();

        context.ntt(&mut x, l, k);
        context.inv_ntt(&mut x, l, k);
        assert_eq!(x, x_original);
    }

    #[test]
    fn test_mult_ntt() {
        let l = 2;
        let k = 0;
        let n = 1 << 15;
        let context = Context::new(n as u64, l, k, 1 << 55, 3.2);

        let length = n * (l + k);
        let mut x = vec![0; length];
        x[0] = 1;
        x[1] = 2;
        x[n + 1] = 1;
        x[n + n / 2] = 1;

        let mut x_square_expected = vec![0; length];
        x_square_expected[0] = 1;
        x_square_expected[1] = 4;
        x_square_expected[2] = 4;
        x_square_expected[n] = context.q_vec[1] - 1;
        x_square_expected[n + 2] = 1;
        x_square_expected[n + n / 2 + 1] = 2;

        context.ntt(&mut x, l, k);
        let mut x_square = context.mul(&x, &x, l, k);
        context.inv_ntt(&mut x_square, l, k);

        assert_eq!(x_square, x_square_expected);
    }

    #[test]
    fn mult_msg_then_decode_single() {
        let l = 2;
        let k = 1;
        let context = Context::new(1 << 14, l, k, 1 << 55, 3.2);
        let v1 = Complex64::new(0.47, 0.3);
        let msg1 = context.encode_single(v1, l);
        let v2 = Complex64::new(0.5, 0.74);
        let msg2 = context.encode_single(v2, l);

        let msg_mul = context.mul(&msg1, &msg2, l, 0);
        let msg_mul = context.rescale(&msg_mul, l);

        let v_mul_decoded = context.decode_single(&msg_mul);
        assert!(equal_up_to_epsilon(&[v1 * v2], &[v_mul_decoded], 0.001));
    }

    #[test]
    fn mult_msg_then_decode_batch() {
        let l = 2;
        let n = 1 << 15;
        let k = 0;
        let slots = 8;
        let context = Context::new(n as u64, l, k, 1 << 55, 3.2);

        let v1 = gen_random_complex_vector(slots);
        let v2 = gen_random_complex_vector(slots);
        let v_mul: Vec<_> = v1.iter().zip(v2.iter()).map(|(c1, c2)| c1 * c2).collect();

        let msg1 = context.encode(&v1, l);
        let msg2 = context.encode(&v2, l);

        let msg_mul = context.mul(&msg1, &msg2, l, 0);
        let msg_mul = context.rescale(&msg_mul, l);

        let v_mul_decoded = context.decode(&msg_mul, slots);
        assert!(equal_up_to_epsilon(&v_mul, &v_mul_decoded, 0.00001));
    }
}
