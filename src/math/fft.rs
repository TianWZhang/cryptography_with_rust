use num_complex::Complex64;

use super::ntt::bitrev;

pub fn fft_special(x: &mut [Complex64], pow_table: &[Complex64], rot_group: &[usize]) {
    let n = x.len();
    let m = rot_group.len() * 4;

    // bit reversal
    bitrev(x);

    let mut len = 2;
    while len <= n {
        for i in 0..n / len {
            for j in 0..len / 2 {
                let idx = (rot_group[j] % (4 * len)) * m / (4 * len);
                let u = x[i * len + j];
                let mut v = x[i * len + j + len / 2];
                v *= pow_table[idx];
                x[i * len + j + len / 2] = u - v;
                x[i * len + j] = u + v;
            }
        }
        len *= 2;
    }
}

pub fn fft_special_inv(x: &mut [Complex64], pow_table: &[Complex64], rot_group: &[usize]) {
    let n = x.len();
    let m = rot_group.len() * 4;

    let mut len = n;
    while len >= 1 {
        let lenq = len << 2;
        for i in 0..n / len {
            for j in 0..len / 2 {
                let idx = (lenq - (rot_group[j] % lenq)) * m / lenq;
                let u = x[i * len + j] + x[i * len + j + len / 2];
                let mut v = x[i * len + j] - x[i * len + j + len / 2];
                v *= pow_table[idx];
                x[i * len + j] = u;
                x[i * len + j + len / 2] = v;
            }
        }
        len /= 2;
    }
    bitrev(x);
    for i in 0..n {
        x[i] /= n as f64;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    use std::f64::consts::PI;

    fn generate_rot_group(n: usize) -> Vec<usize> {
        let mut five_pow = 1;
        let mut rot_group = Vec::with_capacity(n / 2);
        for _ in 0..n / 2 {
            rot_group.push(five_pow);
            five_pow = (5 * five_pow) % (2 * n);
        }
        rot_group
    }

    fn generate_ksi_powers(n: usize) -> Vec<Complex64> {
        let m = 2 * n;
        let mut ksi_pows = Vec::with_capacity(m + 1);
        for i in 0..m {
            let angle = 2.0 * PI * i as f64 / m as f64;
            ksi_pows.push(Complex64::new(angle.cos(), angle.sin()));
        }
        ksi_pows.push(ksi_pows[0]);
        ksi_pows
    }

    fn generate_random_complex_numbers(n: usize) -> Vec<Complex64> {
        let mut rng = rand::thread_rng();
        let mut res = Vec::with_capacity(n);
        for _ in 0..n {
            let re: f64 = rng.gen_range(-10.0..10.0);
            let img: f64 = rng.gen_range(-10.0..10.0);
            res.push(Complex64::new(re, img));
        }
        res
    }

    fn equal_up_to_epsilon(nums1: &[Complex64], nums2: &[Complex64], epsilon: f64) -> bool {
        for (x1, x2) in nums1.iter().zip(nums2.iter()) {
            if (x1.re - x2.re).abs() > epsilon || (x1.im - x2.im).abs() > epsilon {
                return false;
            }
        }
        true
    }

    #[test]
    fn test_fft_special() {
        let n = 1 << 15;
        let rot_group = generate_rot_group(n);
        let ksi_powers = generate_ksi_powers(n);
        let mut x = generate_random_complex_numbers(8);
        let x_original = x.clone();
        fft_special(&mut x, &ksi_powers, &rot_group);
        fft_special_inv(&mut x, &ksi_powers, &rot_group);
        assert!(equal_up_to_epsilon(&x, &x_original, 0.000000000001));
    }
}
