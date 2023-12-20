use std::ops::{Mul, Sub};
use super::finite_field::FiniteField;

fn is_power_of_two(n: usize) -> bool {
    1 << (63 - n.leading_zeros()) == n
}

/// bit reversal
/// the length of x should be a power of two
fn bitrev<T: Copy>(x: &mut [T]) {
    let n = x.len();
    if !is_power_of_two(n) {
        panic!("The length n of x must be a power of two");
    }

    let mut rho = vec![0usize; n];
    let mut k = 2;

    while k <= n {
        // compute rho_k(0: k-1)
        for i in 0..k / 2 {
            rho[i + k / 2] = 2 * rho[i] + 1;
            rho[i] = 2 * rho[i];
        }
        k *= 2;
    }

    for i in 0..n {
        if i < rho[i] {
            let tem = x[i];
            x[i] = x[rho[i]];
            x[rho[i]] = tem;
        }
    }
}

/// Computes the forward number-theoretic transform of the given vector x in place,
/// with respect to the given primitive nth root of unity under the given modulus.
/// The length of the x must be a power of 2.
pub fn ntt_radix2<T>(x: &mut [T], pow_table: &[T])
where
    T: FiniteField + Copy + Sub<Output = T> + Mul<Output = T>,
{
    let n = x.len();
    if !is_power_of_two(n) {
        panic!("Length is not a power of 2");
    }

    // bit reversal
    bitrev(x);

    let mut k = 2;
    while k <= n {
        for r in 0..n / k {
            for j in 0..k / 2 {
                let tau = pow_table[n / k * j] * x[r * k + j + k / 2];
                x[r * k + j + k / 2] = x[r * k + j] - tau;
                x[r * k + j] += tau;
            }
        }
        k *= 2;
    }
}

// Returns the inverse number-theoretic transform of the given vector x.
pub fn inv_ntt_radix2<T>(x: &mut [T], inv_pow_table: &[T])
where
    T: FiniteField + Copy + Sub<Output = T> + Mul<Output = T> + From<u64>,
{
    let n = x.len();
    ntt_radix2(x, inv_pow_table);
    let n_inv = T::from(n as u64).inv().unwrap();
    for i in 0..n {
        x[i] = x[i] * n_inv;
    }
}

pub fn generate_power_table<T>(root: &T, n: usize) -> Vec<T>
where
    T: FiniteField + Copy + Sub<Output = T> + Mul<Output = T> + From<u64>,
{
    let mut pow_table = vec![T::ZERO; n / 2];
    let mut temp = T::ONE;
    for i in 0..n / 2 {
        pow_table[i] = temp;
        temp = temp * (*root);
    }
    pow_table
}

#[cfg(test)]
mod tests {
    use crate::math::finite_field::Fp;

    use super::*;

    #[test]
    fn test_is_power_of_two() {
        let mut n = 8;
        assert!(is_power_of_two(n));

        n = 256;
        assert!(is_power_of_two(n));

        n = 78;
        assert!(!is_power_of_two(n));
    }

    #[test]
    fn test_bitrev() {
        let mut x = [0, 1, 2, 3, 4, 5, 6, 7];
        bitrev(&mut x);
        assert_eq!(x, [0, 4, 2, 6, 1, 5, 3, 7]);
    }

    #[test]
    fn test_ntt() {
        let mut x: Vec<_> = [6, 0, 10, 7, 2, 8, 7, 4]
            .iter()
            .map(|val| Fp::from(*val))
            .collect();
        let original_x = x.clone();
        let pow_table = generate_power_table(&Fp::<17>::get_generator(), x.len());
        let inv_pow_table = generate_power_table(&Fp::<17>::get_generator().inv().unwrap(), x.len());
        ntt_radix2(&mut x, &pow_table);
        inv_ntt_radix2(&mut x, &inv_pow_table);
        assert_eq!(x, original_x);

        let mut x: Vec<_> = [8, 2, 104, 57, 42, 18, 37, 46, 33, 5, 62, 15, 31, 88, 3, 108]
            .iter()
            .map(|val| Fp::from(*val))
            .collect();
        let original_x = x.clone();
        let pow_table = generate_power_table(&Fp::<113>::get_generator(), x.len());
        let inv_pow_table = generate_power_table(&Fp::<113>::get_generator().inv().unwrap(), x.len());
        ntt_radix2(&mut x, &pow_table);
        inv_ntt_radix2(&mut x, &inv_pow_table);
        assert_eq!(x, original_x);
    }
}
