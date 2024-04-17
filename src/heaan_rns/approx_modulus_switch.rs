use num_bigint::BigUint;

use crate::math::finite_field::{inv_mod, mul_mod};

/// C = {q_0, \ldots, q_{l-1}}, B = {p_0, \ldots, p_{k-1}}
/// Conv_{C->B}([a]_C) = [a + Q * e]_B, the input and output are of
/// coefficient form.
pub(crate) fn fast_basis_conversion(
    a: &[u64],
    basis_c: &[u64],
    basis_b: &[u64],
    ring_dim: usize,
    c_hat_inv_mod_c: &Vec<u64>,
    c_hat_mod_b: &Vec<Vec<u64>>,
) -> Vec<u64> {
    // the length of a is basis_c.len() * ring_dim
    let mut res = vec![];

    let tmp: Vec<_> = (0..basis_c.len())
        .map(|i| {
            (0..ring_dim)
                .map(move |n| mul_mod(a[n + i * ring_dim], c_hat_inv_mod_c[i], basis_c[i]) as u128)
        })
        .flatten()
        .collect();

    for k in 0..basis_b.len() {
        res.extend((0..ring_dim).map(|n| {
            let mut sum: u128 = 0;
            for i in 0..basis_c.len() {
                sum += tmp[n + i * ring_dim] * (c_hat_mod_b[i][k] as u128);
            }
            (sum % (basis_b[k] as u128)) as u64
        }));
    }

    // the length of res is basis_b.len() * ring_dim
    res
}

#[allow(dead_code)]
pub(crate) fn crt_reconstruct(remainders: &[u64], basis: &Vec<u64>) -> BigUint {
    let prod: BigUint = basis.iter().product();
    let mut res = BigUint::from(0u64);
    for (a, p) in remainders.iter().zip(basis.iter()) {
        let p_hat = &prod / p;
        let p_hat_mod_p = (&p_hat % p).try_into().unwrap();
        let p_hat_inv_mod_p = inv_mod(p_hat_mod_p, *p);
        res += BigUint::from(mul_mod(*a, p_hat_inv_mod_p, *p)) * p_hat;
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fast_basis_conversion() {
        let ring_dim = 1;
        let basis_b = vec![23, 29];
        let basis_c = vec![31, 41];
        let c_hat_inv_mod_c = vec![28, 4];
        let mut c_hat_mod_b = vec![];
        c_hat_mod_b.push(vec![18, 12]);
        c_hat_mod_b.push(vec![8, 2]);
        let a_num = 672;
        let a = vec![a_num % basis_c[0], a_num % basis_c[1]];
        let a_num = crt_reconstruct(&a, &basis_c);

        let b1 = fast_basis_conversion(
            &a,
            &basis_c,
            &basis_b,
            ring_dim,
            &c_hat_inv_mod_c,
            &c_hat_mod_b,
        );
        let b2 = vec![
            (&a_num % basis_b[0]).try_into().unwrap(),
            (&a_num % basis_b[1]).try_into().unwrap(),
        ];
        assert_eq!(b1, b2);
    }
}
