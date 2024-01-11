use super::{
    finite_field::{add_mod, div_mod, mul_mod, sub_mod, FiniteField, FiniteRing},
    Poly3329, PolyMatrix3329, PolyVec3329, F3329,
};

fn is_power_of_two(n: usize) -> bool {
    1 << (63 - n.leading_zeros()) == n
}

/// bit reversal
/// the length of x should be a power of two
pub(crate) fn bitrev<T: Copy>(x: &mut [T]) {
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
/// This is used for polynomials over Z_q[X]/(X^n - 1)
pub fn ntt_radix2<T: FiniteField>(x: &mut [T], pow_table: &[T]) {
    let n = x.len();
    if !is_power_of_two(n) {
        panic!("Length is not a power of 2");
    }

    // bit reversal
    bitrev(x);

    let mut len = 2;
    while len <= n {
        // compute x = (I_{n/len}\bigtimes B_{len})x
        for r in 0..n / len {
            // compute x[r * len..(r + 1) * len] = B_k * x[r * len..(r + 1) * len]
            for j in 0..len / 2 {
                // w_{len}^j, w_{len} = w_n^{n/len} is the len-th primitive root of unity
                let tau = pow_table[n / len * j] * x[r * len + j + len / 2];
                x[r * len + j + len / 2] = x[r * len + j] - tau;
                x[r * len + j] += tau;
            }
        }
        len *= 2;
    }
}

// Returns the inverse number-theoretic transform of the given vector x.
pub fn inv_ntt_radix2<T>(x: &mut [T], inv_pow_table: &[T])
where
    T: FiniteField + From<u64>,
{
    let n = x.len();
    ntt_radix2(x, inv_pow_table);
    let n_inv = T::from(n as u64).inv().unwrap();
    for i in 0..n {
        x[i] = x[i] * n_inv;
    }
}

// To ensure correctness, the maximal value of x should be less than p.
pub fn ntt_radix2_u64(x: &mut [u64], pow_table: &[u64], p: u64) {
    let n = x.len();
    if !is_power_of_two(n) {
        panic!("Length is not a power of 2");
    }

    bitrev(x);

    let mut len = 2;
    while len <= n {
        for r in 0..n / len {
            for j in 0..len / 2 {
                let tau = mul_mod(pow_table[n / len * j], x[r * len + j + len / 2], p);
                x[r * len + j + len / 2] = sub_mod(x[r * len + j], tau, p);
                x[r * len + j] = add_mod(x[r * len + j], tau, p);
            }
        }
        len *= 2;
    }
}

pub fn inv_ntt_radix2_u64(x: &mut [u64], inv_pow_table: &[u64], p: u64) {
    let n = x.len();
    ntt_radix2_u64(x, inv_pow_table, p);
    for i in 0..n {
        x[i] = div_mod(x[i], n as u64, p);
    }
}

pub fn generate_power_table<T>(root: &T, n: usize) -> Vec<T>
where
    T: FiniteField + From<u64>,
{
    let mut pow_table = vec![T::ZERO; n / 2];
    let mut temp = T::ONE;
    for i in 0..n / 2 {
        pow_table[i] = temp;
        temp = temp * (*root);
    }
    pow_table
}

/// 17^128 = -1 mod KYBER_Q, 17 = 256-th root of unity mod KYBER_Q
/// ZETAS[i] = 17^{brv_{7}(i)}
const ZETAS: [u64; 128] = [
    1, 1729, 2580, 3289, 2642, 630, 1897, 848, 1062, 1919, 193, 797, 2786, 3260, 569, 1746, 296,
    2447, 1339, 1476, 3046, 56, 2240, 1333, 1426, 2094, 535, 2882, 2393, 2879, 1974, 821, 289, 331,
    3253, 1756, 1197, 2304, 2277, 2055, 650, 1977, 2513, 632, 2865, 33, 1320, 1915, 2319, 1435,
    807, 452, 1438, 2868, 1534, 2402, 2647, 2617, 1481, 648, 2474, 3110, 1227, 910, 17, 2761, 583,
    2649, 1637, 723, 2288, 1100, 1409, 2662, 3281, 233, 756, 2156, 3015, 3050, 1703, 1651, 2789,
    1789, 1847, 952, 1461, 2687, 939, 2308, 2437, 2388, 733, 2337, 268, 641, 1584, 2298, 2037,
    3220, 375, 2549, 2090, 1645, 1063, 319, 2773, 757, 2099, 561, 2466, 2594, 2804, 1092, 403,
    1026, 1143, 2150, 2775, 886, 1722, 1212, 1874, 1029, 2110, 2935, 885, 2154,
];

/// For i in [0, 64), ZETAS_INV[i] = ZETAS[64 + i]^{-1} mod 3329          1175-3312
/// For i in [64, 96), ZETAS_INV[i] = ZETAS[32 + i - 64]^{-1} mod 3329    2419-3040
/// For i in [96, 112), ZETAS_INV[i] = ZETAS[16 + i - 96]^{-1} mod 3329
/// For i in [112, 120), ZETAS_INV[i] = ZETAS[8 + i - 112]^{-1} mod 3329
/// For i in [120, 124), ZETAS_INV[i] = ZETAS[4 + i - 120]^{-1} mod 3329
/// For i in [124, 126), ZETAS_INV[i] = ZETAS[2 + i - 124]^{-1} mod 3329
/// ZETAS_INV[127] = 128^{-1} mod 3329
const ZETAS_INV: [u64; 128] = [
    1175, 2444, 394, 1219, 2300, 1455, 2117, 1607, 2443, 554, 1179, 2186, 2303, 2926, 2237, 525,
    735, 863, 2768, 1230, 2572, 556, 3010, 2266, 1684, 1239, 780, 2954, 109, 1292, 1031, 1745,
    2688, 3061, 992, 2596, 941, 892, 1021, 2390, 642, 1868, 2377, 1482, 1540, 540, 1678, 1626, 279,
    314, 1173, 2573, 3096, 48, 667, 1920, 2229, 1041, 2606, 1692, 680, 2746, 568, 3312, 2419, 2102,
    219, 855, 2681, 1848, 712, 682, 927, 1795, 461, 1891, 2877, 2522, 1894, 1010, 1414, 2009, 3296,
    464, 2697, 816, 1352, 2679, 1274, 1052, 1025, 2132, 1573, 76, 2998, 3040, 2508, 1355, 450, 936,
    447, 2794, 1235, 1903, 1996, 1089, 3273, 283, 1853, 1990, 882, 3033, 1583, 2760, 69, 543, 2532,
    3136, 1410, 2267, 2481, 1432, 2699, 687, 40, 749, 1600, 3303,
];

/// PRIMITIVE_ZETAS[i] = KYBER_ROOT_OF_UNITY^{2 * brv_7(i) + 1}
const PRIMITIVE_ZETAS: [u64; 128] = [
    17, 3312, 2761, 568, 583, 2746, 2649, 680, 1637, 1692, 723, 2606, 2288, 1041, 1100, 2229, 1409,
    1920, 2662, 667, 3281, 48, 233, 3096, 756, 2573, 2156, 1173, 3015, 314, 3050, 279, 1703, 1626,
    1651, 1678, 2789, 540, 1789, 1540, 1847, 1482, 952, 2377, 1461, 1868, 2687, 642, 939, 2390,
    2308, 1021, 2437, 892, 2388, 941, 733, 2596, 2337, 992, 268, 3061, 641, 2688, 1584, 1745, 2298,
    1031, 2037, 1292, 3220, 109, 375, 2954, 2549, 780, 2090, 1239, 1645, 1684, 1063, 2266, 319,
    3010, 2773, 556, 757, 2572, 2099, 1230, 561, 2768, 2466, 863, 2594, 735, 2804, 525, 1092, 2237,
    403, 2926, 1026, 2303, 1143, 2186, 2150, 1179, 2775, 554, 886, 2443, 1722, 1607, 1212, 2117,
    1874, 1455, 1029, 2300, 2110, 1219, 2935, 394, 885, 2444, 2154, 1175,
];

/*
// Code to generate zetas and zetas_inv used in the NTT:
const TREE: [usize; 128] = [
    0, 64, 32, 96, 16, 80, 48, 112, 8, 72, 40, 104, 24, 88, 56, 120, 4, 68, 36, 100, 20, 84, 52,
    116, 12, 76, 44, 108, 28, 92, 60, 124, 2, 66, 34, 98, 18, 82, 50, 114, 10, 74, 42, 106, 26, 90,
    58, 122, 6, 70, 38, 102, 22, 86, 54, 118, 14, 78, 46, 110, 30, 94, 62, 126, 1, 65, 33, 97, 17,
    81, 49, 113, 9, 73, 41, 105, 25, 89, 57, 121, 5, 69, 37, 101, 21, 85, 53, 117, 13, 77, 45, 109,
    29, 93, 61, 125, 3, 67, 35, 99, 19, 83, 51, 115, 11, 75, 43, 107, 27, 91, 59, 123, 7, 71, 39,
    103, 23, 87, 55, 119, 15, 79, 47, 111, 31, 95, 63, 127,
];
const KYBER_ROOT_OF_UNITY: u64 = 17;
const KYBER_Q: u64 = 3329;

fn init_ntt() {
    let mut tmp = [0; 256];
    let mut zetas = [0; 128];
    let mut zetas_inv = [0; 128];

    tmp[0] = 1;
    for i in 1..256 {
        tmp[i] = (tmp[i - 1] * KYBER_ROOT_OF_UNITY) % KYBER_Q;
    }

    for i in 0..128 {
        zetas[i] = tmp[TREE[i]];
    }

    let mut k = 0;
    let mut i = 64;
    while i >= 1 {
        for j in i..2 * i {
            zetas_inv[k] = KYBER_Q - tmp[128 - TREE[j]];
            k += 1;
        }
        i >>= 1;
    }
}
*/

/// This is used for polynomials over Z_{3329}[X]/(X^256 + 1)
pub fn ntt_kyber(x: &Poly3329<256>) -> Poly3329<256> {
    let mut k = 1;
    let mut len = 128;
    let mut start;
    let mut x = x.coefficients;

    while len >= 2 {
        start = 0;
        while start < 256 {
            let zeta_k = F3329::from(ZETAS[k]);
            k += 1;
            for j in start..start + len {
                let t = zeta_k * x[j + len];
                x[j + len] = x[j] - t;
                x[j] += t;
            }
            start += 2 * len;
        }
        len >>= 1;
    }

    Poly3329::with_coefficients(x)
}

pub fn inv_ntt_kyber(x: &Poly3329<256>) -> Poly3329<256> {
    let mut k = 0;
    let mut len = 2;
    let mut start;
    let mut x = x.coefficients;

    while len <= 128 {
        start = 0;
        while start < 256 {
            let zeta_inv_k = F3329::from(ZETAS_INV[k]);
            k += 1;
            for j in start..start + len {
                let t = x[j];
                x[j] += x[j + len];
                x[j + len] = (t - x[j + len]) * zeta_inv_k;
            }
            start += 2 * len;
        }
        len <<= 1;
    }

    for i in 0..256 {
        x[i] = x[i] * F3329::from(ZETAS_INV[127]);
    }
    Poly3329::with_coefficients(x)
}

/// Basecase multiplication between polynomials
fn base_mul(a: &Poly3329<256>, b: &Poly3329<256>) -> Poly3329<256> {
    // BCM with the zero polynomial is the zero polynomial
    if a.is_zero() || b.is_zero() {
        return Poly3329::ZERO;
    }

    let mut res = Poly3329::ZERO;

    let a = a.coefficients.clone();
    let b = b.coefficients.clone();

    for i in 0..128 {
        let zeta = F3329::from(PRIMITIVE_ZETAS[i]); //zeta = KYBER_ROOT_OF_UNITY^{2*brv_7(i) + 1}

        let mut r0 = a[2 * i] * b[2 * i];
        r0 += a[2 * i + 1] * b[2 * i + 1] * zeta;

        let mut r1 = a[2 * i] * b[2 * i + 1];
        r1 += a[2 * i + 1] * b[2 * i];

        res.set_coefficient(2 * i, r0);
        res.set_coefficient(2 * i + 1, r1);
    }
    res
}

fn vec_inner_prod<const D: usize>(
    a_hat: &PolyVec3329<256, D>,
    b_hat: &PolyVec3329<256, D>,
) -> Poly3329<256> {
    let mut res = base_mul(&a_hat[0], &b_hat[0]);
    for i in 1..D {
        res += base_mul(&a_hat[i], &b_hat[i]);
    }
    res
}

/// Computes a^T*b as NTT^-1(a_hat^T o b_hat)
pub fn inv_ntt_inner_prod<const D: usize>(
    a_hat: &PolyVec3329<256, D>,
    b_hat: &PolyVec3329<256, D>,
) -> Poly3329<256> {
    inv_ntt_kyber(&vec_inner_prod(a_hat, b_hat))
}

/// Computes A_hat * b_hat
pub fn product_matvec<const X: usize, const Y: usize>(
    a_hat: &PolyMatrix3329<256, X, Y>,
    b_hat: &PolyVec3329<256, Y>,
) -> PolyVec3329<256, X> {
    let mut res = PolyVec3329::zero();
    for i in 0..X {
        res[i] = vec_inner_prod(&a_hat.row(i), b_hat);
    }
    res
}

/// Computes A * b as NTT^-1(A_hat * b_hat)
pub fn inv_ntt_product_matvec<const X: usize, const Y: usize>(
    a_hat: &PolyMatrix3329<256, X, Y>,
    b_hat: &PolyVec3329<256, Y>,
) -> PolyVec3329<256, X> {
    let mut res = PolyVec3329::zero();
    for i in 0..X {
        res[i] = vec_inner_prod(&a_hat.row(i), b_hat);
    }
    inv_ntt_vec(&res)
}

/// Number theoretic Transform on vectors
pub fn ntt_vec<const D: usize>(p: &PolyVec3329<256, D>) -> PolyVec3329<256, D> {
    let mut coeffs = [Poly3329::ZERO; D];
    for i in 0..D {
        coeffs[i] = ntt_kyber(&p[i]);
    }
    PolyVec3329::from(coeffs)
}

/// Reverse NTT on vectors
fn inv_ntt_vec<const D: usize>(p_hat: &PolyVec3329<256, D>) -> PolyVec3329<256, D> {
    let mut coeffs = [Poly3329::ZERO; D];
    for i in 0..D {
        coeffs[i] = inv_ntt_kyber(&p_hat[i]);
    }
    PolyVec3329::from(coeffs)
}

#[cfg(test)]
mod tests {
    use crate::math::{finite_field::Fp, prime::get_primitive_root_of_unity};

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
        let mut n = x.len();
        let pow_table = generate_power_table(
            &Fp::<17>::from(get_primitive_root_of_unity(n as u64, 17)),
            n,
        );
        let inv_pow_table = generate_power_table(
            &Fp::<17>::from(get_primitive_root_of_unity(n as u64, 17))
                .inv()
                .unwrap(),
            n,
        );
        ntt_radix2(&mut x, &pow_table);
        let expected_x: Vec<_> = [10, 16, 3, 8, 6, 2, 13, 7]
            .iter()
            .map(|val| Fp::from(*val))
            .collect();
        assert_eq!(x, expected_x);
        inv_ntt_radix2(&mut x, &inv_pow_table);
        assert_eq!(x, original_x);

        let mut x: Vec<_> = [8, 2, 104, 57, 42, 18, 37, 46, 33, 5, 62, 15, 31, 88, 3, 108]
            .iter()
            .map(|val| Fp::from(*val))
            .collect();
        let original_x = x.clone();
        n = x.len();
        let pow_table = generate_power_table(
            &Fp::<113>::from(get_primitive_root_of_unity(n as u64, 113)),
            n,
        );
        let inv_pow_table = generate_power_table(
            &Fp::<113>::from(get_primitive_root_of_unity(n as u64, 113))
                .inv()
                .unwrap(),
            n,
        );
        ntt_radix2(&mut x, &pow_table);
        let expected_x: Vec<_> = [
            94, 75, 17, 29, 21, 78, 105, 105, 94, 99, 94, 39, 21, 5, 108, 48,
        ]
        .iter()
        .map(|val| Fp::from(*val))
        .collect();
        assert_eq!(x, expected_x);
        inv_ntt_radix2(&mut x, &inv_pow_table);
        assert_eq!(x, original_x);
    }

    #[test]
    fn test_ntt_kyber() {
        let mut x = [F3329::ZERO; 256];
        for i in 0..256 {
            x[i] = F3329::from(i as u64);
        }
        let poly = Poly3329::with_coefficients(x);
        let ntt_poly = ntt_kyber(&poly);
        assert_eq!(poly, inv_ntt_kyber(&ntt_poly));
    }

    #[test]
    fn test_kyber_ntt_product() {
        let mut x = [F3329::ZERO; 256];
        for i in 0..256 {
            x[i] = F3329::from(i as u64);
        }
        let a = Poly3329::with_coefficients(x);
        let a_hat = ntt_kyber(&a);

        for i in 0..256 {
            x[i] = F3329::from((i + 1) as u64);
        }
        let b = Poly3329::with_coefficients(x);
        let b_hat = ntt_kyber(&b);

        assert_eq!(a * b, inv_ntt_kyber(&base_mul(&a_hat, &b_hat)));
    }
}
