use crate::math::{finite_field::FiniteRing, Poly3329, PolyVec3329, F3329};

fn compress_int(x: u64, d: usize, q: usize) -> u64 {
    let m = 1 << d;
    ((((m as f64) / (q as f64)) * (x as f64)).round() as u64) % m
}

fn decompress_int(x: u64, d: usize, q: usize) -> u64 {
    let m = 1 << d;
    (((q as f64) / (m as f64)) * (x as f64)).round() as u64
}

pub(crate) fn compress_poly<const N: usize>(x: Poly3329<N>, d: usize, q: usize) -> Poly3329<N> {
    let mut coeffs = [F3329::ZERO; N];
    for i in 0..N {
        coeffs[i] = F3329::from(compress_int(x.coefficients[i].into(), d, q));
    }
    Poly3329::with_coefficients(coeffs)
}

pub(crate) fn decompress_poly<const N: usize>(x: Poly3329<N>, d: usize, q: usize) -> Poly3329<N> {
    let mut coeffs = [F3329::ZERO; N];
    for i in 0..N {
        coeffs[i] = F3329::from(decompress_int(x.coefficients[i].into(), d, q));
    }
    Poly3329::with_coefficients(coeffs)
}

pub(crate) fn compress_polyvec<const N: usize, const D: usize>(
    x: PolyVec3329<N, D>,
    d: usize,
    q: usize,
) -> PolyVec3329<N, D> {
    let mut res = [Poly3329::ZERO; D];
    for i in 0..D {
        res[i] = compress_poly(x[i], d, q);
    }
    PolyVec3329::from(res)
}

pub(crate) fn decompress_polyvec<const N: usize, const D: usize>(
    x: PolyVec3329<N, D>,
    d: usize,
    q: usize,
) -> PolyVec3329<N, D> {
    let mut res = [Poly3329::ZERO; D];
    for i in 0..D {
        res[i] = decompress_poly(x[i], d, q);
    }
    PolyVec3329::from(res)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compress_decompress_poly() {
        let mut x = [F3329::ZERO; 256];
        for i in 0..256 {
            x[i] = F3329::from(i as u64);
        }
        let f = Poly3329::with_coefficients(x);
        let compressed_f = compress_poly(f, 12, 3329);
        let uncompressed_f = decompress_poly(compressed_f, 12, 3329);
        assert_eq!(f, uncompressed_f);
    }
}
