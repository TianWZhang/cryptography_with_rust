use crate::math::{finite_field::FiniteRing, Poly3329, PolyVec3329, F3329};

/// Deserialize an array of bytes into a polynomial f
pub(crate) fn decode_to_poly<const N: usize>(data: &[u8]) -> Poly3329<N> {
    assert_eq!(data.len(), N * 8);
    let mut f = [F3329::ZERO; N];
    for i in 0..N {
        f[i] = F3329::from(u64::from_le_bytes(
            data[8 * i..8 * (i + 1)].try_into().unwrap(),
        ));
    }
    Poly3329::with_coefficients(f)
}

/// Serialize Poly into byte array
pub(crate) fn encode_poly<const N: usize>(p: Poly3329<N>) -> Vec<u8> {
    let mut res = vec![];
    for i in 0..N {
        let bytes = p.coefficients[i].to_le_bytes();
        res.extend(bytes)
    }
    res
}

pub(crate) fn decode_to_polyvec<const N: usize, const D: usize>(data: &[u8]) -> PolyVec3329<N, D> {
    assert!(data.len() == 8 * D * N);
    let mut res = PolyVec3329::from([Poly3329::ZERO; D]);
    for i in 0..D {
        res[i] = decode_to_poly(&data[8 * N * i..8 * N * (i + 1)])
    }
    res
}

pub(crate) fn encode_polyvec<const N: usize, const D: usize>(p_vec: PolyVec3329<N, D>) -> Vec<u8> {
    let mut res = vec![];
    for i in 0..D {
        res.extend(encode_poly(p_vec[i]));
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_poly() {
        let mut x = [F3329::ZERO; 256];
        for i in 0..256 {
            x[i] = F3329::from(i as u64);
        }
        let f = Poly3329::with_coefficients(x);
        let encoded_f = encode_poly(f);
        let decoded_f = decode_to_poly::<256>(&encoded_f);
        assert_eq!(f, decoded_f);
    }

    #[test]
    fn test_encode_decode_polyvec() {
        let mut x = [F3329::ZERO; 256];
        for i in 0..256 {
            x[i] = F3329::from(i as u64);
        }
        let f1 = Poly3329::with_coefficients(x);

        for i in 0..256 {
            x[i] = F3329::from((i + 1) as u64);
        }
        let f2 = Poly3329::with_coefficients(x);

        let polys = PolyVec3329::from([f1, f2]);
        let encoded_polys = encode_polyvec(polys);
        let decoded_polys = decode_to_polyvec::<256, 2>(&encoded_polys);
        assert_eq!(polys[0], decoded_polys[0]);
        assert_eq!(polys[1], decoded_polys[1]);
    }
}
