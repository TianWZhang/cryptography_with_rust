use crate::math::{finite_field::FiniteRing, Poly3329, PolyVec3329, F3329};

use super::utils::get_bit;

/// Deserialize an array of bytes into a polynomial f
pub(crate) fn decode_to_poly(data: &[u8], l: usize) -> Poly3329<256> {
    let mut f = [F3329::ZERO; 256];
    for i in 0..256 {
        for j in 0..l {
            if get_bit(data, i * l + j) {
                f[i] += F3329::from(1 << j);
            }
        }
    }
    Poly3329::with_coefficients(f)
}

/// Serialize Poly into byte array
pub(crate) fn encode_poly(p: Poly3329<256>, l: usize) -> Vec<u8> {
    let mut res = vec![];
    let mut c: u8 = 0;

    for i in 0..256 {
        let mut v: u64 = p.coefficients[i].into();
        for j in 0..l {
            let s = (i * l + j) % 8;
            if s == 0 && !(i == 0 && j == 0) {
                res.push(c);
                c = 0;
            }
            if v & 1 == 1 {
                let a = 1 << s;
                c += a as u8;
            }
            v >>= 1;
        }
    }
    res.push(c);
    res
}

pub(crate) fn decode_to_polyvec<const D: usize>(data: &[u8], ell: usize) -> PolyVec3329<256, D> {
    let mut res = PolyVec3329::from([Poly3329::ZERO; D]);
    for i in 0..D {
        res[i] = decode_to_poly(&data[32 * ell * i..32 * ell * (i + 1)], ell)
    }
    res
}

pub(crate) fn encode_polyvec<const D: usize>(p_vec: PolyVec3329<256, D>, ell: usize) -> Vec<u8> {
    let mut res = vec![];
    for i in 0..D {
        res.extend(encode_poly(p_vec[i], ell));
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
        let encoded_f = encode_poly(f, 12);
        let decoded_f = decode_to_poly(&encoded_f, 12);
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
        let encoded_polys = encode_polyvec(polys, 12);
        let decoded_polys = decode_to_polyvec::<2>(&encoded_polys, 12);
        assert_eq!(polys[0], decoded_polys[0]);
        assert_eq!(polys[1], decoded_polys[1]);
    }
}
