use crate::math::{finite_field::FiniteRing, Poly3329, F3329};

/// KYBER uses a deterministic approach to sample elements in R_q that
/// are statistically close to a uniformly random distribution.
/// Outputs the NTT-representation a_hat of a \in R_q.
pub(crate) fn parse(byte_stream: &[u8]) -> Poly3329<256> {
    let mut res = Poly3329::ZERO;
    let mut i = 0;
    let mut j = 0;

    while j < 256 {
        let d1 = byte_stream[i] as u64 + 256 * (byte_stream[i + 1] % 16) as u64;
        let d2 = (byte_stream[i + 1] / 16) as u64 + 16 * (byte_stream[i + 2] as u64);
        if d1 < 3329 {
            res.coefficients[j] = d1.into();
            j += 1;
        }
        if d2 < 3329 && j < 256 {
            res.coefficients[j] = d2.into();
            j += 1;
        }
        i += 3;
    }
    res
}

fn get_bit(bytes: &[u8], pos: usize) -> bool {
    let (index, offset) = (pos / 8, pos % 8);
    let mask = 1 << offset;
    (bytes[index] & mask) != 0
}

/// Centered binomial distribution
/// Takes as input an array of 64 * eta bytes
pub(crate) fn cbd(byte_stream: &[u8], eta: usize) -> Poly3329<256> {
    let mut res = Poly3329::ZERO;
    for i in 0..256 {
        let mut a = 0;
        let mut b = 0;
        for j in 0..eta {
            if get_bit(byte_stream, 2 * i * eta + j) {
                a += 1;
            }
            if get_bit(byte_stream, 2 * i * eta + eta + j) {
                b += 1;
            }
        }
        let (a, b) = (F3329::from(a), F3329::from(b));
        res.coefficients[i] = a - b;
    }
    res
}
