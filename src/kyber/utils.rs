use rand::RngCore;

use crate::{
    hash,
    math::{finite_field::FiniteRing, Poly3329, F3329},
};

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

pub(crate) fn get_bit(bytes: &[u8], pos: usize) -> bool {
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

/// Pseduo random function => SHAKE-256(s||b)
pub(crate) fn prf(s: &[u8], b: usize, len: usize) -> Vec<u8> {
    let mut input = s.to_vec();
    input.extend(b.to_be_bytes());
    hash::shake_256(&input, len)
}

/// Extendable output function => SHAKE-128(rho||i||j) with output of length `len`
pub(crate) fn xof(rho: &[u8], i: usize, j: usize, len: usize) -> Vec<u8> {
    let mut input = rho.to_vec();
    input.extend(i.to_be_bytes());
    input.extend(j.to_be_bytes());
    hash::shake_128(&input, len)
}

/// Hash function SHA3-256
pub(crate) fn h(r: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let hash = hash::sha3_256(r);
    let (p0, p1) = hash.split_at(16);
    (p0.to_vec(), p1.to_vec())
}

/// Hash function SHA3-512
pub(crate) fn g(r: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let hash = hash::sha3_512(r);
    let (p0, p1) = hash.split_at(32);
    (p0.to_vec(), p1.to_vec())
}

/// Hash function SHA3-256
pub(crate) fn kdf(r: &[u8], len: usize) -> Vec<u8> {
    hash::shake_256(r, len)
}

pub(crate) fn random_bytes(len: usize) -> Vec<u8> {
    let mut data = vec![0; len];
    let mut rng = rand::thread_rng();
    rng.fill_bytes(&mut data);
    data
}

pub(crate) fn concat_bytes(items: &[&[u8]]) -> Vec<u8> {
    let len = items.iter().map(|item| item.len()).sum();
    let mut data = Vec::with_capacity(len);
    for item in items {
        data.extend_from_slice(item);
    }
    data
}
