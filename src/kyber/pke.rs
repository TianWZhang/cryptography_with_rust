use crate::math::{
    ntt::{inv_ntt_inner_prod, inv_ntt_product_matvec, ntt_vec, product_matvec},
    PolyMatrix3329, PolyVec3329,
};

use super::{
    compress::{compress_poly, compress_polyvec, decompress_poly, decompress_polyvec},
    encode::{decode_to_poly, decode_to_polyvec, encode_poly, encode_polyvec},
    utils::{cbd, g, parse, prf, random_bytes, xof},
};

/// Default length used for XOF
const XOF_LEN: usize = 4000;
const Q: usize = 3329;
const N: usize = 256;

pub struct PKE<const K: usize> {
    eta1: usize,
    eta2: usize,
    du: usize,
    dv: usize,
}

impl<const K: usize> PKE<K> {
    pub const fn new(eta1: usize, eta2: usize, du: usize, dv: usize) -> Self {
        Self { eta1, eta2, du, dv }
    }

    pub fn key_gen(&self) -> (Vec<u8>, Vec<u8>) {
        let d = random_bytes(32);
        let (rho, sigma) = g(&d);

        let mut a_hat = PolyMatrix3329::<N, K, K>::zero();
        for i in 0..K {
            for j in 0..K {
                a_hat[(i, j)] = parse(&xof(&rho, j, i, XOF_LEN));
            }
        }

        let (mut s, mut e) = (PolyVec3329::<N, K>::zero(), PolyVec3329::<N, K>::zero());
        let prf_len = 64 * self.eta1;

        for i in 0..K {
            s[i] = cbd(&prf(&sigma, i, prf_len), self.eta1);
            e[i] = cbd(&prf(&sigma, K + i, prf_len), self.eta1);
        }

        let (s_hat, e_hat) = (ntt_vec(&s), ntt_vec(&e));
        let t_hat = product_matvec(&a_hat, &s_hat) + e_hat;

        let mut pk = encode_polyvec(t_hat, 12);
        pk.extend_from_slice(&rho);
        let sk = encode_polyvec(s_hat, 12);
        (sk, pk)
    }

    pub fn encrypt(&self, pk: &[u8], msg: &[u8], random_coins: &[u8]) -> Vec<u8> {
        let offset = 12 * K * N / 8;
        let (t, rho) = pk.split_at(offset);

        let t_hat = decode_to_polyvec(t, 12);
        let mut a_t = PolyMatrix3329::zero();

        // generate A_hat \in R_q^{K \times K} in NTT domain
        for i in 0..K {
            for j in 0..K {
                a_t[(i, j)] = parse(&xof(rho, i, j, XOF_LEN));
            }
        }

        // sample r \in R_q^K from B_{eta_1}, sample e_1 \in R_q^K from B_{eta_2}
        let (mut r, mut e1) = (PolyVec3329::<N, K>::zero(), PolyVec3329::<N, K>::zero());
        for i in 0..K {
            r[i] = cbd(&prf(random_coins, i, 64 * self.eta1), self.eta1);
            e1[i] = cbd(&prf(random_coins, K + i, 64 * self.eta2), self.eta2);
        }
        // sample e_1 \in R_q from B_{eta_2}
        let e2 = cbd(&prf(random_coins, 2 * K, 64 * self.eta2), self.eta2);

        let r_hat = ntt_vec(&r);

        // u = A^T * r + e_1
        let u = inv_ntt_product_matvec(&a_t, &r_hat) + e1;
        // v = t^T * r + e_2 + Decompress_q(msg, 1)
        let v =
            inv_ntt_inner_prod(&t_hat, &r_hat) + e2 + decompress_poly(decode_to_poly(msg, 1), 1, Q);

        let mut c1 = encode_polyvec(compress_polyvec(u, self.du, Q), self.du);
        let c2 = encode_poly(compress_poly(v, self.dv, Q), self.dv);
        c1.extend(c2);
        c1
    }

    pub fn decrypt(&self, sk: &[u8], ciphertext: &[u8]) -> Vec<u8> {
        let offset = (self.du as usize) * K * N / 8;
        let (c1, c2) = ciphertext.split_at(offset);

        let u = decompress_polyvec(decode_to_polyvec::<K>(c1, self.du), self.du, Q);
        let v = decompress_poly(decode_to_poly(c2, self.dv), self.dv, Q);
        let s = decode_to_polyvec(sk, 12);

        let u_hat = ntt_vec(&u);
        let x = inv_ntt_inner_prod(&s, &u_hat);
        let p = v - x;
        encode_poly(compress_poly(p, 1, Q), 1)
    }
}

#[cfg(test)]
mod tests {
    use crate::kyber::{KYBER1024PKE, KYBER512PKE, KYBER768PKE};

    use super::*;

    #[test]
    fn test_pke512() {
        let pke = KYBER512PKE;
        let (sk, pk) = pke.key_gen();

        let msg = random_bytes(32);
        let random_coins = random_bytes(32);

        let ciphertext = pke.encrypt(&pk, &msg, &random_coins);
        let dec_msg = pke.decrypt(&sk, &ciphertext);

        assert_eq!(msg, dec_msg);
    }

    #[test]
    fn test_pke768() {
        let pke = KYBER768PKE;
        let (sk, pk) = pke.key_gen();

        let msg = random_bytes(32);
        let random_coins = random_bytes(32);

        let ciphertext = pke.encrypt(&pk, &msg, &random_coins);
        let dec_msg = pke.decrypt(&sk, &ciphertext);

        assert_eq!(msg, dec_msg);
    }

    #[test]
    fn test_pke1024() {
        let pke = KYBER1024PKE;
        let (sk, pk) = pke.key_gen();

        let msg = random_bytes(32);
        let random_coins = random_bytes(32);

        let ciphertext = pke.encrypt(&pk, &msg, &random_coins);
        let dec_msg = pke.decrypt(&sk, &ciphertext);

        assert_eq!(msg, dec_msg);
    }
}
