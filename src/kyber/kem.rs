use super::{
    pke::PKE,
    utils::{concat_bytes, g, h, kdf, random_bytes},
};

pub struct KEM<const K: usize> {
    pke: PKE<K>,
    sk_size: usize,
}

impl<const K: usize> KEM<K> {
    pub fn key_gen(&self) -> (Vec<u8>, Vec<u8>) {
        let z = random_bytes(32);

        let (sk_prime, pk) = self.pke.key_gen();
        let (h1, h2) = h(&pk);
        let sk = concat_bytes(&[&sk_prime, &pk, &h1, &h2, &z]);

        (sk, pk)
    }

    /// Encryption : public key  => ciphertext, Shared Key
    pub fn encapsulate(&self, pk: &[u8]) -> (Vec<u8>, Vec<u8>) {
        let m = random_bytes(32);
        let (mut m1, m2) = h(&m);
        let (h1, h2) = h(pk);
        let (k_bar, r) = g(&concat_bytes(&[&m1, &m2, &h1, &h2]));

        m1.extend(m2);
        let c = self.pke.encrypt(pk, &m1, &r);

        let (h1, h2) = h(&c);
        let k = kdf(&concat_bytes(&[&k_bar, &h1, &h2]), self.sk_size);

        (c, k)
    }

    /// Decryption : secret key, ciphertext => Shared Key
    pub fn decapsulate(&self, c: &[u8], sk: &[u8]) -> Vec<u8> {
        // Spliting sk = (sk'||pk||H(pk)||z)
        let (sk_prime, rem) = sk.split_at(12 * K * 256 / 8);
        let (pk, rem) = rem.split_at(12 * K * 256 / 8 + 32);
        let (hash, z) = rem.split_at(32);

        let m_prime = self.pke.decrypt(&sk_prime, c);
        let (k_bar, r_prime) = g(&concat_bytes(&[&m_prime, &hash]));
        let c_prime = self.pke.encrypt(&pk, &m_prime, &r_prime);

        let (h1, h2) = h(c);
        if *c == c_prime {
            kdf(&concat_bytes(&[&k_bar, &h1, &h2]), self.sk_size)
        } else {
            kdf(&concat_bytes(&[&z, &h1, &h2]), self.sk_size)
        }
    }

    pub const fn new(pke: PKE<K>, sk_size: usize) -> Self {
        Self { pke, sk_size }
    }
}

#[cfg(test)]
mod tests {
    use crate::kyber::{KYBER1024KEM, KYBER512KEM, KYBER768KEM};

    #[test]
    fn encapsulate_then_decapsulate_ccakem_512() {
        let kem = KYBER512KEM;

        let (sk, pk) = kem.key_gen();
        let (ctx, shk) = kem.encapsulate(&pk);
        let shk_dec = kem.decapsulate(&ctx, &sk);
        assert_eq!(shk, shk_dec);
    }

    #[test]
    fn encapsulate_then_decapsulate_ccakem_768() {
        let kem = KYBER768KEM;

        let (sk, pk) = kem.key_gen();
        let (ctx, shk) = kem.encapsulate(&pk);
        let shk_dec = kem.decapsulate(&ctx, &sk);
        assert_eq!(shk, shk_dec);
    }

    #[test]
    fn encapsulate_then_decapsulate_ccakem_1024() {
        let kem = KYBER1024KEM;

        let (sk, pk) = kem.key_gen();
        let (ctx, shk) = kem.encapsulate(&pk);
        let shk_dec = kem.decapsulate(&ctx, &sk);
        assert_eq!(shk, shk_dec);
    }
}
