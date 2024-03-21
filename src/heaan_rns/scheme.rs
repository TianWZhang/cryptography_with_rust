use super::{
    context::Context,
    key::{Key, KeyType, Sk},
    plaintext::{Ciphertext, Plaintext},
};
use num_complex::Complex64;
use std::collections::HashMap;

pub struct HeaanRnsScheme {
    pub context: Context,
    /// contains encryption, multiplication and conjugation keys
    pub key_map: HashMap<KeyType, Key>,
    /// contains left rotation keys
    pub left_rot_key_map: HashMap<KeyType, Key>,
    pub sk: Sk,
}

impl HeaanRnsScheme {
    pub fn new(context: Context) -> Self {
        let sk = Sk::new(&context);
        let mut res = Self {
            context,
            key_map: HashMap::new(),
            left_rot_key_map: HashMap::new(),
            sk: sk.clone(),
        };
        res.add_encryption_key(&sk);
        res.add_multiplication_key(&sk);
        res
    }

    fn add_encryption_key(&mut self, sk: &Sk) {
        let l = self.context.max_level;
        // let ax = self.context.sample_uniform(l, 0);
        let ax = vec![0; (self.context.n as usize) * l];

        let mut ex = self.context.sample_guass(l, 0);
        self.context.ntt(&mut ex, l, 0);

        // bx = ex - ax * sx
        let mut bx = self.context.mul(&ax, &sk.sx, l, 0);
        self.context.sub2_inplace(&ex, &mut bx, l, 0);

        self.key_map.insert(KeyType::Encryption, Key::new(ax, bx));
    }

    fn add_multiplication_key(&mut self, sk: &Sk) {
        let l = self.context.max_level;
        let k = self.context.num_special_modulus;

        // TODO: eval and equal?
        let sx_sqaure = self.context.mul(&sk.sx, &sk.sx, l, 0);

        let mut ex = self.context.sample_guass(l, k);
        self.context.ntt(&mut ex, l, k);
        self.context.add_inplace(&mut ex, &sx_sqaure, l, 0);

        let ax = self.context.sample_uniform(l, k);
        let mut bx = self.context.mul(&ax, &sk.sx, l, k);
        self.context.sub2_inplace(&ex, &mut bx, l, k);

        self.key_map
            .insert(KeyType::Multiplication, Key::new(ax, bx));
    }

    fn add_conjugation_key(&mut self, sk: &Sk) {
        let l = self.context.max_level;
        let k = self.context.num_special_modulus;

        // TODO: eval and equal?
        let sx_conj = self.context.conjugate(&sk.sx, l);

        let mut ex = self.context.sample_guass(l, k);
        self.context.ntt(&mut ex, l, k);
        self.context.add_inplace(&mut ex, &sx_conj, l, 0);

        let ax = self.context.sample_uniform(l, k);
        let mut bx = self.context.mul(&ax, &sk.sx, l, k);
        self.context.sub2_inplace(&ex, &mut bx, l, k);

        self.key_map.insert(KeyType::Conjugation, Key::new(ax, bx));
    }

    pub fn encode(&self, v: &Vec<Complex64>, l: usize) -> Plaintext {
        let msg = self.context.encode(v, l);
        Plaintext {
            msg,
            slots: v.len(),
            l,
        }
    }

    pub fn decode(&self, plaintext: &Plaintext) -> Vec<Complex64> {
        self.context.decode(&plaintext.msg, plaintext.slots)
    }

    pub fn encrypt(&self, msg: &Plaintext) -> Ciphertext {
        let key = self.key_map.get(&KeyType::Encryption).unwrap();

        let mut vx = self.context.sample_zo(msg.l, 0);
        self.context.ntt(&mut vx, msg.l, 0);

        let mut e1x = self.context.sample_guass(msg.l, 0);
        self.context.ntt(&mut e1x, msg.l, 0);

        // ct.ax = vx * key.ax + e1x
        let mut ax = self.context.mul(&vx, &key.ax, msg.l, 0);
        self.context.add_inplace(&mut ax, &e1x, msg.l, 0);

        // ct.bx =  vx * key.bx + e2x + msg
        let mut bx = self.context.mul(&vx, &key.bx, msg.l, 0);
        let mut e2x = self.context.sample_guass(msg.l, 0);
        self.context.ntt(&mut e2x, msg.l, 0);
        self.context.add_inplace(&mut bx, &e2x, msg.l, 0);
        self.context.add_inplace(&mut bx, &msg.msg, msg.l, 0);

        Ciphertext {
            ax,
            bx,
            n: self.context.n as usize,
            slots: msg.slots,
            l: msg.l,
        }
    }

    pub fn decrypt(&self, sk: &Sk, ciphertext: &Ciphertext) -> Plaintext {
        let mut msg = self.context.mul(&ciphertext.ax, &sk.sx, 1, 0);
        self.context.add_inplace(&mut msg, &ciphertext.bx, 1, 0);
        Plaintext {
            msg,
            slots: ciphertext.slots,
            l: 1,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::heaan_rns::utils::equal_up_to_epsilon;
    use super::*;

    #[test]
    fn test_encrypt_then_decrypt() {
        let l = 1;
        let k = 1;
        let context = Context::new(1 << 15, l, k, 1 << 55, 3.2);
        let scheme = HeaanRnsScheme::new(context);

        let v = vec![
            Complex64::new(0.47, 0.97),
            Complex64::new(0.12, 0.77),
            Complex64::new(-0.45, 0.37),
            Complex64::new(0.08, 0.39),
            Complex64::new(0.44, -0.98),
            Complex64::new(0.19, 0.98),
            Complex64::new(-0.12, -0.44),
            Complex64::new(0.20, 0.24),
        ];
        let plaintext = scheme.encode(&v, l);
        let ciphertext = scheme.encrypt(&plaintext);
        let decrypted = scheme.decrypt(&scheme.sk, &ciphertext);
        let v_decoded = scheme.decode(&decrypted);
        assert!(equal_up_to_epsilon(&v, &v_decoded, 0.000000000001));
    }
}
