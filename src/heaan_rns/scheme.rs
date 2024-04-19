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
        let ax = self.context.sample_uniform(l, 0);

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

        let mut sx_sqaure = self.context.mul(&sk.sx, &sk.sx, l, 0);
        self.context.mul_by_bigp(&mut sx_sqaure, l);

        let mut ex = self.context.sample_guass(l, k);
        self.context.ntt(&mut ex, l, k);
        // here we set k = 0 because after multiplying sx_square with P,
        // it's zero mod pi for every special prime
        self.context.add_inplace(&mut ex, &sx_sqaure, l, 0);

        let ax = self.context.sample_uniform(l, k);
        let mut bx = self.context.mul(&ax, &sk.sx, l, k);
        self.context.sub2_inplace(&ex, &mut bx, l, k);

        self.key_map
            .insert(KeyType::Multiplication, Key::new(ax, bx));
    }

    #[allow(dead_code)]
    fn add_conjugation_key(&mut self, sk: &Sk) {
        let l = self.context.max_level;
        let k = self.context.num_special_modulus;

        let mut sx_conj = self.context.conjugate(&sk.sx, l);
        self.context.mul_by_bigp(&mut sx_conj, l);

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

    pub fn encode_single(&self, v: Complex64, l: usize) -> Plaintext {
        let msg = self.context.encode_single(v, l);
        Plaintext { msg, slots: 1, l }
    }

    pub fn decode_single(&self, plaintext: &Plaintext) -> Complex64 {
        self.context.decode_single(&plaintext.msg)
    }

    pub fn encrypt_plaintext(&self, plaintext: &Plaintext) -> Ciphertext {
        let key = self.key_map.get(&KeyType::Encryption).unwrap();

        let mut vx = self.context.sample_zo(plaintext.l, 0);
        self.context.ntt(&mut vx, plaintext.l, 0);

        let mut e1x = self.context.sample_guass(plaintext.l, 0);
        self.context.ntt(&mut e1x, plaintext.l, 0);

        // ct.ax = vx * key.ax + e1x
        let mut ax = self.context.mul(&vx, &key.ax, plaintext.l, 0);
        self.context.add_inplace(&mut ax, &e1x, plaintext.l, 0);

        // ct.bx =  vx * key.bx + e2x + msg
        let mut bx = self.context.mul(&vx, &key.bx, plaintext.l, 0);
        let mut e2x = self.context.sample_guass(plaintext.l, 0);
        self.context.ntt(&mut e2x, plaintext.l, 0);
        self.context.add_inplace(&mut bx, &e2x, plaintext.l, 0);
        self.context
            .add_inplace(&mut bx, &plaintext.msg, plaintext.l, 0);

        Ciphertext {
            ax,
            bx,
            n: self.context.n as usize,
            slots: plaintext.slots,
            l: plaintext.l,
        }
    }

    pub fn decrypt_ciphertext(&self, sk: &Sk, ciphertext: &Ciphertext) -> Plaintext {
        let mut msg = self.context.mul(&ciphertext.ax, &sk.sx, 1, 0);
        self.context.add_inplace(&mut msg, &ciphertext.bx, 1, 0);
        Plaintext {
            msg,
            slots: ciphertext.slots,
            l: 1,
        }
    }

    pub fn encrypt(&self, vals: &Vec<Complex64>, l: usize) -> Ciphertext {
        let plaintext = self.encode(vals, l);
        self.encrypt_plaintext(&plaintext)
    }

    pub fn encrypt_single(&self, v: Complex64, l: usize) -> Ciphertext {
        let plaintext = self.encode_single(v, l);
        self.encrypt_plaintext(&plaintext)
    }

    pub fn decrypt(&self, sk: &Sk, ciphertext: &Ciphertext) -> Vec<Complex64> {
        let plaintext = self.decrypt_ciphertext(sk, ciphertext);
        self.decode(&plaintext)
    }

    pub fn decrypt_single(&self, sk: &Sk, ciphertext: &Ciphertext) -> Complex64 {
        let plaintext = self.decrypt_ciphertext(sk, ciphertext);
        self.decode_single(&plaintext)
    }

    pub fn add(&self, ct1: &Ciphertext, ct2: &Ciphertext) -> Ciphertext {
        assert_eq!(ct1.l, ct2.l);
        assert_eq!(ct1.slots, ct2.slots);

        Ciphertext {
            ax: self.context.add(&ct1.ax, &ct2.ax, ct1.l, 0),
            bx: self.context.add(&ct1.bx, &ct2.bx, ct1.l, 0),
            n: ct1.n,
            slots: ct1.slots,
            l: ct1.l,
        }
    }

    pub fn add_inplace(&self, ct1: &mut Ciphertext, ct2: &Ciphertext) {
        assert_eq!(ct1.l, ct2.l);
        assert_eq!(ct1.slots, ct2.slots);
        self.context.add_inplace(&mut ct1.ax, &ct2.ax, ct1.l, 0);
        self.context.add_inplace(&mut ct1.bx, &ct2.bx, ct1.l, 0);
    }

    pub fn negate(&self, ct: &Ciphertext) -> Ciphertext {
        Ciphertext {
            ax: self.context.negate(&ct.ax, ct.l, 0),
            bx: self.context.negate(&ct.bx, ct.l, 0),
            n: ct.n,
            slots: ct.slots,
            l: ct.l,
        }
    }

    pub fn negate_inplace(&self, ct: &mut Ciphertext) {
        self.context.negate_inplace(&mut ct.ax, ct.l, 0);
        self.context.negate_inplace(&mut ct.bx, ct.l, 0);
    }

    pub fn sub(&self, ct1: &Ciphertext, ct2: &Ciphertext) -> Ciphertext {
        assert_eq!(ct1.l, ct2.l);
        assert_eq!(ct1.slots, ct2.slots);

        Ciphertext {
            ax: self.context.sub(&ct1.ax, &ct2.ax, ct1.l, 0),
            bx: self.context.sub(&ct1.bx, &ct2.bx, ct1.l, 0),
            n: ct1.n,
            slots: ct1.slots,
            l: ct1.l,
        }
    }

    pub fn sub_inplace(&self, ct1: &mut Ciphertext, ct2: &Ciphertext) {
        assert_eq!(ct1.l, ct2.l);
        assert_eq!(ct1.slots, ct2.slots);
        self.context.sub_inplace(&mut ct1.ax, &ct2.ax, ct1.l, 0);
        self.context.sub_inplace(&mut ct1.bx, &ct2.bx, ct1.l, 0);
    }

    pub fn sub2_inplace(&self, ct1: &Ciphertext, ct2: &mut Ciphertext) {
        assert_eq!(ct1.l, ct2.l);
        assert_eq!(ct1.slots, ct2.slots);
        self.context.sub2_inplace(&ct1.ax, &mut ct2.ax, ct1.l, 0);
        self.context.sub2_inplace(&ct1.bx, &mut ct2.bx, ct1.l, 0);
    }

    pub fn rescale_by(&self, ct: &Ciphertext, dl: usize) -> Ciphertext {
        let mut res = ct.clone();
        for _ in 0..dl {
            let ax = self.context.rescale(&res.ax, res.l);
            let bx = self.context.rescale(&res.bx, res.l);
            res = Ciphertext {
                ax,
                bx,
                n: res.n,
                slots: res.slots,
                l: res.l - 1,
            };
        }
        res
    }

    pub fn rescale_to(&self, ct: &Ciphertext, l: usize) -> Ciphertext {
        let dl = ct.l - l;
        self.rescale_by(ct, dl)
    }

    pub fn mul(&self, ct1: &Ciphertext, ct2: &Ciphertext) -> Ciphertext {
        assert_eq!(ct1.l, ct2.l);
        let ax1bx2 = self.context.mul(&ct1.ax, &ct2.bx, ct1.l, 0);
        let ax2bx1 = self.context.mul(&ct2.ax, &ct1.bx, ct2.l, 0);

        let axax = self.context.mul(&ct1.ax, &ct2.ax, ct1.l, 0);
        let axax = self.context.mod_up(&axax, ct1.l);
        let bxbx = self.context.mul(&ct1.bx, &ct2.bx, ct1.l, 0);

        let mul_key = self.key_map.get(&KeyType::Multiplication).unwrap();
        // ax_mul.len() == ct1.n * (ct1.l + self.context.num_special_modulus)
        let ax_mul = self.context.mul_key(&axax, &mul_key.ax, ct1.l);
        let bx_mul = self.context.mul_key(&axax, &mul_key.bx, ct1.l);

        let mut ax_mul = self.context.approx_modulus_reduction(&ax_mul, ct1.l);
        // ax_mul.len() == ct1.n * ct1.l
        let mut bx_mul = self.context.approx_modulus_reduction(&bx_mul, ct1.l);

        // ax_mul = ax1 * bx2 + ax2 * bx1 + axax * mul_key.ax
        self.context.add_inplace(&mut ax_mul, &ax1bx2, ct1.l, 0);
        self.context.add_inplace(&mut ax_mul, &ax2bx1, ct1.l, 0);

        // bx_mul = bx1 * bx2 + axax * mul_key.bx
        self.context.add_inplace(&mut bx_mul, &bxbx, ct1.l, 0);

        Ciphertext {
            ax: ax_mul,
            bx: bx_mul,
            n: ct1.n,
            slots: ct1.slots,
            l: ct1.l,
        }
    }

    pub fn mul_inplace(&self, ct1: &mut Ciphertext, ct2: &Ciphertext) {
        assert_eq!(ct1.l, ct2.l);
        let ax1bx2 = self.context.mul(&ct1.ax, &ct2.bx, ct1.l, 0);
        let ax2bx1 = self.context.mul(&ct2.ax, &ct1.bx, ct2.l, 0);

        let axax = self.context.mul(&ct1.ax, &ct2.ax, ct1.l, 0);
        let axax = self.context.mod_up(&axax, ct1.l);
        let bxbx = self.context.mul(&ct1.bx, &ct2.bx, ct1.l, 0);

        let mul_key = self.key_map.get(&KeyType::Multiplication).unwrap();
        let ax_mul = self.context.mul_key(&axax, &mul_key.ax, ct1.l);
        let bx_mul = self.context.mul_key(&axax, &mul_key.bx, ct1.l);
        let mut ax_mul = self.context.approx_modulus_reduction(&ax_mul, ct1.l);
        let mut bx_mul = self.context.approx_modulus_reduction(&bx_mul, ct1.l);

        // ax_mul = ax1 * bx2 + ax2 * bx1 + axax * mul_key.ax
        self.context.add_inplace(&mut ax_mul, &ax1bx2, ct1.l, 0);
        self.context.add_inplace(&mut ax_mul, &ax2bx1, ct1.l, 0);

        // bx_mul = bx1 * bx2 + axax * mul_key.bx
        self.context.add_inplace(&mut bx_mul, &bxbx, ct1.l, 0);

        ct1.ax = ax_mul;
        ct1.bx = bx_mul;
    }

    pub fn square(&self, ct: &Ciphertext) -> Ciphertext {
        let ct1 = ct.clone();
        self.mul(&ct, &ct1)
    }

    pub fn square_inplace(&self, ct: &mut Ciphertext) {
        let ct1 = ct.clone();
        self.mul_inplace(ct, &ct1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::heaan_rns::utils::{equal_up_to_epsilon, gen_random_complex_vector};

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
        let ciphertext = scheme.encrypt_plaintext(&plaintext);
        let decrypted = scheme.decrypt_ciphertext(&scheme.sk, &ciphertext);
        let v_decoded = scheme.decode(&decrypted);
        assert!(equal_up_to_epsilon(&v, &v_decoded, 0.000000000001));
    }

    #[test]
    fn test_homomorphic_mul() {
        let l = 2;
        let k = 2;
        let slots = 8;
        let context = Context::new(1 << 15, l, k, 1 << 55, 3.2);
        let scheme = HeaanRnsScheme::new(context);

        let v1 = gen_random_complex_vector(slots);
        let v2 = gen_random_complex_vector(slots);
        let v_mul: Vec<_> = v1.iter().zip(v2.iter()).map(|(c1, c2)| c1 * c2).collect();

        let ct1 = scheme.encrypt(&v1, l);
        let ct2 = scheme.encrypt(&v2, l);
        let ct_mul = scheme.mul(&ct1, &ct2);
        let ct_mul = scheme.rescale_by(&ct_mul, 1);

        let v_mul_decrypted = scheme.decrypt(&scheme.sk, &ct_mul);
        assert!(equal_up_to_epsilon(&v_mul, &v_mul_decrypted, 0.0000001));
    }

    #[test]
    fn test_homomorphic_add() {
        let l = 3;
        let k = 4;
        let slots = 8;
        let context = Context::new(1 << 15, l, k, 1 << 55, 3.2);
        let scheme = HeaanRnsScheme::new(context);

        let v1 = gen_random_complex_vector(slots);
        let v2 = gen_random_complex_vector(slots);
        let v_add: Vec<_> = v1.iter().zip(v2.iter()).map(|(c1, c2)| c1 + c2).collect();

        let ct1 = scheme.encrypt(&v1, l);
        let ct2 = scheme.encrypt(&v2, l);
        let ct_add = scheme.add(&ct1, &ct2);

        let v_add_decrypted = scheme.decrypt(&scheme.sk, &ct_add);
        assert!(equal_up_to_epsilon(&v_add, &v_add_decrypted, 0.0000001));
    }
}
