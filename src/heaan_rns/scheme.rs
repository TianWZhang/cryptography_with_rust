use super::{
    context::Context,
    key::{Key, KeyType, Sk},
    plaintext::{Ciphertext, Plaintext},
};
use num_complex::Complex64;
use std::collections::HashMap;

pub const SIGMOID: [f64; 11] = [ 1.0 / 2.0, 1.0 / 4.0, 0.0, -1.0 / 48.0, 0.0, 1.0 / 480.0, 0.0, -17.0 / 80640.0, 0.0, 31.0 / 1451520.0, 0.0 ];
pub const EXPONENT: [f64; 11] = [ 1.0, 1.0, 0.5, 1.0 / 6.0, 1.0 / 24.0, 1.0 / 120.0, 1.0 / 720.0, 1.0 / 5040.0, 1.0 / 40320.0, 1.0 / 362880.0, 1.0 / 3628800.0 ];

pub struct HeaanRnsScheme {
    pub context: Context,
    /// contains encryption, multiplication and conjugation keys
    pub key_map: HashMap<KeyType, Key>,
    /// contains left rotation keys
    pub left_rot_key_map: HashMap<usize, Key>,
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

    pub fn add_conjugation_key(&mut self) {
        let l = self.context.max_level;
        let k = self.context.num_special_modulus;

        let mut sx_conj = self.context.conjugate(&self.sk.sx, l);
        self.context.mul_by_bigp(&mut sx_conj, l);

        let mut ex = self.context.sample_guass(l, k);
        self.context.ntt(&mut ex, l, k);
        self.context.add_inplace(&mut ex, &sx_conj, l, 0);

        let ax = self.context.sample_uniform(l, k);
        let mut bx = self.context.mul(&ax, &self.sk.sx, l, k);
        self.context.sub2_inplace(&ex, &mut bx, l, k);

        self.key_map.insert(KeyType::Conjugation, Key::new(ax, bx));
    }

    pub fn add_left_rot_key(&mut self, rot: usize) {
        let l = self.context.max_level;
        let k = self.context.num_special_modulus;

        let mut sx_rot = self.context.left_rot(&self.sk.sx, l, rot);
        self.context.mul_by_bigp(&mut sx_rot, l);

        let mut ex = self.context.sample_guass(l, k);
        self.context.ntt(&mut ex, l, k);
        self.context.add_inplace(&mut ex, &sx_rot, l, 0);

        let ax = self.context.sample_uniform(l, k);
        let mut bx = self.context.mul(&ax, &self.sk.sx, l, k);
        self.context.sub2_inplace(&ex, &mut bx, l, k);

        self.left_rot_key_map.insert(rot, Key::new(ax, bx));
    }

    /// add left rotation keys for rot = 1, 2, 2^2, 2^3, ..., slots
    pub fn add_left_rot_keys(&mut self, slots: usize) {
        let mut i = 1;
        while i <= slots {
            if self.left_rot_key_map.get(&i).is_none() {
                self.add_left_rot_key(i);
            }
            i *= 2;
        }
    }

    pub fn add_right_rot_keys(&mut self, slots: usize) {
        let mut i = 1;
        while i <= slots {
            let idx = slots - i;
            if self.left_rot_key_map.get(&idx).is_none() {
                self.add_left_rot_key(idx);
            }
            i *= 2;
        }
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

    // (ax, bx) bx + ax * sx ->
    // (ax_rot, bx_rot) bx_rot + ax_rot * sx_rot -> keyswitch
    // ax_rot * (key.ax, key.bx) + (0, bx_rot) ->phase: bx_rot + ax_rot * key.bx + ax_rot * key.ax * sx
    // = bx_rot + ax_rot * (key.bx + key.ax * sx) = bx_rot + ax_rot * sx_rot
    fn left_rotate_by_pow2_inplace(&self, ct: &mut Ciphertext, rot_slots: usize) {
        let key = self
            .left_rot_key_map
            .get(&rot_slots)
            .expect("please add left rotation key first");

        let ax_rot = self.context.left_rot(&ct.ax, ct.l, rot_slots);
        let bx_rot = self.context.left_rot(&ct.bx, ct.l, rot_slots);
        let ax_rot = self.context.mod_up(&ax_rot, ct.l);

        let ax = self.context.mul_key(&ax_rot, &key.ax, ct.l);
        let ax = self.context.approx_modulus_reduction(&ax, ct.l);
        let bx = self.context.mul_key(&ax_rot, &key.bx, ct.l);
        let mut bx = self.context.approx_modulus_reduction(&bx, ct.l);
        self.context.add_inplace(&mut bx, &bx_rot, ct.l, 0);

        ct.ax = ax;
        ct.bx = bx;
    }

    pub fn left_rotate_inplace(&self, ct: &mut Ciphertext, rot_slots: usize) {
        let rot_slots = rot_slots % ct.slots;
        let log_rot_slots = (rot_slots as f64).log2() as usize;
        for i in 0..=log_rot_slots {
            if (rot_slots >> i) & 1 == 1 {
                self.left_rotate_by_pow2_inplace(ct, 1 << i);
            }
        }
    }

    pub fn right_rotate_inplace(&self, ct: &mut Ciphertext, rot_slots: usize) {
        let rot_slots = ct.slots - (rot_slots % ct.slots);
        self.left_rotate_inplace(ct, rot_slots);
    }

    pub fn left_rotate(&self, ct: &Ciphertext, rot_slots: usize) -> Ciphertext {
        let mut res = ct.clone();
        self.left_rotate_inplace(&mut res, rot_slots);
        res
    }

    pub fn right_rotate(&self, ct: &Ciphertext, rot_slots: usize) -> Ciphertext {
        let rot_slots = ct.slots - (rot_slots % ct.slots);
        self.left_rotate(ct, rot_slots)
    }

    // (ax, bx) ->phase bx + ax * sx
    // (ax_conj, bx_conj)
    // ax_conj * key + (0, bx_conj) ->phase ax_conj * (key.bx + key.ax * sx_conj) + bx_conj = bx_conj + ax_conj * sx_conj
    pub fn conjugate_inplace(&self, ct: &mut Ciphertext) {
        self.context.conjugate_inplace(&mut ct.ax, ct.l);
        self.context.conjugate_inplace(&mut ct.bx, ct.l);
        let key = self
            .key_map
            .get(&KeyType::Conjugation)
            .expect("please add conjugation key first");

        let ax_conj = self.context.mod_up(&ct.ax, ct.l);
        let ax_conj_times_key_ax = self.context.mul_key(&ax_conj, &key.ax, ct.l);
        let ax_conj_times_key_bx = self.context.mul_key(&ax_conj, &key.bx, ct.l);

        let ax = self
            .context
            .approx_modulus_reduction(&ax_conj_times_key_ax, ct.l);
        let bx = self
            .context
            .approx_modulus_reduction(&ax_conj_times_key_bx, ct.l);
        ct.ax = ax;
        self.context.add_inplace(&mut ct.bx, &bx, ct.l, 0);
    }

    pub fn conjugate(&self, ct: &Ciphertext) -> Ciphertext {
        let ax_conj = self.context.conjugate(&ct.ax, ct.l);
        let bx_conj = self.context.conjugate(&ct.bx, ct.l);
        let key = self
            .key_map
            .get(&KeyType::Conjugation)
            .expect("please add conjugation key first");

        let ax_conj = self.context.mod_up(&ax_conj, ct.l);
        let ax_conj_times_key_ax = self.context.mul_key(&ax_conj, &key.ax, ct.l);
        let ax_conj_times_key_bx = self.context.mul_key(&ax_conj, &key.bx, ct.l);

        let ax = self
            .context
            .approx_modulus_reduction(&ax_conj_times_key_ax, ct.l);
        let mut bx = self
            .context
            .approx_modulus_reduction(&ax_conj_times_key_bx, ct.l);
        self.context.add_inplace(&mut bx, &bx_conj, ct.l, 0);

        Ciphertext {
            ax,
            bx,
            n: ct.n,
            slots: ct.slots,
            l: ct.l,
        }
    }

    pub fn imult(&self, ct: &Ciphertext) -> Ciphertext {
        let ax = self
            .context
            .mul_by_xpow(&ct.ax, ct.l, self.context.n as usize / 2);
        let bx = self
            .context
            .mul_by_xpow(&ct.bx, ct.l, self.context.n as usize / 2);
        Ciphertext {
            ax,
            bx,
            n: ct.n,
            slots: ct.slots,
            l: ct.l,
        }
    }

    pub fn add_const(&self, ct: &Ciphertext, constant: f64) -> Ciphertext {
        let const_encoded = (constant.abs() * self.context.p as f64) as u64;
        let ax = ct.ax.clone();
        let bx = if constant >= 0.0 {
            self.context.add_const(&ct.bx, const_encoded, ct.l, 0)
        } else {
            self.context.sub_const(&ct.bx, const_encoded, ct.l, 0)
        };
        Ciphertext {
            ax,
            bx,
            n: ct.n,
            slots: ct.slots,
            l: ct.l,
        }
    }

    pub fn add_const_inplace(&self, ct: &mut Ciphertext, constant: f64) {
        let const_encoded = (constant.abs() * self.context.p as f64) as u64;
        if constant >= 0.0 {
            self.context
                .add_const_inplace(&mut ct.bx, const_encoded, ct.l, 0);
        } else {
            self.context
                .sub_const_inplace(&mut ct.bx, const_encoded, ct.l, 0);
        }
    }

    pub fn mul_const(&self, ct: &Ciphertext, constant: f64) -> Ciphertext {
        let const_encoded = (constant.abs() * self.context.p as f64) as u64;
        let mut ax = self.context.mul_const(&ct.ax, const_encoded, ct.l, 0);
        let mut bx = self.context.mul_const(&ct.bx, const_encoded, ct.l, 0);
        if constant <= 0.0 {
            self.context.negate_inplace(&mut ax, ct.l, 0);
            self.context.negate_inplace(&mut bx, ct.l, 0);
        }

        Ciphertext {
            ax,
            bx,
            n: ct.n,
            slots: ct.slots,
            l: ct.l,
        }
    }

    pub fn mul_const_inplace(&self, ct: &mut Ciphertext, constant: f64) {
        let const_encoded = (constant.abs() * self.context.p as f64) as u64;
        self.context
            .mul_const_inplace(&mut ct.ax, const_encoded, ct.l, 0);
        self.context
            .mul_const_inplace(&mut ct.bx, const_encoded, ct.l, 0);
        if constant <= 0.0 {
            self.context.negate_inplace(&mut ct.ax, ct.l, 0);
            self.context.negate_inplace(&mut ct.bx, ct.l, 0);
        }
    }

    pub fn mod_down_by(&self, ct: &Ciphertext, dl: usize) -> Ciphertext {
        assert!(dl < ct.l);
        let ax = self.context.mod_down_by(&ct.ax, ct.l, dl);
        let bx = self.context.mod_down_by(&ct.bx, ct.l, dl);
        Ciphertext {
            ax,
            bx,
            n: ct.n,
            slots: ct.slots,
            l: ct.l - dl,
        }
    }

    pub fn mod_down_by_inplace(&self, ct: &mut Ciphertext, dl: usize) {
        assert!(dl < ct.l);
        self.context.mod_down_by_inplace(&mut ct.ax, ct.l, dl);
        self.context.mod_down_by_inplace(&mut ct.bx, ct.l, dl);
        ct.l -= dl;
    }

    // (a0, a1, a2, a3) -> (a0 + a1, a1 + a2, a2 + a3, a3 + a0) -> (sum, sum, sum, sum)
    pub fn slots_sum(&self, ct: &Ciphertext) -> Ciphertext {
        let mut res = ct.clone();
        let mut i = 1;
        while i < ct.slots {
            let rot = self.left_rotate(&res, i);
            self.add_inplace(&mut res, &rot);
            i <<= 1;
        }
        res
    }

    fn power_of_2s(&self, ct: &Ciphertext, log_degree: usize) -> Vec<Ciphertext> {
        let mut res = Vec::with_capacity(log_degree + 1);
        res.push(ct.clone());
        for i in 1..=log_degree {
            let tmp = self.square(&res[i - 1]);
            let tmp = self.rescale_by(&tmp, 1);
            res.push(tmp);
        }
        res
    }

    fn powers(&self, ct: &Ciphertext, degree: usize) -> Vec<Ciphertext> {
        let mut res = Vec::with_capacity(degree);
        let log_degree = (degree as f64).log2() as usize;
        let cpow2s = self.power_of_2s(ct, log_degree);

        for i in 0..log_degree {
            // initially, idx = 2^i - 1
            // idx \in [2^i - 1, 2^{i + 1} - 1)
            let powi = 1 << i;
            res.push(cpow2s[i].clone());
            for j in 0..(powi - 1) {
                // res[j] = ct^{j+1}
                let mut tmp = self.mod_down_by(&res[j], res[j].l - cpow2s[i].l);
                self.mul_inplace(&mut tmp, &cpow2s[i]);
                // res[idx] = ct^{j+1+2^i}, idx == 2^i + j
                tmp = self.rescale_by(&tmp, 1);
                res.push(tmp);
            }
        }
        res.push(cpow2s[log_degree].clone());
        let deg2 = 1 << log_degree;
        for i in 0..(degree - deg2) {
            let mut tmp = self.mod_down_by(&res[i], res[i].l - cpow2s[log_degree].l);
            self.mul_inplace(&mut tmp, &cpow2s[log_degree]);
            tmp = self.rescale_by(&tmp, 1);
            res.push(tmp);
        }
        res
    }

    pub fn homomorphic_function(&self, ct: &Ciphertext, coeffs: &[f64], degree: usize) -> Ciphertext {
        let cpows = self.powers(ct, degree);
        let mut res = self.mul_const(&cpows[0], coeffs[1]);

        for i in 1..degree {
            if coeffs[i + 1].abs() > 1e-27 {
                let aixi = self.mul_const(&cpows[i], coeffs[i + 1]);
                let oldl = res.l;
                self.mod_down_by_inplace(&mut res, oldl - aixi.l);
                self.add_inplace(&mut res, &aixi);
            }
        }
        let mut res = self.rescale_by(&res, 1);
        self.add_const_inplace(&mut res, coeffs[0]);
        res
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

    #[test]
    fn test_homomorphic_rotate() {
        let l = 3;
        let k = l;
        let slots = 8;
        let context = Context::new(1 << 15, l, k, 1 << 55, 3.2);
        let mut scheme = HeaanRnsScheme::new(context);
        scheme.add_left_rot_keys(slots);
        scheme.add_right_rot_keys(slots);

        let rot_slots = 9;
        let mut v = gen_random_complex_vector(slots);
        let ct = scheme.encrypt(&v, l);
        let ct_rotate = scheme.left_rotate(&ct, rot_slots);
        let v_rot_decrypted = scheme.decrypt(&scheme.sk, &ct_rotate);
        v.rotate_left(rot_slots % slots);
        assert!(equal_up_to_epsilon(&v, &v_rot_decrypted, 0.0000001));

        let rot_slots = 3;
        let mut v = gen_random_complex_vector(slots);
        let mut ct = scheme.encrypt(&v, l);
        scheme.right_rotate_inplace(&mut ct, rot_slots);
        let v_rot_decrypted = scheme.decrypt(&scheme.sk, &ct);
        v.rotate_right(rot_slots % slots);
        assert!(equal_up_to_epsilon(&v, &v_rot_decrypted, 0.0000001));
    }

    #[test]
    fn test_homomorphic_conjugation() {
        let l = 2;
        let k = l;
        let slots = 8;
        let context = Context::new(1 << 15, l, k, 1 << 55, 3.2);
        let mut scheme = HeaanRnsScheme::new(context);
        scheme.add_conjugation_key();

        let v = gen_random_complex_vector(slots);
        let v_conj: Vec<_> = v.iter().map(|c| c.conj()).collect();

        let ct = scheme.encrypt(&v, l);
        let ct_conj = scheme.conjugate(&ct);
        let v_conj_decrypted = scheme.decrypt(&scheme.sk, &ct_conj);
        assert!(equal_up_to_epsilon(&v_conj, &v_conj_decrypted, 0.0000001));

        let v = gen_random_complex_vector(slots);
        let v_conj: Vec<_> = v.iter().map(|c| c.conj()).collect();
        let mut ct = scheme.encrypt(&v, l);
        scheme.conjugate_inplace(&mut ct);
        let v_conj_decrypted = scheme.decrypt(&scheme.sk, &ct);
        assert!(equal_up_to_epsilon(&v_conj, &v_conj_decrypted, 0.0000001));
    }

    #[test]
    fn test_imult() {
        let l = 2;
        let k = l;
        let slots = 8;
        let context = Context::new(1 << 15, l, k, 1 << 55, 3.2);
        let scheme = HeaanRnsScheme::new(context);

        let v = gen_random_complex_vector(slots);
        let i = Complex64::new(0.0, 1.0);
        let v_times_i: Vec<_> = v.iter().map(|c| c * i).collect();

        let ct = scheme.encrypt(&v, l);
        let ct_times_i = scheme.imult(&ct);
        let v_times_i_decrypted = scheme.decrypt(&scheme.sk, &ct_times_i);
        assert!(equal_up_to_epsilon(
            &v_times_i,
            &v_times_i_decrypted,
            0.0000001
        ));
    }

    #[test]
    fn test_add_const() {
        let l = 2;
        let k = l;
        let slots = 8;
        let context = Context::new(1 << 15, l, k, 1 << 55, 3.2);
        let scheme = HeaanRnsScheme::new(context);

        let v = gen_random_complex_vector(slots);
        let constant = 0.67;
        let v_add_const: Vec<_> = v.iter().map(|c| c + constant).collect();
        let ct = scheme.encrypt(&v, l);
        let ct_add_const = scheme.add_const(&ct, constant);
        let v_add_const_decrypted = scheme.decrypt(&scheme.sk, &ct_add_const);
        assert!(equal_up_to_epsilon(
            &v_add_const,
            &v_add_const_decrypted,
            0.0000001
        ));

        let v = gen_random_complex_vector(slots);
        let constant = -0.34;
        let v_add_const: Vec<_> = v.iter().map(|c| c + constant).collect();
        let mut ct = scheme.encrypt(&v, l);
        scheme.add_const_inplace(&mut ct, constant);
        let v_add_const_decrypted = scheme.decrypt(&scheme.sk, &ct);
        assert!(equal_up_to_epsilon(
            &v_add_const,
            &v_add_const_decrypted,
            0.0000001
        ));
    }

    #[test]
    fn test_mul_const() {
        let l = 2;
        let k = l;
        let slots = 8;
        let context = Context::new(1 << 15, l, k, 1 << 55, 3.2);
        let scheme = HeaanRnsScheme::new(context);

        let v = gen_random_complex_vector(slots);
        let constant = 0.67;
        let v_mul_const: Vec<_> = v.iter().map(|c| c * constant).collect();
        let ct = scheme.encrypt(&v, l);
        let ct_mul_const = scheme.mul_const(&ct, constant);
        let ct_mul_const = scheme.rescale_by(&ct_mul_const, 1);
        let v_mul_const_decrypted = scheme.decrypt(&scheme.sk, &ct_mul_const);
        assert!(equal_up_to_epsilon(
            &v_mul_const,
            &v_mul_const_decrypted,
            0.0000001
        ));

        let v = gen_random_complex_vector(slots);
        let constant = -0.34;
        let v_mul_const: Vec<_> = v.iter().map(|c| c * constant).collect();
        let mut ct = scheme.encrypt(&v, l);
        scheme.mul_const_inplace(&mut ct, constant);
        let ct = scheme.rescale_by(&ct, 1);
        let v_mul_const_decrypted = scheme.decrypt(&scheme.sk, &ct);
        assert!(equal_up_to_epsilon(
            &v_mul_const,
            &v_mul_const_decrypted,
            0.0000001
        ));
    }

    #[test]
    fn test_slots_sum() {
        let l = 2;
        let k = l;
        let slots = 8;
        let context = Context::new(1 << 15, l, k, 1 << 55, 3.2);
        let mut scheme = HeaanRnsScheme::new(context);
        scheme.add_left_rot_keys(slots);

        let v = gen_random_complex_vector(slots);
        let slots_sum = v.iter().fold(Complex64::new(0.0, 0.0), |acc, x| acc + x);
        let v_slots_sum = vec![slots_sum; slots];

        let ct = scheme.encrypt(&v, l);
        let ct_slots_sum = scheme.slots_sum(&ct);
        let v_slots_sum_decrypted = scheme.decrypt(&scheme.sk, &ct_slots_sum);
        assert!(equal_up_to_epsilon(
            &v_slots_sum,
            &v_slots_sum_decrypted,
            0.0000001
        ));
    }

    #[test]
    fn test_homomorphic_function() {
        let l = 8;
        let k = l + 1;
        let slots = 8;
        let context = Context::new(1 << 15, l, k, 1 << 55, 3.2);
        let scheme = HeaanRnsScheme::new(context);

        let v = gen_random_complex_vector(slots);
        let v_sigmoid: Vec<_> = v.iter().map(|c| c.exp() / (1.0 + c.exp())).collect();
        let ct = scheme.encrypt(&v, l);
        let ct_sigmoid = scheme.homomorphic_function(&ct, &SIGMOID, 5);
        let v_sigmoid_decrypted = scheme.decrypt(&scheme.sk, &ct_sigmoid);
        assert!(equal_up_to_epsilon(
            &v_sigmoid,
            &v_sigmoid_decrypted,
            0.01
        ));

        let v = gen_random_complex_vector(slots);
        let v_exp: Vec<_> = v.iter().map(|c| c.exp()).collect();
        let ct = scheme.encrypt(&v, l);
        let ct_exp = scheme.homomorphic_function(&ct, &EXPONENT, 5);
        let v_exp_decrypted = scheme.decrypt(&scheme.sk, &ct_exp);
        assert!(equal_up_to_epsilon(
            &v_exp,
            &v_exp_decrypted,
            0.01
        ));
    }
}
