use std::ops::{Add, Sub, Mul, Div, Neg, AddAssign};

use crate::math::{prime::modpow, finite_field::NTT};

/// This struct is useful if the value of prime modulus p is only 
/// known at run time.
#[derive(PartialEq, Eq, Clone, Copy, Debug)]
pub struct Fp_ {
    val: u64,
    p: u64
}

impl Add for Fp_ {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        if self.p != rhs.p {
            panic!("The two numbers must have the same modulus p");
        }
        Self {
            val: (self.val + rhs.val) % self.p,
            p: self.p
        }
    }
}

impl Mul for Fp_ {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.p != rhs.p {
            panic!("The two numbers must have the same modulus p");
        }
        Self {
            val: (self.val * rhs.val) % self.p,
            p: self.p
        }
    }
}

impl Neg for Fp_ {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.val == 0 {
            return self;
        }

        Self { val: self.p - self.val, p: self.p }
    }
}

impl Sub for Fp_ {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl AddAssign for Fp_ {
    fn add_assign(&mut self, rhs: Self) {
        if self.p != rhs.p {
            panic!("The two numbers must have the same modulus p");
        }
        self.val = (self.val + rhs.val) % self.p;
    }
}

impl Div for Fp_ {
    type Output = Result<Self, String>;

    fn div(self, rhs: Self) -> Self::Output {
        Ok(self * (rhs.inv()?))
    }
}

impl NTT for Fp_ {
    // fn get_generator() -> Self {
    //     for i in 1..self.p {
    //         if is_primitive_root(i, self.p) {
    //             return Self { val: i, p: self.p };
    //         }
    //     }
    //     unreachable!("P is not a prime number")
    // }

    fn inv(&self) -> Result<Self, String> {
        if self.val == 0 {
            return Err("Zero has no inverse".to_string());
        }

        let val = modpow(self.val, self.p - 2, self.p);
        Ok(Self { val, p: self.p })

    }
}

// impl Fp_ {
//     /// p mod n == 1
//     /// get n-th primitive root of unity
//     pub fn get_primitive_root_of_unity(n: usize) -> Self {
//         if P % (n as u64) != 1 {
//             panic!("p mod n == 1 not satisfied");
//         }

//         if P == 3329 && n == 256 {
//             return Self::from(17);
//         }

//         let g = Self::get_generator().val;
//         let power = (P - 1) / n as u64;
//         let val = modpow(g, power, P);
//         Self { val }
//     }
// }
