use num_bigint::BigUint;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

use super::{constants::INV_3329, prime::{prime_factors, modpow}};

pub trait FiniteRing:
    Sized
    + Eq
    + Add<Output = Self>
    + Neg<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + AddAssign
    + Copy
    + Clone
{
    const ZERO: Self;
    const ONE: Self;
}

pub trait NTT:
    Sized
    + Eq
    + Add<Output = Self>
    + Neg<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + AddAssign
    + Div
    + Copy
    + Clone {
        fn inv(&self) -> Result<Self, String>;
}

pub trait FiniteField: FiniteRing + NTT {
    /// The multiplicative group of the finite field Fp^{*} is cyclic.
    /// Return the generator of Fp^{*}.
    fn get_generator() -> Self;
    // fn inv(&self) -> Result<Self, String>;
    fn order() -> usize;
}

pub struct FpBigUint {
    p: BigUint,
}

impl FpBigUint {
    pub fn new(p: BigUint) -> Self {
        Self { p }
    }

    pub fn add(&self, a: &BigUint, b: &BigUint) -> BigUint {
        (a + b) % &self.p
    }

    pub fn mult(&self, a: &BigUint, b: &BigUint) -> BigUint {
        (a * b) % &self.p
    }

    pub fn neg(&self, a: &BigUint) -> BigUint {
        if *a == BigUint::from(0u32) {
            return a.clone();
        }
        &self.p - a % &self.p
    }

    pub fn subtract(&self, a: &BigUint, b: &BigUint) -> BigUint {
        let b_inv = self.neg(b);
        self.add(a, &b_inv)
    }

    pub fn inv_mult(&self, a: &BigUint) -> BigUint {
        a.modpow(&(&self.p - BigUint::from(2u32)), &self.p)
    }
}


#[derive(PartialEq, Eq, Clone, Copy, Debug)]
pub struct Fp<const P: u64> {
    val: u64,
}

impl<const P: u64> Add for Fp<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            val: (self.val + rhs.val) % P,
        }
    }
}

impl<const P: u64> Mul for Fp<P> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            val: (self.val * rhs.val) % P,
        }
    }
}

impl<const P: u64> Neg for Fp<P> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.val == 0 {
            return Self::ZERO;
        }

        Self { val: P - self.val }
    }
}

impl<const P: u64> Sub for Fp<P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<const P: u64> AddAssign for Fp<P> {
    fn add_assign(&mut self, rhs: Self) {
        self.val = (self.val + rhs.val) % P;
    }
}

impl<const P: u64> Div for Fp<P> {
    type Output = Result<Self, String>;

    fn div(self, rhs: Self) -> Self::Output {
        Ok(self * (rhs.inv()?))
    }
}

impl<const P: u64> FiniteRing for Fp<P> {
    const ZERO: Self = Self { val: 0 };
    const ONE: Self = Self { val: 1 };
}

impl<const P: u64> FiniteField for Fp<P> {
    fn get_generator() -> Self {
        if P == 3329 {
            return Self::from(3);
        }
        for i in 1..P {
            if is_primitive_root(i, P) {
                return Self { val: i };
            }
        }
        unreachable!("P is not a prime number")
    }

    fn order() -> usize {
        P as usize
    }
}

impl<const P: u64> From<u64> for Fp<P> {
    fn from(val: u64) -> Self {
        Self { val: val % P }
    }
}

impl<const P: u64> Into<u64> for Fp<P> {
    fn into(self) -> u64 {
        self.val
    }
}

impl<const P: u64> Fp<P> {
    /// p mod n == 1
    /// get n-th primitive root of unity
    pub fn get_primitive_root_of_unity(n: usize) -> Self {
        if P % (n as u64) != 1 {
            panic!("p mod n == 1 not satisfied");
        }

        if P == 3329 && n == 256 {
            return Self::from(17);
        }

        let g = Self::get_generator().val;
        let power = (P - 1) / n as u64;
        let val = modpow(g, power, P);
        Self { val }
    }
}

impl<const P: u64> NTT for Fp<P> {
    fn inv(&self) -> Result<Self, String> {
        if self == &Self::ZERO {
            return Err("Zero has no inverse".to_string());
        }

        if P == 3329 {
            Ok(INV_3329[self.val as usize].into())
        } else {
            let val = modpow(self.val, P - 2, P);
            Ok(Self { val })
        }
    }
}

fn is_primitive_root(val: u64, p: u64) -> bool {
    if modpow(val, p - 1, p) != 1 {
        return false;
    }
    for i in prime_factors(p - 1) {
        if modpow(val, (p - 1) / i, p) == 1 {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use crate::math::F3329;

    use super::*;

    #[test]
    fn test_add() {
        let mut a = BigUint::from(4u32);
        let b = BigUint::from(10u32);
        let f11 = FpBigUint::new(BigUint::from(11u32));
        let mut res = f11.add(&a, &b);
        assert_eq!(res, BigUint::from(3u32));

        let f41 = FpBigUint::new(BigUint::from(41u32));
        res = f41.add(&a, &b);
        assert_eq!(res, BigUint::from(14u32));

        a = BigUint::from(1u32);
        res = f11.add(&a, &b);
        assert_eq!(res, BigUint::from(0u32));
    }

    #[test]
    fn test_mult() {
        let a = BigUint::from(4u32);
        let b = BigUint::from(10u32);
        let f11 = FpBigUint::new(BigUint::from(11u32));
        let mut res = f11.mult(&a, &b);
        assert_eq!(res, BigUint::from(7u32));

        let f51 = FpBigUint::new(BigUint::from(41u32));
        res = f51.mult(&a, &b);
        assert_eq!(res, BigUint::from(40u32));
    }

    #[test]
    fn test_neg() {
        let mut a = BigUint::from(4u32);
        let f51 = FpBigUint::new(BigUint::from(51u32));
        let mut res = f51.neg(&a);
        assert_eq!(res, BigUint::from(47u32));

        a = BigUint::from(0u32);
        res = f51.neg(&a);
        assert_eq!(res, BigUint::from(0u32));

        a = BigUint::from(52u32);
        res = f51.neg(&a);
        assert_eq!(res, BigUint::from(50u32));
    }

    #[test]
    fn test_subtract() {
        let a = BigUint::from(4u32);
        let f51 = FpBigUint::new(BigUint::from(51u32));
        let res = f51.subtract(&a, &a);
        assert_eq!(res, BigUint::from(0u32));
    }

    #[test]
    fn test_inv_mult() {
        let a = BigUint::from(4u32);
        let f11 = FpBigUint::new(BigUint::from(11u32));
        let a_inv = f11.inv_mult(&a);
        assert_eq!(a_inv, BigUint::from(3u32));

        let res = f11.mult(&a, &a_inv);
        assert_eq!(res, BigUint::from(1u32));
    }

    #[test]
    fn test_generator() {
        assert_eq!(Fp::<17>::get_generator(), Fp::<17>::from(3));
        assert_eq!(Fp::<113>::get_generator(), Fp::<113>::from(3));
        assert_eq!(Fp::<22273>::get_generator(), Fp::<22273>::from(5));
        assert_eq!(Fp::<31489>::get_generator(), Fp::<31489>::from(7));
        assert_eq!(Fp::<26881>::get_generator(), Fp::<26881>::from(11));
        assert_eq!(F3329::get_generator(), F3329::from(3));
    }
}
