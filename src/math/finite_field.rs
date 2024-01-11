use num_bigint::BigUint;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

use super::constants::INV_3329;

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

pub trait FiniteField: FiniteRing + Div {
    fn inv(&self) -> Result<Self, String>;
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
            val: add_mod(self.val, rhs.val, P),
        }
    }
}

impl<const P: u64> Mul for Fp<P> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            val: mul_mod(self.val, rhs.val, P),
        }
    }
}

impl<const P: u64> Neg for Fp<P> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.val == 0 {
            return Self::ZERO;
        }

        Self {
            val: neg_mod(self.val, P),
        }
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
    fn inv(&self) -> Result<Self, String> {
        if self == &Self::ZERO {
            return Err("Zero has no inverse".to_string());
        }

        if P == 3329 {
            Ok(INV_3329[self.val as usize].into())
        } else {
            let val = BigUint::from(self.val)
                .modpow(&BigUint::from(P - 2), &BigUint::from(P))
                .try_into()
                .unwrap();
            Ok(Self { val })
        }
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

#[inline]
pub(crate) fn add_mod(a1: u64, a2: u64, p: u64) -> u64 {
    ((a1 as u128 + a2 as u128) % p as u128) as u64
}

#[inline]
pub(crate) fn mul_mod(a1: u64, a2: u64, p: u64) -> u64 {
    ((a1 as u128 * a2 as u128) % p as u128) as u64
}

#[inline]
pub(crate) fn neg_mod(mut a: u64, p: u64) -> u64 {
    a = a % p;
    if a == 0 {
        0
    } else {
        p - a
    }
}

#[inline]
pub(crate) fn sub_mod(a1: u64, a2: u64, p: u64) -> u64 {
    let a2_neg = neg_mod(a2, p);
    add_mod(a1, a2_neg, p)
}

#[inline]
pub(crate) fn inv_mod(mut a: u64, p: u64) -> u64 {
    a = a % p;
    if a == 0 {
        panic!("Zero has no inverse");
    }
    BigUint::from(a)
        .modpow(&BigUint::from(p - 2), &BigUint::from(p))
        .try_into()
        .unwrap()
}

pub(crate) fn div_mod(a1: u64, a2: u64, p: u64) -> u64 {
    let a2_inv = inv_mod(a2, p);
    mul_mod(a1, a2_inv, p)
}

#[cfg(test)]
mod tests {
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
}
