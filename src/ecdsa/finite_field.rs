use std::ops::Rem;
use num_bigint::BigUint;

pub struct FiniteField {
    p: BigUint,
}

impl FiniteField {
    pub fn new(p: BigUint) -> Self {
        Self { p }
    }

    pub fn add(&self, a: &BigUint, b: &BigUint) -> BigUint {
        (a + b).rem(&self.p)
    }

    pub fn mult(&self, a: &BigUint, b: &BigUint) -> BigUint {
        (a * b).rem(&self.p)
    }

    pub fn inv_add(&self, a: &BigUint) -> BigUint {
        if *a == BigUint::from(0u32) {
            return a.clone();
        }
        // &self.p - a.mod_floor(&self.p)
        &self.p - a.rem(&self.p)
    }

    pub fn subtract(&self, a: &BigUint, b: &BigUint) -> BigUint {
        let b_inv = self.inv_add(b);
        self.add(a, &b_inv)
    }

    pub fn inv_mult(&self, a: &BigUint) -> BigUint {
        a.modpow(&(&self.p - BigUint::from(2u32)), &self.p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        let mut a = BigUint::from(4u32);
        let b = BigUint::from(10u32);
        let f11 = FiniteField::new(BigUint::from(11u32));
        let mut res = f11.add(&a, &b);
        assert_eq!(res, BigUint::from(3u32));

        let f41 = FiniteField::new(BigUint::from(41u32));
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
        let f11 = FiniteField::new(BigUint::from(11u32));
        let mut res = f11.mult(&a, &b);
        assert_eq!(res, BigUint::from(7u32));

        let f51 = FiniteField::new(BigUint::from(41u32));
        res = f51.mult(&a, &b);
        assert_eq!(res, BigUint::from(40u32));
    }

    #[test]
    fn test_inv_add() {
        let mut a = BigUint::from(4u32);
        let f51 = FiniteField::new(BigUint::from(51u32));
        let mut res = f51.inv_add(&a);
        assert_eq!(res, BigUint::from(47u32));

        a = BigUint::from(0u32);
        res = f51.inv_add(&a);
        assert_eq!(res, BigUint::from(0u32));

        a = BigUint::from(52u32);
        res = f51.inv_add(&a);
        assert_eq!(res, BigUint::from(50u32));
    }

    #[test]
    fn test_subtract() {
        let a = BigUint::from(4u32);
        let f51 = FiniteField::new(BigUint::from(51u32));
        let res = f51.subtract(&a, &a);
        assert_eq!(res, BigUint::from(0u32));
    }

    #[test]
    fn test_inv_mult() {
        let a = BigUint::from(4u32);
        let f11 = FiniteField::new(BigUint::from(11u32));
        let a_inv = f11.inv_mult(&a);
        assert_eq!(a_inv, BigUint::from(3u32));

        let res = f11.mult(&a, &a_inv);
        assert_eq!(res, BigUint::from(1u32));
    }
}
