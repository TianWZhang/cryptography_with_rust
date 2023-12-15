/*!
A generic elliptic curve is defined as `y^2 = x^3 + ax + b mod p`, and in this crate we require that
- `p` is a prime number larger than 3
- `4 a^3 + 27 b^2 != 0`
*/

use crate::ecdsa::finite_field;
use num::bigint::BigUint;

use super::finite_field::FiniteField;

#[derive(PartialEq, Clone)]
pub enum Point {
    Coor(BigUint, BigUint),
    Identity,
}

pub enum EllipticCurveError {
    InvalidPoint(Point),
}

#[derive(PartialEq, Clone, Debug)]
pub struct EllipticCurve {
    pub a: BigUint,
    pub b: BigUint,
    pub p: BigUint,
}

impl EllipticCurve {
    pub fn is_on_curve(&self, a: &Point) -> bool {
        match a {
            Point::Identity => true,
            Point::Coor(x, y) => {
                let y2 = y.modpow(&BigUint::from(2u32), &self.p);
                let x3 = x.modpow(&BigUint::from(3u32), &self.p);
                let fp = FiniteField::new(self.p.clone());
                let ax = fp.mult(&self.a, x);
                let x3_plus_ax = fp.add(&x3, &ax);
                y2 == fp.add(&x3_plus_ax, &self.b)
            }
        }
    }

    fn compute_add_coor(
        &self,
        x1: &BigUint,
        y1: &BigUint,
        x2: &BigUint,
        s: &BigUint,
    ) -> (BigUint, BigUint) {
        let s2 = s.modpow(&2u32.into(), &self.p);
        let fp = FiniteField::new(self.p.clone());
        let mut x3 = fp.subtract(&s2, x1);
        x3 = fp.subtract(&x3, x2);

        let mut y3 = fp.subtract(x1, &x3);
        y3 = fp.mult(s, &y3);
        y3 = fp.subtract(&y3, y1);

        (x3, y3)
    }

    pub fn add(&self, a: &Point, b: &Point) -> Result<Point, EllipticCurveError> {
        if !self.is_on_curve(a) {
            return Err(EllipticCurveError::InvalidPoint(a.clone()));
        }

        if !self.is_on_curve(b) {
            return Err(EllipticCurveError::InvalidPoint(b.clone()));
        }

        if *a == *b {
            return Ok(self.double(a));
        }

        Ok(Point::Identity)
    }

    fn double(&self, a: &Point) -> Point {
        match a {
            Point::Identity => Point::Identity,
            Point::Coor(x1, y1) => Point::Identity,
        }
    }
}
