/*!
A generic elliptic curve is defined as `y^2 = x^3 + ax + b mod p`, and in this crate we require that
- `p` is a prime number larger than 3
- `4 a^3 + 27 b^2 != 0`
*/

use std::ops::Rem;
use super::finite_field::FpBigUint;
use num_bigint::BigUint;

#[derive(PartialEq, Clone, Debug)]
pub enum Point {
    Coor(BigUint, BigUint),
    Identity
}

#[derive(Debug)]
pub enum PointError {
    CompressIdentity,
    UnsupportedPrime,
    InvalidCompressedPoint
}

#[derive(Debug)]
pub enum EllipticCurveError {
    InvalidPoint(Point),
    InvalidScalar(BigUint)
}

#[derive(PartialEq, Clone, Debug)]
pub struct EllipticCurve {
    pub a: BigUint,
    pub b: BigUint,
    pub p: BigUint,
}

impl EllipticCurve {
    fn discriminant(a: &BigUint, b: &BigUint) -> BigUint {
        let a_cubic  = a.pow(3);
        let b_sqaure = b.pow(2);
        BigUint::from(4u32) * a_cubic + BigUint::from(27u32) * b_sqaure
    }

    pub fn new(a: BigUint, b: BigUint, p: BigUint) -> Self {
        if Self::discriminant(&a, &b) == BigUint::from(0u32) {
            panic!("The elliptic curve is not smooth");
        }
        Self { a, b, p }
    }

    pub fn compress(point: &Point) -> Result<String, PointError> {
        match point {
            Point::Identity => Err(PointError::CompressIdentity),
            Point::Coor(x, y) => {
                let mut prefix = if y.rem(2u8) == 0u32.into() {
                    String::from("02")
                } else {
                    String::from("03")
                };
                let hex_x = format!("{:x}", x);
                prefix.push_str(&hex_x);
                Ok(prefix)
            }
        }
    }

    /// We require p mod 4 == 3 as
    /// there is an efficient sqaqure root algorithm for such prime p
    pub fn uncompress(&self, point: &str) -> Result<Point, PointError> {
        if (&self.p).rem(4u8) != 3u32.into() {
            return Err(PointError::UnsupportedPrime);
        }
        if !point.starts_with("02") && !point.starts_with("03") {
            return Err(PointError::InvalidCompressedPoint);
        }

        let x = BigUint::parse_bytes(
            &point[2..].as_bytes(),
            16,
        );
        if x.is_none() {
            return Err(PointError::InvalidCompressedPoint);
        }

        let x = x.unwrap();
        let x_cubic = x.modpow(&3u32.into(), &self.p);
        let fp = FpBigUint::new(self.p.clone());
        let ax = fp.mult(&self.a, &x);
        let mut y_square = fp.add(&x_cubic, &ax);
        y_square = fp.add(&y_square, &self.b);
        let pe = (&self.p + BigUint::from(1u32)) / 4u32;
        let mut y = y_square.modpow(&pe, &self.p);

        let sign_y = BigUint::parse_bytes(&point[..2].as_bytes(), 16).unwrap() - 2u32;
        // if the parity does not match, y should be the other root
        if (&y).rem(2u32) != sign_y {
            y = fp.neg(&y);
        }
        Ok(Point::Coor(x, y))
    }

    pub fn is_on_curve(&self, a: &Point) -> bool {
        match a {
            Point::Identity => true,
            Point::Coor(x, y) => {
                let y2 = y.modpow(&BigUint::from(2u32), &self.p);
                let x3 = x.modpow(&BigUint::from(3u32), &self.p);
                let fp = FpBigUint::new(self.p.clone());
                let ax = fp.mult(&self.a, x);
                let x3_plus_ax = fp.add(&x3, &ax);
                y2 == fp.add(&x3_plus_ax, &self.b)
            }
        }
    }

    /// x_3 = s^2 - x_1 - x_2
    /// y_3 = s * (x_1 - x_3) - y_1
    fn compute_add_coor(
        &self,
        x1: &BigUint,
        y1: &BigUint,
        x2: &BigUint,
        s: &BigUint,
    ) -> (BigUint, BigUint) {
        let s_sqaure = s.modpow(&2u32.into(), &self.p);
        let fp = FpBigUint::new(self.p.clone());
        let mut x3 = fp.subtract(&s_sqaure, x1);
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

        match (a, b) {
            (Point::Identity, _) => Ok(b.clone()),
            (_, Point::Identity) => Ok(a.clone()),
            (Point::Coor(x1, y1), Point::Coor(x2, y2)) => {
                let fp = FpBigUint::new(self.p.clone());
                // Check whether they are additive inverse to each other
                if x1 == x2 && fp.add(y1, y2) == 0u32.into() {
                    return Ok(Point::Identity);
                }

                // s = (y2 - y1) / (x2 - x1) mod p
                let numerator = fp.subtract(y2, y1);
                let mut denominator_inv = fp.subtract(x2, x1);
                denominator_inv = fp.inv_mult(&denominator_inv);
                let s = fp.mult(&numerator, &denominator_inv);
                let (x3, y3) = self.compute_add_coor(x1, y1, x2, &s);
                Ok(Point::Coor(x3, y3))
            }
        }
    }

    fn double(&self, a: &Point) -> Point {
        match a {
            Point::Identity => Point::Identity,
            Point::Coor(x1, y1) => {
                if *y1 == 0u32.into() {
                    return Point::Identity;
                }

                let fp = FpBigUint::new(self.p.clone());
                // s = (3 * x1 ^ 2 + a) / (2 * y1) mod p
                let mut numerator = fp.mult(x1, x1);
                numerator = fp.mult(&3u32.into(), &numerator);
                numerator = fp.add(&self.a, &numerator);

                let denominator = fp.mult(&2u32.into(), y1);
                let denominator_inv = fp.inv_mult(&denominator);
                let s = fp.mult(&numerator, &denominator_inv);
                let (x3, y3) = self.compute_add_coor(x1, y1, x1, &s);
                Point::Coor(x3, y3)
            }
        }
    }

    pub fn scalar_mult(&self, a: &Point, d: &BigUint) -> Result<Point, EllipticCurveError> {
        if *d == 0u32.into() {
            return Err(EllipticCurveError::InvalidScalar(d.clone()));
        }

        let mut res = a.clone();
        for i in (0..d.bits() - 1).rev() {
            res = self.double(&res);
            if d.bit(i) {
                res = self.add(&res, a)?;
            }
        }
        Ok(res)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_on_curve() {
        // y^2 = x^3 + 2x + 2 mod 17
        let ec = EllipticCurve::new(2u32.into(), 2u32.into(), 17u32.into());
        // (6, 3) + (5, 1) = (10, 6)
        let p1 = Point::Coor(6u32.into(), 3u32.into());
        let p2 = Point::Coor(5u32.into(), 1u32.into());
        let p3 = Point::Coor(10u32.into(), 6u32.into());
        assert!(ec.is_on_curve(&p1));
        assert!(ec.is_on_curve(&p2));
        assert!(ec.is_on_curve(&p3));

        let p4 = Point::Coor(4u32.into(), 1u32.into());
        let p5 = Point::Coor(1u32.into(), 1u32.into());
        let p6 = Point::Coor(0u32.into(), 1u32.into());
        assert!(!ec.is_on_curve(&p4));
        assert!(!ec.is_on_curve(&p5));
        assert!(!ec.is_on_curve(&p6));
    }

    #[test]
    fn test_point_addition() {
        // y^2 = x^3 + 2x + 2 mod 17
        let ec = EllipticCurve::new(2u32.into(), 2u32.into(), 17u32.into());

        // (6,3) + (5,1) = (10,6)
        let mut p1 = Point::Coor(BigUint::from(6u32), BigUint::from(3u32));
        let mut p2 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        let mut pr = Point::Coor(BigUint::from(10u32), BigUint::from(6u32));

        let mut res = ec.add(&p1, &p2).unwrap();
        assert_eq!(res, pr);
        res = ec.add(&p2, &p1).unwrap();
        assert_eq!(res, pr);

        // (5,1) + (5,1) = (6,3)
        p1 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        p2 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        pr = Point::Coor(BigUint::from(6u32), BigUint::from(3u32));
        res = ec.add(&p1, &p2).unwrap();
        assert_eq!(res, pr);

        // (10, 6) + (5, 1) = (3,1)
        p1 = Point::Coor(BigUint::from(10u32), BigUint::from(6u32));
        p2 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        pr = Point::Coor(BigUint::from(3u32), BigUint::from(1u32));
        res = ec.add(&p1, &p2).unwrap();
        assert_eq!(res, pr);

        // (16, 13) + (5, 1) = (0, 6)
        p1 = Point::Coor(BigUint::from(16u32), BigUint::from(13u32));
        p2 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        pr = Point::Coor(BigUint::from(0u32), BigUint::from(6u32));
        res = ec.add(&p1, &p2).unwrap();
        assert_eq!(res, pr);

        // (6,3) + I = (6,3)
        p1 = Point::Coor(BigUint::from(6u32), BigUint::from(3u32));
        res = ec.add(&p1, &Point::Identity).unwrap();
        assert_eq!(res, p1.clone());

        res = ec.add(&Point::Identity, &p1).unwrap();
        assert_eq!(res, p1.clone());
        // I + I = 2 * I = I
        res = ec.add(&Point::Identity, &Point::Identity).unwrap();
        assert_eq!(res, Point::Identity);

        // (5,16) + (5,1) = I
        p1 = Point::Coor(BigUint::from(5u32), BigUint::from(16u32));
        p2 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        res = ec.add(&p1, &p2).unwrap();
        assert_eq!(res, Point::Identity);

        res = ec.add(&p2, &p1).unwrap();
        assert_eq!(res, Point::Identity);
    }

    #[test]
    fn test_point_doubling() {
        // y^2 = x^3 + 2x + 2 mod 17
        let ec = EllipticCurve::new(2u32.into(), 2u32.into(), 17u32.into());

        // (5,1) + (5,1) = 2 (5, 1) = (6,3)
        let p1 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        let pr = Point::Coor(BigUint::from(6u32), BigUint::from(3u32));
        let mut res = ec.double(&p1);
        assert_eq!(res, pr);

        // I + I = 2 * I = I
        res = ec.double(&Point::Identity);
        assert_eq!(res, Point::Identity);
    }

    #[test]
    fn test_scalar_mult() {
        // y^2 = x^3 + 2x + 2 mod 17   |G| = 19  19 * A = I
        let ec = EllipticCurve::new(2u32.into(), 2u32.into(), 17u32.into());

        let a = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));

        // 2 * (5, 1) = (6,3)
        let mut pr = Point::Coor(BigUint::from(6u32), BigUint::from(3u32));
        let mut res = ec.scalar_mult(&a, &BigUint::from(2u32)).unwrap();
        assert_eq!(res, pr);

        // 10 * (5, 1) = (7,11)
        pr = Point::Coor(BigUint::from(7u32), BigUint::from(11u32));
        res = ec.scalar_mult(&a, &BigUint::from(10u32)).unwrap();
        assert_eq!(res, pr);

        // 15 * (5, 1) = (3,16)
        pr = Point::Coor(BigUint::from(3u32), BigUint::from(16u32));
        res = ec.scalar_mult(&a, &BigUint::from(15u32)).unwrap();
        assert_eq!(res, pr);

        // 16 * (5, 1) = (10,11)
        pr = Point::Coor(BigUint::from(10u32), BigUint::from(11u32));
        res = ec.scalar_mult(&a, &BigUint::from(16u32)).unwrap();
        assert_eq!(res, pr);

        // 17 * (5, 1) = (6,14)
        pr = Point::Coor(BigUint::from(6u32), BigUint::from(14u32));
        res = ec.scalar_mult(&a, &BigUint::from(17u32)).unwrap();
        assert_eq!(res, pr);

        // 18 * (5, 1) = (5,16)
        pr = Point::Coor(BigUint::from(5u32), BigUint::from(16u32));
        res = ec.scalar_mult(&a, &BigUint::from(18u32)).unwrap();
        assert_eq!(res, pr);

        // 19 * (5, 1) = I
        pr = Point::Identity;
        res = ec.scalar_mult(&a, &BigUint::from(19u32)).unwrap();
        assert_eq!(res, pr);

        // 2 * (10, 6) = (16,13)
        let p1 = Point::Coor(BigUint::from(10u32), BigUint::from(6u32));
        pr = Point::Coor(BigUint::from(16u32), BigUint::from(13u32));
        res = ec.double(&p1);
        assert_eq!(res, pr);
    }

    #[test]
    fn test_ec_secp256k1() {
        /*
          y^2 = x^3 + 7 mod p (large)

          p = FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFE FFFFFC2F 
          p mod 4 == 3
          n = FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFE BAAEDCE6 AF48A03B BFD25E8C D0364141
        G = (
            x = 79BE667E F9DCBBAC 55A06295 CE870B07 029BFCDB 2DCE28D9 59F2815B 16F81798,
            y = 483ADA77 26A3C465 5DA4FBFC 0E1108A8 FD17B448 A6855419 9C47D08F FB10D4B8
        )
        a = 00000000 00000000 00000000 00000000 00000000 00000000 00000000 00000000
        b = 00000000 00000000 00000000 00000000 00000000 00000000 00000000 00000007
        */

        // n * G = I
        let p = BigUint::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .expect("could not convert p");

        let n = BigUint::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
            16,
        )
        .expect("could not convert n");

        let gx = BigUint::parse_bytes(
            b"79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
            16,
        )
        .expect("could not convert gx");

        let gy = BigUint::parse_bytes(
            b"483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
            16,
        )
        .expect("could not convert gy");

        let ec = EllipticCurve {
            a: BigUint::from(0u32),
            b: BigUint::from(7u32),
            p,
        };

        let g = Point::Coor(gx, gy);

        let res = ec.scalar_mult(&g, &n).unwrap();
        assert_eq!(res, Point::Identity);

        // p = 1201 * G -> it is also a generator
        let p = ec.scalar_mult(&g, &BigUint::from(1201u32)).unwrap();

        let res = ec.scalar_mult(&p, &n).unwrap();
        assert_eq!(res, Point::Identity);
    }

    #[test]
    fn test_compress_and_uncompress() {
        let p = BigUint::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .expect("could not convert p");

        let gx = BigUint::parse_bytes(
            b"79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
            16,
        )
        .expect("could not convert gx");

        let gy = BigUint::parse_bytes(
            b"483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
            16,
        )
        .expect("could not convert gy");

        let ec = EllipticCurve {
            a: BigUint::from(0u32),
            b: BigUint::from(7u32),
            p,
        };

        let g = Point::Coor(gx, gy);

        let compressed_point = EllipticCurve::compress(&g).unwrap();
        let uncompressed_g = ec.uncompress(&compressed_point).unwrap();
        assert_eq!(g, uncompressed_g);
    }
}
