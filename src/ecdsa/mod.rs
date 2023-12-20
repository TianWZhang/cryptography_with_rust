use crate::math::{
    elliptic_curve::{EllipticCurve, Point},
    finite_field::FpBigUint,
};
use crate::sha256::Sha256;
use num_bigint::{BigUint, RandBigInt};
use std::ops::Rem;

pub struct ECDSA {
    elliptic_curve: EllipticCurve,
    g: Point,   // group generator
    q: BigUint, // group order
}

#[derive(Debug)]
pub enum ECDSAError {
    BadArgument(String),
    OperationFailure(String),
}

impl ECDSA {
    pub fn new(a: BigUint, b: BigUint, p: BigUint, g: Point, q: BigUint) -> Self {
        let elliptic_curve = EllipticCurve::new(a, b, p);
        Self {
            elliptic_curve,
            g,
            q,
        }
    }

    pub fn generate_private_key(&self) -> BigUint {
        let mut rng = rand::thread_rng();
        rng.gen_biguint_range(&BigUint::from(1u32), &self.q)
    }

    pub fn generate_public_key(&self, private_key: &BigUint) -> Result<Point, ECDSAError> {
        self.elliptic_curve
            .scalar_mult(&self.g, private_key)
            .map_err(|_| ECDSAError::OperationFailure("GenPublicKeyError".into()))
    }

    /// Generates: sk, pk where pk = sk * g
    pub fn generate_key_pair(&self) -> Result<(BigUint, Point), ECDSAError> {
        let sk = self.generate_private_key();
        let pk = self.generate_public_key(&sk)?;
        Ok((sk, pk))
    }

    fn gen_hash_less_than(&self, msg: &str) -> BigUint {
        let mut hasher = Sha256::new();
        hasher.update(msg.as_bytes());
        let digest = hasher.digest();
        let mut digest_bytes = vec![];
        for i in digest.iter() {
            digest_bytes.extend_from_slice(&i.to_be_bytes());
        }
        let hash = BigUint::from_bytes_be(&digest_bytes);
        hash.rem(&self.q)
    }

    /// R = k_e * g -> r = x_R
    /// s = (hash(msg) + sk * r) * k_e^(-1) mod q
    pub fn sign(&self, msg: &str, sk: &BigUint) -> Result<(BigUint, BigUint), ECDSAError> {
        let hash = self.gen_hash_less_than(msg);

        let mut rng = rand::thread_rng();
        let k_e = rng.gen_biguint_range(&BigUint::from(1u32), &self.q); // random ephemeral key
        let r_point = self
            .elliptic_curve
            .scalar_mult(&self.g, &k_e)
            .map_err(|_| {
                // R = k_e * g
                ECDSAError::OperationFailure("Error computing Random Point R".into())
            })?;
        match r_point {
            Point::Identity => Err(ECDSAError::OperationFailure(
                "Random Point R is the identity".into(),
            )),
            Point::Coor(r, _) => {
                // r = x_R
                let fq = FpBigUint::new(self.q.clone());
                let mut s = fq.mult(&r, &sk);
                s = fq.add(&s, &hash);
                let k_e_inv = fq.inv_mult(&k_e);
                s = fq.mult(&s, &k_e_inv); // s = (hash(x) + sk * r) k_e^{-1} mod q
                Ok((r.rem(&self.q), s))
            }
        }
    }

    ///
    /// Verifies if a signature is valid for a particular message and public key.
    ///
    /// (r, s) = signature
    /// u1 = s^(-1) * hash(msg) mod q
    /// u2 = s^(-1) * r mod q
    /// P = u1 * g + u2 * pk mod q = (xp, yp)
    /// if r == xp then verified!
    ///
    pub fn verify(
        &self,
        msg: &str,
        pk: &Point,
        signature: &(BigUint, BigUint),
    ) -> Result<bool, ECDSAError> {
        let hash = self.gen_hash_less_than(msg);
        let (r, s) = signature;
        let fq = FpBigUint::new(self.q.clone());
        let s_inv = fq.inv_mult(s);
        let u1 = fq.mult(&s_inv, &hash);
        let u2 = fq.mult(&s_inv, r);
        let u1g = self
            .elliptic_curve
            .scalar_mult(&self.g, &u1)
            .map_err(|_| ECDSAError::OperationFailure("Error computing u1 * g".into()))?;
        let u2pk = self
            .elliptic_curve
            .scalar_mult(pk, &u2)
            .map_err(|_| ECDSAError::OperationFailure("Error computing u2 * pk".into()))?;
        let p_point = self
            .elliptic_curve
            .add(&u1g, &u2pk)
            .map_err(|_| ECDSAError::OperationFailure("Error computing u1 * g + u2 * pk".into()))?;

        match p_point {
            Point::Identity => Err(ECDSAError::OperationFailure(
                "Point P is the identity".into(),
            )),
            Point::Coor(xp, _) => {
                // xp = x_P
                Ok(xp.rem(&self.q) == *r) // xp == r mod q
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_sign_verify() {
        let a = BigUint::from(2u32);
        let b = BigUint::from(2u32);
        let p = BigUint::from(17u32);
        let g = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        let q = BigUint::from(19u32);
        let ecdsa = ECDSA::new(a, b, p, g, q);

        let sk = BigUint::from(7u32);
        let pk = ecdsa
            .generate_public_key(&sk)
            .expect("Could not compute PubKey");

        let message = "Bob -> 1 BTC -> Alice";
        let signature = ecdsa.sign(message, &sk).expect("Could not sign");

        let verify_result = ecdsa
            .verify(message, &pk, &signature)
            .expect("Could not verify");

        assert!(verify_result, "Verification should success");
    }

    #[test]
    fn test_sign_verify_tempered_message() {
        let a = BigUint::from(2u32);
        let b = BigUint::from(2u32);
        let p = BigUint::from(17u32);
        let g = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        let q = BigUint::from(19u32);
        let ecdsa = ECDSA::new(a, b, p, g, q);

        let sk = BigUint::from(7u32);
        let pk = ecdsa
            .generate_public_key(&sk)
            .expect("Could not compute PubKey");

        let message = "Bob -> 1 BTC -> Alice";
        let signature = ecdsa.sign(message, &sk).expect("Could not sign");

        let another_message = "Bob -> 2 BTC -> Alice";

        let verify_result = ecdsa
            .verify(another_message, &pk, &signature)
            .expect("Could not verify");

        assert!(
            !verify_result,
            "Verification should fail when message is tempered"
        );
    }

    #[test]
    fn test_sign_verify_tempered_signature() {
        let a = BigUint::from(2u32);
        let b = BigUint::from(2u32);
        let p = BigUint::from(17u32);
        let g = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        let q = BigUint::from(19u32);
        let ecdsa = ECDSA::new(a, b, p, g, q);

        let sk = BigUint::from(7u32);
        let pk = ecdsa
            .generate_public_key(&sk)
            .expect("Could not compute PubKey");

        let message = "Bob -> 1 BTC -> Alice";
        let (r, s) = ecdsa.sign(message, &sk).expect("Could not sign");
        let tempered_signature = (
            (r + BigUint::from(1u32)).modpow(&BigUint::from(1u32), &ecdsa.q),
            s,
        );

        let verify_result = ecdsa
            .verify(message, &pk, &tempered_signature)
            .expect("Could not verify");

        assert!(
            !verify_result,
            "Verification should fail when signature is tempered"
        );
    }

    #[test]
    fn test_secp256_sign_verify() {
        let p = BigUint::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .expect("could not convert p");

        let q = BigUint::parse_bytes(
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

        let a = BigUint::from(0u32);
        let b = BigUint::from(7u32);
        let g = Point::Coor(gx, gy);
        let ecdsa = ECDSA::new(a, b, p, g, q);

        let sk = BigUint::parse_bytes(
            b"483ADB7726A3C4655DA4FBFC0E1208A8F017B448A68554199C47D08FFB10E4B9",
            16,
        )
        .expect("Could not convert hex to private key");
        let pk = ecdsa
            .generate_public_key(&sk)
            .expect("Could not compute PubKey");

        let message = "Bob -> 1 BTC -> Alice";
        let signature = ecdsa.sign(message, &sk).expect("Could not sign");

        let verify_result = ecdsa
            .verify(&message, &pk, &signature)
            .expect("Could not verify");

        assert!(verify_result, "Verification should have succeed");
    }

    #[test]
    fn test_secp256_sign_verify_tempered_message() {
        let p = BigUint::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .expect("could not convert p");

        let q = BigUint::parse_bytes(
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

        let a = BigUint::from(0u32);
        let b = BigUint::from(7u32);
        let g = Point::Coor(gx, gy);
        let ecdsa = ECDSA::new(a, b, p, g, q);

        let sk = BigUint::parse_bytes(
            b"483ADB7726A3C4655DA4FBFC0E1208A8F017B448A68554199C47D08FFB10E4B9",
            16,
        )
        .expect("Could not convert hex to private key");
        let pk = ecdsa
            .generate_public_key(&sk)
            .expect("Could not compute PubKey");

        let message = "Bob -> 1 BTC -> Alice";
        let signature = ecdsa.sign(message, &sk).expect("Could not sign");

        let another_message = "Bob -> 2 BTC -> Alice";
        let verify_result = ecdsa
            .verify(another_message, &pk, &signature)
            .expect("Could not verify");

        assert!(
            !verify_result,
            "Verification should have failed due to tempered message"
        );
    }

    #[test]
    fn test_secp256_random_key_pairs() {
        let p = BigUint::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .expect("could not convert p");

        let q = BigUint::parse_bytes(
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

        let a = BigUint::from(0u32);
        let b = BigUint::from(7u32);
        let g = Point::Coor(gx, gy);
        let ecdsa = ECDSA::new(a, b, p, g, q);

        let (sk, pk) = ecdsa.generate_key_pair().unwrap();

        let message = "Bob -> 1 BTC -> Alice";
        let signature = ecdsa.sign(message, &sk).expect("Could not sign");

        let verify_result = ecdsa
            .verify(message, &pk, &signature)
            .expect("Could not verify");

        assert!(
            verify_result,
            "Verification should have failed due to tempered message"
        );
    }
}
