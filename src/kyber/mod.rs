use self::{kem::KEM, pke::PKE};

mod compress;
mod encode;
pub mod kem;
pub mod pke;
mod utils;

pub const KYBER512PKE: PKE<2> = PKE::new(3, 2, 10, 4);
pub const KYBER768PKE: PKE<3> = PKE::new(2, 2, 10, 4);
pub const KYBER1024PKE: PKE<4> = PKE::new(2, 2, 11, 5);
pub const KYBER512KEM: KEM<2> = KEM::new(KYBER512PKE, 1632);
pub const KYBER768KEM: KEM<3> = KEM::new(KYBER768PKE, 2400);
pub const KYBER1024KEM: KEM<4> = KEM::new(KYBER1024PKE, 3168);
