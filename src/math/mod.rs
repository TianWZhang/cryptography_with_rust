use self::{
    finite_field::Fp,
    matrix::{Matrix, Vector},
    poly::Polynomial,
};

pub mod elliptic_curve;
pub mod finite_field;
pub mod matrix;
pub mod ntt;
pub mod poly;

/// Polynomial Ring Rq = Zq[X]/(X^N+1), q = 3329
pub type F3329 = Fp<3329>;
pub type Poly3329<const N: usize> = Polynomial<F3329, N>;
pub type PolyVec3329<const N: usize, const D: usize> = Vector<Poly3329<N>, D>;
pub type PolyMatrix3329<const N: usize, const X: usize, const Y: usize> = Matrix<Poly3329<N>, X, Y>;
