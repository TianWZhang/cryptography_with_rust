use super::finite_field::FiniteRing;
use std::ops::{Add, Index, IndexMut, Neg, Sub};

pub struct Matrix<T, const M: usize, const N: usize> {
    entries: [[T; N]; M],
}

#[derive(Clone, Copy)]
pub struct Vector<T, const D: usize> {
    entries: [T; D],
}

impl<T, const D: usize> Neg for Vector<T, D>
where
    T: FiniteRing + Clone + Default + Copy,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut entries = [T::default(); D];
        for i in 0..D {
            entries[i] = -self.entries[i];
        }
        Self { entries }
    }
}

impl<T, const D: usize> Add for Vector<T, D>
where
    T: FiniteRing + Clone + Default + Copy,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut entries = [T::default(); D];
        for i in 0..D {
            entries[i] = self.entries[i] + rhs.entries[i];
        }
        Self { entries }
    }
}

impl<T, const D: usize> Sub for Vector<T, D>
where
    T: FiniteRing + Clone + Default + Copy,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<T, const D: usize> Vector<T, D>
where
    T: FiniteRing + Clone + Copy + Default,
{
    pub fn zero() -> Self {
        Self {
            entries: [T::ZERO; D],
        }
    }

    pub fn dimension() -> usize {
        D
    }

    pub fn dot(&self, other: &Self) -> T {
        let mut res = T::ZERO;

        for i in 0..D {
            res += self.entries[i] * other.entries[i];
        }
        res
    }

    pub fn basis_vector(index: usize) -> Self {
        let mut res = Self::zero();
        res[index] = T::ONE;
        res
    }

    pub fn scalar_mult(&self, scalar: &T) -> Self {
        let mut res = Self::zero();
        for i in 0..D {
            res.entries[i] = self.entries[i] * scalar.clone();
        }
        res
    }
}

impl<T, const D: usize> Index<usize> for Vector<T, D> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.entries[index]
    }
}

impl<T, const D: usize> IndexMut<usize> for Vector<T, D> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.entries[index]
    }
}

impl<T, const D: usize> From<[T; D]> for Vector<T, D> {
    fn from(value: [T; D]) -> Self {
        Vector { entries: value }
    }
}

impl<T, const M: usize, const N: usize> Matrix<T, M, N>
where
    T: FiniteRing + Clone + Default + Copy,
{
    pub fn zero() -> Self {
        Self {
            entries: [[T::ZERO; N]; M],
        }
    }

    pub fn dimensions() -> (usize, usize) {
        (M, N)
    }

    pub fn row(&self, index: usize) -> Vector<T, N> {
        Vector {
            entries: self.entries[index],
        }
    }

    pub fn column(&self, index: usize) -> Vector<T, M> {
        let mut entries = [T::ZERO; M];

        for i in 0..M {
            entries[i] = self.entries[i][index].clone();
        }

        Vector { entries }
    }

    pub fn mult_vec(&self, v: Vector<T, N>) -> Vector<T, M> {
        let mut entries = [T::ZERO; M];

        for i in 0..M {
            entries[i] = v.dot(&self.row(i));
        }
        Vector { entries }
    }
}

impl<T, const M: usize, const N: usize> Index<(usize, usize)> for Matrix<T, M, N> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (i, j) = index;
        &self.entries[i][j]
    }
}

impl<T, const M: usize, const N: usize> IndexMut<(usize, usize)> for Matrix<T, M, N> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (i, j) = index;
        &mut self.entries[i][j]
    }
}
