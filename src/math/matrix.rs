use std::{ops::{Neg, Sub, Add, Mul, Index, IndexMut}, process::Output};

use super::finite_field::FiniteRing;

pub struct Matrix<T, const X: usize, const Y: usize> 
where 
    T: FiniteRing + Clone + Default
{
    entries: [[T; X]; Y]
}

#[derive(Clone, Copy)]
pub struct Vector<T, const D: usize> 
where 
    T: FiniteRing + Clone + Default + Copy
{
    entries: [T; D]
}

impl<T, const D: usize> Neg for Vector<T, D>
where 
    T: FiniteRing + Clone + Default + Copy + Neg<Output = T>
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut entries = [T::default(); D];
        for i in 0..D {
            entries[i] = -self.entries[i]; 
        }
        Self {
            entries
        }
    }
}

impl<T, const D: usize> Add for Vector<T, D>
where 
    T: FiniteRing + Clone + Default + Copy + Add<Output = T>
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut entries = [T::default(); D];
        for i in 0..D {
            entries[i] = self.entries[i] + rhs.entries[i]; 
        }
        Self {
            entries
        }
    }
}

impl<T, const D: usize> Sub for Vector<T, D>
where 
    T: FiniteRing + Clone + Default + Copy + Add<Output = T> + Neg<Output = T>
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<T, const D: usize> Vector<T, D>
where 
    T: FiniteRing + Clone +  Default + Copy + Mul<Output = T>
{
    pub fn zero() -> Self {
        Self {
            entries: [T::default(); D]
        }
    }

    pub fn dimension() -> usize {
        D
    }

    pub fn dot(&self, other: &Self) -> T {
        let mut res = T::zero();

        for i in 0..D {
            res += self.entries[i] * other.entries[i];
        }
        res
    }

    pub fn basis_vector(index: usize) -> Self {
        let mut res = Self::zero();
        res[index] = T::one();
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

impl<T, const D: usize> Index<usize> for Vector<T, D>
where 
    T: FiniteRing + Clone +  Default + Copy
{
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.entries[index]
    }
}

impl<T, const D: usize> IndexMut<usize> for Vector<T, D>
where 
    T: FiniteRing + Clone +  Default + Copy
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.entries[index]
    }
}


impl<K, const X: usize, const Y: usize> Matrix<K, X, Y>
where 
    K: FiniteRing + Clone +  Default + Copy
{
    pub fn new() -> Self {
        Self {
            entries: [[K::default(); X]; Y]
        }
    }

    pub fn dimensions() -> (usize, usize) {
        (X, Y)
    }
}

