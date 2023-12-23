use super::finite_field::FiniteRing;
use std::ops::{Add, AddAssign, Mul, Neg, Sub};

/// Polynomial in the ring T[X]/(X^N + 1)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Polynomial<T, const N: usize>
where
    T: Clone + PartialEq + Eq,
{
    pub coefficients: [T; N],
    /// The zero polynomial has degree < 0
    pub degree: Option<usize>,
}

impl<T: FiniteRing, const N: usize> Add for Polynomial<T, N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        if self.is_zero() {
            return rhs.clone();
        }
        if rhs.is_zero() {
            return self.clone();
        }

        let mut degree = self.degree.unwrap().max(rhs.degree.unwrap());
        let mut coefficients = [T::ZERO; N];
        for i in 0..N {
            coefficients[i] = self.coefficients[i] + rhs.coefficients[i];
        }

        // Diminish degree if the leading coefficient is zero
        let mut leading = coefficients[degree];
        while degree > 0 && leading.eq(&T::ZERO) {
            degree -= 1;
            leading = coefficients[degree];
        }

        if degree == 0 && leading == T::ZERO {
            return Self::ZERO;
        }

        Self {
            coefficients,
            degree: Some(degree),
        }
    }
}

impl<T: FiniteRing, const N: usize> Neg for Polynomial<T, N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.degree.is_none() {
            return Self::ZERO;
        }

        let mut coefficients = [T::ZERO; N];
        for (i, c) in self.coefficients.iter().enumerate() {
            coefficients[i] = -(*c);
        }

        Self {
            coefficients,
            degree: self.degree,
        }
    }
}

impl<T: FiniteRing, const N: usize> Sub for Polynomial<T, N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<T: FiniteRing, const N: usize> AddAssign for Polynomial<T, N> {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone() + rhs;
    }
}

impl<T: FiniteRing, const N: usize> Mul for Polynomial<T, N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.is_zero() || rhs.is_zero() {
            return Self::ZERO;
        }
        let mut coefficients = [T::ZERO; N];

        for i in 0..N {
            for j in 0..N {
                let c = self.coefficients[i] * rhs.coefficients[j];
                let k = i + j;
                if k < N {
                    coefficients[k] += c;
                } else {
                    // X^n + 1 = 0
                    coefficients[k % N] += -c;
                }
            }
        }

        let mut degree = N - 1;
        while degree > 0 && coefficients[degree] == T::ZERO {
            degree -= 1;
        }

        if degree == 0 && coefficients[0] == T::ZERO {
            return Self::ZERO;
        }

        Self {
            coefficients,
            degree: Some(degree),
        }
    }
}

impl<T: FiniteRing, const N: usize> FiniteRing for Polynomial<T, N> {
    const ZERO: Self = Self {
        coefficients: [T::ZERO; N],
        degree: None,
    };

    const ONE: Self = {
        let mut coefficients = [T::ZERO; N];
        coefficients[0] = T::ONE;
        Self {
            coefficients,
            degree: Some(0),
        }
    };
}

impl<T: FiniteRing, const N: usize> Polynomial<T, N> {
    pub fn is_zero(&self) -> bool {
        self.degree.is_none()
    }

    pub fn ring_dimension() -> usize {
        N
    }

    pub fn with_coefficients(coefficients: [T; N]) -> Self {
        let mut degree = N - 1;
        while degree > 0 && coefficients[degree] == T::ZERO {
            degree -= 1;
        }

        if degree == 0 && coefficients[0] == T::ZERO {
            return Self {
                coefficients,
                degree: None,
            };
        }

        Self {
            coefficients,
            degree: Some(degree),
        }
    }

    pub fn scalar_mult(&self, scalar: &T) -> Self {
        if self.is_zero() || scalar == &T::ZERO {
            return Self {
                coefficients: [T::ZERO; N],
                degree: None,
            };
        }

        let degree = self.degree.unwrap();
        let mut coefficients = [T::ZERO; N];
        for i in 0..degree {
            coefficients[i] = self.coefficients[i] * (*scalar);
        }
        Self::with_coefficients(coefficients)
    }

    pub fn set_coefficient(&mut self, index: usize, val: T) {
        if index < N {
            if val != T::ZERO {
                self.degree = match self.degree {
                    Some(d) if d < index => Some(index),
                    Some(d) => Some(d),
                    None => Some(index),
                }
            } else {
                self.degree = match self.degree {
                    Some(d) if d == index && d == 0 => None,
                    Some(d) if d == index => {
                        // now d >= 1
                        let mut degree = d - 1;
                        while degree > 0 && self.coefficients[degree] == T::ZERO {
                            degree -= 1;
                        }
                        if degree == 0 && self.coefficients[0] == T::ZERO {
                            None
                        } else {
                            Some(degree)
                        }
                    }
                    Some(d) => Some(d),
                    None => None,
                }
            }
            self.coefficients[index] = val;
        }
    }
}

impl<T: FiniteRing, const N: usize> Default for Polynomial<T, N> {
    fn default() -> Self {
        Self {
            coefficients: [T::ZERO; N],
            degree: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{Poly3329, F3329};

    #[test]
    fn test_scalar_mult() {
        let mut poly = Poly3329::<3>::with_coefficients([6.into(), 3324.into(), 3.into()]);
        poly.set_coefficient(2, F3329::ZERO);
        assert_eq!(poly.degree.unwrap(), 1);

        poly.set_coefficient(2, F3329::from(5));
        assert_eq!(poly.degree.unwrap(), 2);

        poly.set_coefficient(1, F3329::ZERO);
        assert_eq!(poly.degree.unwrap(), 2);

        poly.set_coefficient(2, F3329::ZERO);
        assert_eq!(poly.degree.unwrap(), 0);

        poly.set_coefficient(0, F3329::ZERO);
        assert!(poly.degree.is_none());
    }
}
