use std::fmt::Debug;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Rem, RemAssign, Sub, SubAssign};

#[derive(Eq, PartialEq, Hash, Clone, Copy, Default)]
pub struct PolyGF2 {
    pub poly: u32,
}

impl Debug for PolyGF2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:b}", self.poly)
    }
}

impl Add for PolyGF2 {
    type Output = PolyGF2;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn add(self, rhs: Self) -> Self::Output {
        PolyGF2 {
            poly: self.poly ^ rhs.poly,
        }
    }
}
impl AddAssign for PolyGF2 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for PolyGF2 {
    type Output = PolyGF2;

    fn sub(self, rhs: Self) -> Self::Output {
        self.add(rhs)
    }
}

impl SubAssign for PolyGF2 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}

impl Mul for PolyGF2 {
    type Output = PolyGF2;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = 0;
        let mut a = self.poly;
        let mut b = rhs.poly;
        while b > 0 {
            if b & 1 > 0 {
                result ^= a;
            }
            a <<= 1;
            b >>= 1;
        }
        PolyGF2 { poly: result }
    }
}

impl MulAssign for PolyGF2 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl PolyGF2 {
    fn divmod(self, rhs: Self) -> (Self, Self) {
        if rhs == PolyGF2::default() {
            panic!("division by zero");
        }
        if self == PolyGF2::default() {
            return (PolyGF2::default(), PolyGF2::default());
        }
        if self.poly < rhs.poly {
            return (PolyGF2::default(), self);
        }

        let mut result = 0;
        let mut a = self.poly;
        let b = rhs.poly;
        while a >= b {
            let shift = b.leading_zeros().saturating_sub(a.leading_zeros());
            result ^= 1 << shift;
            a ^= b << shift;
        }
        (PolyGF2 { poly: result }, PolyGF2 { poly: a })
    }
}

impl Div for PolyGF2 {
    type Output = PolyGF2;

    fn div(self, rhs: Self) -> Self::Output {
        self.divmod(rhs).0
    }
}

impl DivAssign for PolyGF2 {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl Rem for PolyGF2 {
    type Output = PolyGF2;

    fn rem(self, rhs: Self) -> Self::Output {
        self.divmod(rhs).1
    }
}

impl RemAssign for PolyGF2 {
    fn rem_assign(&mut self, rhs: Self) {
        *self = *self % rhs;
    }
}

impl From<u32> for PolyGF2 {
    fn from(poly: u32) -> Self {
        PolyGF2 { poly }
    }
}

impl From<PolyGF2> for u32 {
    fn from(value: PolyGF2) -> Self {
        value.poly
    }
}

impl PolyGF2 {
    pub fn new(poly: u32) -> Self {
        PolyGF2 { poly }
    }

    pub fn irreducible(degree: u32) -> Self {
        for poly in 1 << degree..=1 << (degree + 1) {
            let mut is_irreducible = true;
            for i in 2..1 << degree {
                if (PolyGF2::new(poly) % PolyGF2::new(i)) == PolyGF2::default() {
                    is_irreducible = false;
                    break;
                }
            }
            if is_irreducible {
                return PolyGF2::new(poly);
            }
        }
        PolyGF2::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poly_add() {
        let a = PolyGF2::new(0b101);
        let b = PolyGF2::new(0b110);
        assert_eq!(a + b, PolyGF2::new(0b011));
        let b = PolyGF2::new(0b01);
        let a = PolyGF2::new(0b10);
        assert_eq!(a + b, PolyGF2::new(0b11));
    }

    #[test]
    fn test_poly_mul() {
        let a = PolyGF2::new(0b101);
        let b = PolyGF2::new(0b110);
        assert_eq!(a * b, PolyGF2::new(0b11110));
        let a = PolyGF2::new(0b101);
        let b = PolyGF2::new(0b111);
        assert_eq!(a * b, PolyGF2::new(0b11011));
    }

    #[test]
    fn test_poly_mul_fuzzy() {
        for mut a in 1..100 {
            for mut b in 1..100 {
                let result_poly = PolyGF2::new(a) * PolyGF2::new(b);
                let mut result = 0;
                while b != 0 {
                    if b & 1 > 0 {
                        result ^= a;
                    }
                    a <<= 1;
                    b >>= 1;
                }
                assert_eq!(result_poly.poly, result);
            }
        }
    }

    #[test]
    fn test_poly_divmod_fuzzy() {
        for a in 0..100 {
            for b in 1..100 {
                let (q, r) = PolyGF2::new(a).divmod(PolyGF2::new(b));
                assert_eq!(q * PolyGF2::new(b) + r, PolyGF2::new(a));
                assert!(r.poly.leading_zeros() >= b.leading_zeros())
            }
        }
    }

    #[test]
    fn test_poly_irreducible() {
        let degree = 1;
        assert_eq!(PolyGF2::irreducible(degree), PolyGF2::new(0b10));
        let degree = 2;
        assert_eq!(PolyGF2::irreducible(degree), PolyGF2::new(0b111));
        let degree = 3;
        assert_eq!(PolyGF2::irreducible(degree), PolyGF2::new(0b1011));
        let degree = 4;
        assert_eq!(PolyGF2::irreducible(degree), PolyGF2::new(0b10011));
        let degree = 5;
        assert_eq!(PolyGF2::irreducible(degree), PolyGF2::new(0b100101));
        let degree = 6;
        assert_eq!(PolyGF2::irreducible(degree), PolyGF2::new(0b1000011));
        let degree = 7;
        assert_eq!(PolyGF2::irreducible(degree), PolyGF2::new(0b10000011));
    }
}
