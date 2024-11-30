use crate::GF2TM;
use num_traits::Zero;
use polynomial::Polynomial;
use std::fmt::Debug;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign,
};

#[derive(Eq, PartialEq, Hash, Clone, Copy, Default, PartialOrd)]
pub struct PolyGF2 {
    pub poly: u32,
}

impl Zero for PolyGF2 {
    fn zero() -> Self {
        PolyGF2 { poly: 0 }
    }

    fn is_zero(&self) -> bool {
        self.poly == 0
    }
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
    pub fn divmod(self, rhs: Self) -> (Self, Self) {
        if rhs.is_zero() {
            panic!("division by zero");
        }

        let mut quotient = 0;
        let mut remainder = self.poly;
        let divisor = rhs.poly;

        while remainder != 0 && remainder.leading_zeros() <= divisor.leading_zeros() {
            // Compute the shift needed to align leading terms
            let shift = divisor.leading_zeros() - remainder.leading_zeros();
            // Update the quotient by adding the shifted term
            quotient ^= 1 << shift;
            // Update the remainder by subtracting (XOR) the shifted divisor
            remainder ^= divisor << shift;
        }

        (
            PolyGF2::new(quotient),  // The quotient
            PolyGF2::new(remainder), // The remainder
        )
    }

    pub fn pow(&self, rhs: u32) -> Self {
        let mut result = PolyGF2::new(1);
        for _ in 0..rhs {
            result *= *self;
        }
        result
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

impl Neg for PolyGF2 {
    type Output = PolyGF2;

    fn neg(self) -> Self::Output {
        PolyGF2::default() - self
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

impl<const M: u32> From<Polynomial<GF2TM<M>>> for PolyGF2 {
    fn from(poly: Polynomial<GF2TM<M>>) -> Self {
        poly.data()
            .iter()
            .rev()
            .fold(0, |acc, x| acc << 1 | x.value().poly & 1)
            .into()
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

    pub fn degree(&self) -> usize {
        (self.poly.leading_zeros() ^ 31) as usize
    }

    pub fn gcd(&self, rhs: Self) -> Self {
        let mut a = *self;
        let mut b = rhs;
        while !b.is_zero() {
            let r = a % b;
            a = b;
            b = r;
        }
        a
    }

    pub fn lcm(&self, rhs: Self) -> Self {
        let gcd = self.gcd(rhs);
        *self * rhs / gcd
    }

    pub fn eval(&self, x: u32) -> u32 {
        let x = x & 1;
        let mut poly = self.poly;
        let mut result = poly & 1;
        poly >>= 1;
        while poly > 0 {
            result ^= poly & x;
            poly >>= 1;
        }
        result
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

    #[test]
    fn test_reduction() {
        assert_eq!(PolyGF2::new(0b110) % PolyGF2::new(0b111), PolyGF2::new(0b1),);
    }

    #[test]
    fn test_degree() {
        assert_eq!(PolyGF2::new(0b110).degree(), 2);
        assert_eq!(PolyGF2::new(0b111).degree(), 2);
        assert_eq!(PolyGF2::new(0b1011).degree(), 3);
        assert_eq!(PolyGF2::new(0b10011).degree(), 4);
        assert_eq!(PolyGF2::new(0b100101).degree(), 5);
    }

    #[test]
    fn test_gcd() {
        let a = PolyGF2::new(0b1001); // (x+1)(x^2+x+1)
        let b = PolyGF2::new(0b11101); // (x+1)(x^3+x+1)
        assert_eq!(a.gcd(b), PolyGF2::new(0b11)); // (x+1)

        let a = PolyGF2::new(0b11011); // (x+1)(x+1)(x^2+x+1)
        let b = PolyGF2::new(0b100111); // (x+1)(x+1)(x^3+x+1)
        assert_eq!(a.gcd(b), PolyGF2::new(0b101)); // (x+1)(x+1)
    }

    #[test]
    fn test_lcm() {
        let a = PolyGF2::new(0b1001); // (x+1)(x^2+x+1)
        let b = PolyGF2::new(0b11101); // (x+1)(x^3+x+1)
        assert_eq!(a.lcm(b), PolyGF2::new(0b1010011)); // (x+1)(x^2+x+1)(x^3+x+1)
        let a = PolyGF2::new(0b11011); // (x+1)(x+1)(x^2+x+1)
        let b = PolyGF2::new(0b100111); // (x+1)(x+1)(x^3+x+1)
        assert_eq!(a.lcm(b), PolyGF2::new(0b11110101)); // (x+1)(x+1)(x^2+x+1)(x^3+x+1)
    }

    #[test]
    fn test_eval() {
        let a = PolyGF2::new(0b110);
        assert_eq!(a.eval(0), 0);
        assert_eq!(a.eval(1), 0);
        let a = PolyGF2::new(0b101);
        assert_eq!(a.eval(0), 1);
        assert_eq!(a.eval(1), 0);
    }

    #[test]
    fn test_from_poly_over_gf2m() {
        let poly = Polynomial::new(vec![
            GF2TM::<3>::new(PolyGF2::new(1)),
            GF2TM::<3>::new(PolyGF2::new(0)),
            GF2TM::<3>::new(PolyGF2::new(1)),
            GF2TM::<3>::new(PolyGF2::new(1)),
        ]);
        assert_eq!(PolyGF2::from(poly), PolyGF2::new(0b1101));
    }
}
