use crate::matrix::MatrixElement;
use crate::PolyGF2;
use num_traits::{One, Zero};
use polynomial::Polynomial;
use std::cmp::Ordering;
use std::collections::HashSet;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign,
};

#[derive(Eq, Hash, Clone, Copy, Debug, PartialEq)]
pub struct GF2TM<const M: u32> {
    value: PolyGF2,
    irr: PolyGF2,
}

impl<const M: u32> Default for GF2TM<M> {
    fn default() -> GF2TM<M> {
        GF2TM {
            value: PolyGF2::default(),
            irr: PolyGF2::irreducible(M),
        }
    }
}

impl<const M: u32> GF2TM<M> {
    pub fn new(value: PolyGF2) -> GF2TM<M> {
        GF2TM {
            value: value % PolyGF2::new(1 << M),
            irr: PolyGF2::irreducible(M),
        }
    }
    pub fn one() -> GF2TM<M> {
        GF2TM {
            value: PolyGF2::new(1),
            irr: PolyGF2::irreducible(M),
        }
    }

    pub fn value(&self) -> PolyGF2 {
        self.value
    }
}

impl<const M: u32> PartialOrd for GF2TM<M> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.value.partial_cmp(&other.value)
    }
}

impl<const M: u32> Add for GF2TM<M> {
    type Output = GF2TM<M>;

    fn add(self, rhs: Self) -> Self::Output {
        GF2TM {
            value: (self.value + rhs.value) % self.irr,
            irr: self.irr,
        }
    }
}

impl<const M: u32> AddAssign for GF2TM<M> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<const M: u32> Sub for GF2TM<M> {
    type Output = GF2TM<M>;

    fn sub(self, rhs: Self) -> Self::Output {
        GF2TM {
            value: (self.value - rhs.value) % self.irr,
            irr: self.irr,
        }
    }
}

impl<const M: u32> SubAssign for GF2TM<M> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<const M: u32> Mul for GF2TM<M> {
    type Output = GF2TM<M>;

    fn mul(self, rhs: Self) -> Self::Output {
        GF2TM {
            value: (self.value * rhs.value) % self.irr,
            irr: self.irr,
        }
    }
}

impl<const M: u32> MulAssign for GF2TM<M> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<const M: u32> Div for GF2TM<M> {
    type Output = GF2TM<M>;

    fn div(self, rhs: Self) -> Self::Output {
        GF2TM {
            value: (self.value * rhs.inv().value) % self.irr,
            irr: self.irr,
        }
    }
}

impl<const M: u32> DivAssign for GF2TM<M> {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl<const M: u32> Rem for GF2TM<M> {
    type Output = GF2TM<M>;

    fn rem(self, rhs: Self) -> Self::Output {
        GF2TM {
            value: (self.value % rhs.value) % self.irr,
            irr: self.irr,
        }
    }
}

impl<const M: u32> RemAssign for GF2TM<M> {
    fn rem_assign(&mut self, rhs: Self) {
        *self = *self % rhs;
    }
}

impl<const M: u32> Neg for GF2TM<M> {
    type Output = GF2TM<M>;

    fn neg(self) -> Self::Output {
        GF2TM {
            value: -self.value,
            irr: self.irr,
        }
    }
}

impl<const M: u32> From<u32> for GF2TM<M> {
    fn from(poly: u32) -> Self {
        GF2TM {
            value: PolyGF2::new(poly),
            irr: PolyGF2::irreducible(M),
        }
    }
}
impl<const M: u32> From<u8> for GF2TM<M> {
    fn from(poly: u8) -> Self {
        Self::from(poly as u32)
    }
}

impl<const M: u32> Zero for GF2TM<M> {
    fn zero() -> Self {
        GF2TM {
            value: PolyGF2::default(),
            irr: PolyGF2::irreducible(M),
        }
    }

    fn is_zero(&self) -> bool {
        self.value.is_zero()
    }
}
impl<const M: u32> One for GF2TM<M> {
    fn one() -> Self {
        GF2TM {
            value: PolyGF2::new(1),
            irr: PolyGF2::irreducible(M),
        }
    }

    fn is_one(&self) -> bool {
        self.value == PolyGF2::new(1)
    }
}

impl<const M: u32> GF2TM<M> {
    pub fn pow(&self, exp: u32) -> GF2TM<M> {
        let mut result = GF2TM::new(PolyGF2::new(1 as u32));
        for _ in 0..exp {
            result *= *self;
        }
        result
    }

    pub fn inv(&self) -> Self {
        self.pow((1 << M) - 2)
    }

    pub fn minimal_poly(&self) -> PolyGF2 {
        let mut conjugates = HashSet::new();
        for i in 0..M {
            let new = self.pow(2_u32.pow(i));
            conjugates.insert(new);
        }
        conjugates
            .iter()
            .map(|a| Polynomial::new(vec![*a, GF2TM::one()]))
            .fold(Polynomial::new(vec![GF2TM::one()]), |acc, e| acc * e)
            .into()
    }

    pub fn is_primitive(&self) -> bool {
        let order = (1 << M) - 1;
        let mut powers = PolyGF2::new(1);

        for _ in 1..order {
            powers = (powers * self.value) % self.irr;
            if powers == PolyGF2::new(1) {
                return false;
            }
        }

        powers = (powers * self.value) % self.irr;
        powers == PolyGF2::new(1)
    }

    pub fn primitive_element() -> GF2TM<M> {
        for candidate in 1..(1 << M) {
            let alpha = GF2TM::new(PolyGF2::new(candidate));
            if alpha.is_primitive() {
                return alpha;
            }
        }
        unreachable!("Invalid irreducible polynomial, use M > 1");
    }
}

impl<const M: u32> MatrixElement for GF2TM<M> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_minimal_poly() {
        assert_eq!(
            GF2TM::<2>::from(0b10u32).minimal_poly(),
            PolyGF2::new(0b111)
        );
        assert_eq!(
            GF2TM::<2>::from(0b11u32).minimal_poly(),
            PolyGF2::new(0b111)
        );

        assert_eq!(PolyGF2::irreducible(3), PolyGF2::new(0b1011));

        let elem = GF2TM::<3>::from(0b1u32);
        assert_eq!(elem.minimal_poly(), PolyGF2::new(0b11));

        let elem = GF2TM::<3>::from(0b11u32);
        assert_eq!(elem.minimal_poly(), PolyGF2::new(0b1101));

        let elem = GF2TM::<3>::from(0b10u32);
        assert_eq!(elem.minimal_poly(), PolyGF2::new(0b1011));

        let elem = GF2TM::<4>::from(2u32);
        assert_eq!(elem.minimal_poly(), PolyGF2::new(0b10011));

        let elem = GF2TM::<4>::from(3u32);
        assert_eq!(elem.minimal_poly(), PolyGF2::new(0b10011));

        let elem = GF2TM::<4>::from(6u32);
        assert_eq!(elem.minimal_poly(), PolyGF2::new(0b111));

        let elem = GF2TM::<4>::from(12u32);
        assert_eq!(elem.minimal_poly(), PolyGF2::new(0b11111));
    }

    #[test]
    fn test_primitive() {
        assert_eq!(GF2TM::<2>::primitive_element(), GF2TM::from(0b10u32));
        assert_eq!(GF2TM::<3>::primitive_element(), GF2TM::from(0b10u32));
        assert!(GF2TM::<2>::from(0b11u32).is_primitive());
        assert!(GF2TM::<3>::from(0b11u32).is_primitive());
    }
}
