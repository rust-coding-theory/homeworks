use crate::PolyGF2;
use num_traits::Zero;
use polynomial::Polynomial;
use std::collections::HashSet;
use std::fmt::Display;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign,
};

#[derive(Eq, PartialEq, Hash, Clone, Copy, Default, Debug)]
pub struct GF2TM<const M: u32> {
    value: PolyGF2,
    irr: PolyGF2,
}

impl<const M: u32> GF2TM<M> {
    pub fn new(value: PolyGF2) -> GF2TM<M> {
        GF2TM {
            value: value % PolyGF2::new(1 << M),
            irr: PolyGF2::irreducible(M),
        }
    }

    pub fn default() -> GF2TM<M> {
        GF2TM {
            value: PolyGF2::default(),
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
            value: (self.value / rhs.value) % self.irr,
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

impl<const M: u32> GF2TM<M> {
    pub fn pow(&self, exp: u32) -> GF2TM<M> {
        GF2TM {
            value: self.value.pow(exp) % self.irr,
            irr: self.irr,
        }
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
            .reduce(|acc, e| acc * e)
            .expect("Empty set of conjugates")
            .data()
            .iter()
            .map(|x| (x.value % self.irr).poly)
            .reduce(|acc, e| acc << 1 ^ (e & 1))
            .expect("Empty polynomial")
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
        unreachable!("Invalid irreducible polynomial");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_minimal_poly() {
        assert_eq!(GF2TM::<2>::from(0b10).minimal_poly(), PolyGF2::new(0b111));
        assert_eq!(GF2TM::<2>::from(0b11).minimal_poly(), PolyGF2::new(0b111));

        const DEGREE: u32 = 3;
        assert_eq!(PolyGF2::irreducible(DEGREE), PolyGF2::new(0b1011));

        let elem = GF2TM::<DEGREE>::from(0b1);
        assert_eq!(elem.minimal_poly(), PolyGF2::new(0b11));

        let elem = GF2TM::<DEGREE>::from(0b11);
        assert_eq!(elem.minimal_poly(), PolyGF2::new(0b1011));

        let elem = GF2TM::<DEGREE>::from(0b10);
        assert_eq!(elem.minimal_poly(), PolyGF2::new(0b1101));
    }

    #[test]
    fn test_primitive() {
        assert_eq!(GF2TM::<2>::primitive_element(), GF2TM::from(0b10));
        assert_eq!(GF2TM::<3>::primitive_element(), GF2TM::from(0b10));
        assert!(GF2TM::<2>::from(0b11).is_primitive());
        assert!(GF2TM::<3>::from(0b11).is_primitive());
    }
}
