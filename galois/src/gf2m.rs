use crate::poly_gf2::PolyGF2;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Rem, RemAssign, Sub, SubAssign};

#[derive(Eq, PartialEq, Hash, Clone, Copy, Default, Debug)]
pub struct GF2m {
    value: PolyGF2,
    m: u32,
    irr: PolyGF2,
}

impl GF2m {
    pub fn new(value: PolyGF2, m: u32, irr: PolyGF2) -> GF2m {
        GF2m {
            value: value % PolyGF2::new(1 << m),
            m,
            irr,
        }
    }
}

impl Add for GF2m {
    type Output = GF2m;

    fn add(self, rhs: Self) -> Self::Output {
        GF2m {
            value: self.value + rhs.value,
            m: self.m,
            irr: self.irr,
        }
    }
}

impl AddAssign for GF2m {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for GF2m {
    type Output = GF2m;

    fn sub(self, rhs: Self) -> Self::Output {
        GF2m {
            value: self.value - rhs.value,
            m: self.m,
            irr: self.irr,
        }
    }
}

impl SubAssign for GF2m {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul for GF2m {
    type Output = GF2m;

    fn mul(self, rhs: Self) -> Self::Output {
        GF2m {
            value: self.value * rhs.value,
            m: self.m,
            irr: self.irr,
        }
    }
}

impl MulAssign for GF2m {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Div for GF2m {
    type Output = GF2m;

    fn div(self, rhs: Self) -> Self::Output {
        GF2m {
            value: self.value / rhs.value,
            m: self.m,
            irr: self.irr,
        }
    }
}

impl DivAssign for GF2m {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl Rem for GF2m {
    type Output = GF2m;

    fn rem(self, rhs: Self) -> Self::Output {
        GF2m {
            value: self.value % rhs.value,
            m: self.m,
            irr: self.irr,
        }
    }
}

impl RemAssign for GF2m {
    fn rem_assign(&mut self, rhs: Self) {
        *self = *self % rhs;
    }
}

impl GF2m {
    pub fn pow(self, mut exp: u32) -> GF2m {
        let mut result = GF2m::new(PolyGF2::new(1), self.m, self.irr);
        let mut base = self;
        while exp > 0 {
            if exp & 1 == 1 {
                result *= base;
            }
            base *= base;
            exp >>= 1;
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        let a = GF2m::new(PolyGF2::new(0b01), 2, PolyGF2::new(0b111));
        let b = GF2m::new(PolyGF2::new(0b10), 2, PolyGF2::new(0b111));
        let c = a + b;
        assert_eq!(c.value, PolyGF2::new(0b11));
    }

    #[test]
    fn test_sub() {
        let a = GF2m::new(PolyGF2::new(0b01), 2, PolyGF2::new(0b111));
        let b = GF2m::new(PolyGF2::new(0b10), 2, PolyGF2::new(0b111));
        let c = a - b;
        assert_eq!(c.value, PolyGF2::new(0b11));
    }

    #[test]
    fn test_mul() {
        let a = GF2m::new(PolyGF2::new(0b01), 2, PolyGF2::new(0b111));
        let b = GF2m::new(PolyGF2::new(0b10), 2, PolyGF2::new(0b111));
        let c = a * b;
        assert_eq!(c.value, PolyGF2::new(0b10));
    }

    #[test]
    fn test_div() {
        let a = GF2m::new(PolyGF2::new(0b01), 2, PolyGF2::new(0b111));
        let b = GF2m::new(PolyGF2::new(0b10), 2, PolyGF2::new(0b111));
        let c = a / b;
        assert_eq!(c.value, PolyGF2::new(0b0));
        let a = GF2m::new(PolyGF2::new(0b10), 2, PolyGF2::new(0b111));
        let b = GF2m::new(PolyGF2::new(0b01), 2, PolyGF2::new(0b111));
        let c = a / b;
        assert_eq!(c.value, PolyGF2::new(0b10));
    }

    #[test]
    fn test_rem() {
        let a = GF2m::new(PolyGF2::new(0b01), 2, PolyGF2::new(0b111));
        let b = GF2m::new(PolyGF2::new(0b10), 2, PolyGF2::new(0b111));
        let c = a % b;
        assert_eq!(c.value, PolyGF2::new(0b1));
    }

    #[test]
    fn test_pow() {
        let a = GF2m::new(PolyGF2::new(0b01), 2, PolyGF2::new(0b111));
        let b = a.pow(3);
        assert_eq!(b.value, PolyGF2::new(0b01));
    }
}
