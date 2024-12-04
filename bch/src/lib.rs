use galois::{PolyGF2, GF2TM};
use polynomial::Polynomial;

use galois::Matrix;
use num_traits::Zero;

#[derive(Debug, Clone, Copy)]
pub struct BCH<const M: u32> {
    primitive_element: GF2TM<M>,
    distance: usize,
    code_length: usize,
    message_length: usize,
    generator_poly: PolyGF2,
}

impl<const M: u32> BCH<M> {
    pub fn from_distance(distance: usize) -> Self {
        Self::from_max_errors((distance - 1) / 2)
    }

    pub fn from_max_errors(max_errors: usize) -> Self {
        let distance = 2 * max_errors + 1;
        let primitive_element = GF2TM::<M>::primitive_element();
        let code_length = 2_usize.pow(M) - 1;

        let generator_poly = (1..distance)
            .map(|i| primitive_element.pow(i as u32).minimal_poly())
            .reduce(|acc, e| acc.lcm(e))
            .unwrap();
        let message_length = code_length - generator_poly.degree();
        BCH {
            primitive_element,
            distance,
            code_length,
            message_length,
            generator_poly,
        }
    }

    pub fn encode(&self, message: PolyGF2) -> Result<PolyGF2, &'static str> {
        let message_length = message.degree() + 1;
        if message_length > self.message_length {
            return Err("Message is too long");
        }
        let padded = message * PolyGF2::new(1 << (self.code_length - message_length));
        let remainder = padded % self.generator_poly;
        Ok(padded - remainder)
    }

    pub fn decode(&self, received: PolyGF2) -> Result<PolyGF2, &'static str> {
        let received_length = received.degree() + 1;
        if received_length != self.code_length {
            return Err("Received message has wrong length");
        }

        let mut received_poly_gf2 = received.poly;
        let mut coefficients = vec![];
        for _ in 0..self.code_length {
            coefficients.push(GF2TM::<M>::from(received_poly_gf2 & 1));
            received_poly_gf2 >>= 1;
        }
        let received_poly_gf2m = Polynomial::new(coefficients);
        let syndromes: Vec<_> = (1..self.distance)
            .map(|i| received_poly_gf2m.eval(self.primitive_element.pow(i as u32)))
            .collect();

        let error = if let Some(error_locator) = self.error_locator(syndromes) {
            let error_positions = self.chien_search(error_locator);
            let error_values = error_positions.iter().fold(0, |acc, e| acc ^ (1u32 << e));
            error_values
        } else {
            0
        };
        let corrected = received.poly + error;
        Ok(PolyGF2::new(corrected >> self.generator_poly.degree()))
    }

    fn error_locator(&self, syndromes: Vec<GF2TM<M>>) -> Option<Polynomial<GF2TM<M>>> {
        let t = syndromes.len() / 2;
        for v in (1..=t).rev() {
            let mut matrix = Matrix::<GF2TM<M>>::zero(v, v);
            for i in 0..v {
                for j in 0..v {
                    matrix[[i, j]] = *syndromes.get(i + j).unwrap();
                }
            }
            let right_part = syndromes
                .iter()
                .skip(v)
                .take(v)
                .map(|x| -*x)
                .collect::<Vec<GF2TM<M>>>();
            if let Some(mut solution) = matrix.solve(right_part) {
                solution.push(GF2TM::<M>::one());
                return Some(Polynomial::new(solution));
            }
        }
        None
    }

    fn chien_search(&self, error_locator: Polynomial<GF2TM<M>>) -> Vec<usize> {
        (0..self.code_length)
            .filter(|i| {
                error_locator
                    .eval(self.primitive_element.pow(*i as u32))
                    .is_zero()
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use galois::PolyGF2;

    #[test]
    fn test_encode() {
        const M: u32 = 4;
        let bch = BCH::<M>::from_distance(7);
        let message = PolyGF2::new(0b11011);
        let encoded = bch.encode(message);
        assert_eq!(encoded, Ok(PolyGF2::new(0b110111000010100)));
    }

    #[test]
    fn test_decode_2_err() {
        const M: u32 = 4;
        let bch = BCH::<M>::from_distance(7);
        let message = PolyGF2::new(0b11011);
        let encoded = bch.encode(message).unwrap();
        let err = 0b10000000100000;
        let received = PolyGF2::new(encoded.poly ^ err);
        let decoded = bch.decode(received);
        assert_eq!(decoded, Ok(message));
    }

    #[test]
    fn test_decode_3_err() {
        const M: u32 = 4;
        let bch = BCH::<M>::from_distance(7);
        let message = PolyGF2::new(0b11011);
        let encoded = bch.encode(message).unwrap();
        let err = 0b10010000100000;
        let received = PolyGF2::new(encoded.poly ^ err);
        let decoded = bch.decode(received);
        assert_eq!(decoded, Ok(message));
    }

    #[test]
    
    
    fn test_decode_wrong_length() {
        const M: u32 = 4;
        let bch = BCH::<M>::from_distance(7);
        let received = PolyGF2::new(0b0);
        let decoded = bch.decode(received);
        assert_eq!(decoded, Err("Received message has wrong length"));
    }
}
