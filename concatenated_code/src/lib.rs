use bch::BCH;
use galois::{PolyGF2, GF2TM};
use poly_it::Polynomial;
use reed_solomon::ReedSolomon;

pub struct ConcatenatedCode<const M: u32, const N: u32> {
    pub outer_code: ReedSolomon<M>,
    pub inner_code: BCH<N>,
}

impl<const M: u32, const N: u32> ConcatenatedCode<M, N> {
    pub fn encode(&self, message: Polynomial<GF2TM<M>>) -> Vec<PolyGF2> {
        let mut encoded: Vec<PolyGF2> = Vec::new();
        let outer_encoded = self.outer_code.encode(message);
        println!("OUTER ENCODED {outer_encoded:?}");
        for letter in outer_encoded.coeffs().to_vec() {
            encoded.push(
                self.inner_code
                    .encode(
                        letter.value()
                            + PolyGF2 {
                                poly: u32::pow(2, N),
                            },
                    )
                    .unwrap(),
            );
        }
        encoded
    }

    pub fn decode(&self, encoded_message: Vec<PolyGF2>) -> Polynomial<GF2TM<M>> {
        let mut inner_decoded: Vec<GF2TM<M>> = Vec::new();
        for letter in encoded_message {
            // let decoded_letter = self.inner_code.decode(letter).unwrap().poly;
            inner_decoded.push(GF2TM::<M>::new(self.inner_code.decode(letter).unwrap()));
        }
        println!("INNER DECODED {inner_decoded:?}");

        let poly = Polynomial::new(inner_decoded);
        let outer_decoded = self.outer_code.decode(poly);
        outer_decoded
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use galois::PolyGF2;
    #[test]
    fn test_encode() {
        const M: u32 = 8;
        let reed_solomon = ReedSolomon::<M> { distance: 5 };
        const N: u32 = 4;
        let bch: BCH<N> = BCH::<N>::from_distance(7).unwrap();

        let concatenated_code = ConcatenatedCode {
            outer_code: reed_solomon,
            inner_code: bch,
        };

        let poly_msg = Polynomial::new(vec![
            GF2TM::<M>::new(PolyGF2 { poly: 3 as u32 }),
            GF2TM::<M>::new(PolyGF2 { poly: 2 as u32 }),
        ]);

        let encoded = concatenated_code.encode(poly_msg);
        let true_encoded = vec![
            PolyGF2 {
                poly: 0b100110111000010,
            },
            PolyGF2 {
                poly: 0b100011110101100,
            },
            PolyGF2 {
                poly: 0b101110000101001,
            },
            PolyGF2 {
                poly: 0b101011001000111,
            },
            PolyGF2 {
                poly: 0b110111000010100,
            },
            PolyGF2 {
                poly: 0b110010001111010,
            },
            PolyGF2 {
                poly: 0b111111111111111,
            },
        ];
        assert_eq!(encoded, true_encoded);
    }

    #[test]
    fn test_decode_0_err() {
        const M: u32 = 4;
        let reed_solomon = ReedSolomon::<M> { distance: 3 };
        const N: u32 = 4;
        let bch: BCH<N> = BCH::<N>::from_distance(7).unwrap();

        let concatenated_code = ConcatenatedCode {
            outer_code: reed_solomon,
            inner_code: bch,
        };

        let poly_msg = Polynomial::new(vec![
            GF2TM::<M>::new(PolyGF2 { poly: 2 }),
            GF2TM::<M>::new(PolyGF2 { poly: 3 }),
        ]);

        let encoded = concatenated_code.encode(poly_msg.clone());
        let decoded = concatenated_code.decode(encoded);
        assert_eq!(decoded, poly_msg);
    }

    #[test]
    fn test_decode_with_errors() {
        const M: u32 = 4;
        let reed_solomon = ReedSolomon::<M> { distance: 3 };
        const N: u32 = 4;
        let bch: BCH<N> = BCH::<N>::from_distance(7).unwrap();

        let concatenated_code = ConcatenatedCode {
            outer_code: reed_solomon,
            inner_code: bch,
        };

        let poly_msg = Polynomial::new(vec![
            GF2TM::<M>::new(PolyGF2 { poly: 2 }),
            GF2TM::<M>::new(PolyGF2 { poly: 3 }),
        ]);

        let mut encoded = concatenated_code.encode(poly_msg.clone());
        encoded[1] = PolyGF2 {
            poly: 0b100001110100100,
        };
        encoded[3] = PolyGF2 {
            poly: 0b101110000000000,
        };
        let decoded = concatenated_code.decode(encoded);
        assert_eq!(decoded, poly_msg);
    }
}
