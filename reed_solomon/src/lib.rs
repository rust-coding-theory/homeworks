use galois::{Matrix, PolyGF2, GF2TM};
use poly_it::num_traits::Zero;
use poly_it::Polynomial;

pub struct ReedSolomon<const M: u32> {
    pub distance: usize,
}

impl<const M: u32> ReedSolomon<M> {
    pub fn encode(&self, message: Polynomial<GF2TM<M>>) -> Polynomial<GF2TM<M>> {
        let mut encoded: Vec<GF2TM<M>> = Vec::new();
        for i in 0..message.coeffs().len() + self.distance {
            encoded.push(message.eval(GF2TM::<M>::new(PolyGF2 { poly: i as u32 })))
        }
        Polynomial::new(encoded)
    }

    pub fn decode(&self, encoded_message: Polynomial<GF2TM<M>>) -> Polynomial<GF2TM<M>> {
        let mut e = self.max_num_of_errors();
        let (a, rhs) = loop {
            let (lhs, rhs) = self.create_linear_system_with_err_num(e, &encoded_message);
            let a: Matrix<GF2TM<M>> = Matrix::new(
                lhs,
                encoded_message.coeffs().len(),
                encoded_message.coeffs().len(),
            );
            if a.determinant().is_zero() {
                e -= 1;
                continue;
            } else {
                break (a, rhs);
            }
        };
        let result = a.solve(rhs).unwrap();
        let mut e_vec = result[0..e].to_vec();
        e_vec.push(GF2TM::<M>::new(PolyGF2 { poly: 1 }));
        let e_poly = Polynomial::new(e_vec);
        let q_poly = Polynomial::new(result[e..].to_vec());
        let decoded = q_poly / e_poly;
        decoded.0
    }

    pub fn max_num_of_errors(&self) -> usize {
        (self.distance - 1) / 2
    }

    fn create_linear_system_with_err_num(
        &self,
        errors_num: usize,
        encoded_message: &Polynomial<GF2TM<M>>,
    ) -> (Vec<GF2TM<M>>, Vec<GF2TM<M>>) {
        let mut lhs: Vec<GF2TM<M>> = Vec::new();
        let mut rhs: Vec<GF2TM<M>> = Vec::new();
        for i in 0..(encoded_message.coeffs().len()) {
            for j in 0..errors_num {
                lhs.push(
                    encoded_message.coeffs()[i]
                        * ((GF2TM::<M>::new(PolyGF2 { poly: i as u32 })).pow(j as u32)),
                );
            }
            rhs.push(
                -((GF2TM::<M>::new(PolyGF2 { poly: i as u32 }).pow(errors_num as u32))
                    * encoded_message.coeffs()[i]),
            );
            for j in 0..(encoded_message.coeffs().len() - errors_num) {
                lhs.push(-(GF2TM::<M>::new(PolyGF2 { poly: (i) as u32 }).pow(j as u32)))
            }
        }
        (lhs, rhs)
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
        let poly_msg = Polynomial::new(vec![
            GF2TM::new(PolyGF2 { poly: 3 as u32 }),
            GF2TM::new(PolyGF2 { poly: 2 as u32 }),
            GF2TM::new(PolyGF2 { poly: 8 as u32 }),
        ]);

        let encoded = reed_solomon.encode(poly_msg);
        println!("{encoded:?}");
        let true_encoded = Polynomial::new(vec![
            GF2TM::<M>::new(PolyGF2 { poly: 3 as u32 }),
            GF2TM::<M>::new(PolyGF2 { poly: 9 as u32 }),
            GF2TM::<M>::new(PolyGF2 { poly: 39 as u32 }),
            GF2TM::<M>::new(PolyGF2 { poly: 45 as u32 }),
            GF2TM::<M>::new(PolyGF2 { poly: 139 as u32 }),
            GF2TM::<M>::new(PolyGF2 { poly: 129 as u32 }),
            GF2TM::<M>::new(PolyGF2 { poly: 175 as u32 }),
            GF2TM::<M>::new(PolyGF2 { poly: 165 as u32 }),
        ]);
        assert_eq!(encoded, true_encoded);
    }

    #[test]
    fn test_decode_0_err() {
        const M: u32 = 8;
        let reed_solomon = ReedSolomon::<M> { distance: 3 };
        let poly_msg = Polynomial::new(vec![
            GF2TM::new(PolyGF2 { poly: 36 as u32 }),
            GF2TM::new(PolyGF2 { poly: 2 as u32 }),
        ]);

        let encoded = reed_solomon.encode(poly_msg.clone());
        let decoded = reed_solomon.decode(encoded);
        assert_eq!(decoded, poly_msg);
    }

    #[test]
    fn test_decode_1_err() {
        const M: u32 = 8;
        let reed_solomon = ReedSolomon::<M> { distance: 3 };
        let poly_msg = Polynomial::new(vec![
            GF2TM::new(PolyGF2 { poly: 36 as u32 }),
            GF2TM::new(PolyGF2 { poly: 2 as u32 }),
        ]);

        let encoded = reed_solomon.encode(poly_msg.clone());
        let mut encoded_coefs = encoded.coeffs().to_vec();
        encoded_coefs[0] = GF2TM::<M>::new(PolyGF2 { poly: 44 });

        let decoded = reed_solomon.decode(Polynomial::new(encoded_coefs));
        assert_eq!(decoded, poly_msg);
    }
}
