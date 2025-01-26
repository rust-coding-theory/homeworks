# Reed-Solomon Crate

This crate contains original view implementation of Reed-Solomon error correction codes. __Original view implementation__ means that encoder uses letters of message as coefficients of polynomial and get codeword computing the values in some points.

## Usage

```rust
const M: u32 = 8; // 2^M is the power of the alphabet
let reed_solomon = ReedSolomon::<M> { distance: 5 }; // maximum number of errors to be corrected is equal to (distance - 1) / 2

// create polynomial from message
let poly_msg = Polynomial::new(vec![
    GF2TM::new(PolyGF2 { poly: 3 as u32 }),
    GF2TM::new(PolyGF2 { poly: 2 as u32 }),
    GF2TM::new(PolyGF2 { poly: 8 as u32 }),
]);

let encoded = reed_solomon.encode(poly_msg);

// then you can add some errors to encoded message (like it was sent through a noisy channel)

let decoded = reed_solomon.decode(encoded);
```