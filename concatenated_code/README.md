# Concatenated Code Crate

If you want to encode your message using two different codes in order to increase number of errors which can be corrected, you can use the code from this crate. It concatenates Reed-Solomon codes as outer code and BCH as inner code. 

## Usage

```rust
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

// add some noise

let decoded = concatenated_code.decode(encoded);
```

More examples you can find in the tests module in this crate. More detailed documentation about ReedSolomon and BCH usage is contained [here](../reed_solomon/README.md) and [here](../bch/README.md) respectively.