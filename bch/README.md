# BCH Crate

## Usage

```rs
const M: u32 = 4; // or other constant which specifies the field GF(2^M)
let max_errors = 3; // max tolerated number of errors
let bch = BCH::<M>::from_max_errors(3).unwrap() // returns Err if the max_errors specified is too high for the given M

let message = PolyGF2::new(0b11011);
let encoded = bch.encode(message).unwrap(); // returns Err if message is too long. The allowed message length can be obtained with bch.max_message_length()

let err = 0b10010000100000; // error simulation
let received = PolyGF2::new(encoded.poly ^ err);

let decoded = bch.decode(received);
```
