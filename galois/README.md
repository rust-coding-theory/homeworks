# Galois Crate

This crate implements the Galois field arithmetic for $GF(2^m)$.
The implementation is based on the fact that $GF(2^m) \cong GF(2)[x]\Big/(P_m(x))$, where $P_m(x)$ is an irreducible polynomial of degree $m$.

Structs:
- [`PolyGF2`](https://github.com/rust-coding-theory/homeworks/blob/main/galois/src/poly_gf2.rs): Represents a polynomial with coefficients in $GF(2)$.
- [`GF2TM<const M: u32>`](https://github.com/rust-coding-theory/homeworks/blob/main/galois/src/gf2tm.rs): Represents an element in $GF(2^m)$, where $m$ is a const generic parameter.
- [`Matrix<T>`](https://github.com/rust-coding-theory/homeworks/blob/main/galois/src/matrix.rs): Represents a matrix with elements of type `T`. Doesn't need to belong to this crate actually but I put it here for now.

## Usage

For now please refer to `mod tests` in the source code.
