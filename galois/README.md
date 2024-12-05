# Galois Crate

This crate implements the Galois field arithmetic for $GF(2^m)$.
The implementation is based on the fact that $GF(2^m) \cong GF(2)[x]\Big/(P_m(x))$, where $P_m(x)$ is an irreducible polynomial of degree $m$.

Structs:
- `PolyGF2`: Represents a polynomial with coefficients in $GF(2)$.
- `GF2M`: Represents an element in $GF(2^m)$.
- `GF2TM<const M: u32>`: Represents an element in $GF(2^m)$, where $m$ is a const generic parameter. This is the same as the previous one but I worked on this implementation more :)
- `Matrix<T>`: Represents a matrix with elements of type `T`. Doesn't need to belong to this crate actually but I put it here for now.

## Usage

For now please refer to `mod tests` in the source code. I'll be implementing the examples in the crate soon.