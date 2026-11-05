//! Kronecker substitution ("binary segmentation"): the dissertation's Z version,
//! with the 2021 bugs fixed by construction (REFACTOR.md §1.4, CORRECTIONS.md).
//!
//! Each polynomial is packed into one big integer by giving every coefficient
//! a fixed-width slot of `w` bytes, so A = Σ aᵢ·2^(8wi). One GMP multiply (or
//! square) gives C = A·B, whose slots are exactly the coefficients of the
//! linear product, provided no slot overflows. The slots are then unpacked,
//! folded mod xʳ − 1 and reduced mod n.
//!
//! **Slot width.** Every coefficient is in [0, n), so each coefficient of the
//! linear product is a sum of at most min(U, V) products of at most (n − 1)²,
//! where U and V are the operands' term counts. A slot of b bits with
//! 2ᵇ > min(U, V)·(n − 1)² therefore never overflows, and because nothing is
//! negative there is no sign bit and no borrow: plain unsigned digits. Slots
//! are rounded up to whole bytes, as in the corrected C++ (`ZZpXmultiply`), so
//! packing and unpacking are straight byte copies (GMP's `mpz_import`/`export`).

use rug::Integer;
use rug::integer::Order;

use super::{CyclicPoly, PolyMul, Ring};

/// Kronecker substitution with byte-aligned slots and a dedicated square.
#[derive(Debug, Clone, Copy, Default)]
pub struct Kronecker;

/// Number of terms up to and including the highest non-zero coefficient.
fn terms(p: &CyclicPoly) -> usize {
    p.coeffs
        .iter()
        .rposition(|c| !c.is_zero())
        .map_or(0, |i| i + 1)
}

/// Bytes per slot for a product of polynomials with `u` and `v` terms over ℤₙ.
fn slot_bytes(n: &Integer, u: usize, v: usize) -> usize {
    // bound = min(U, V) · (n − 1)², the largest possible product coefficient.
    let nm1 = Integer::from(n - 1u32);
    let bound = Integer::from(nm1.square_ref()) * u.min(v) as u64;
    // 2^bits > bound, so the slot holds every value 0..=bound.
    (bound.significant_bits() as usize).div_ceil(8).max(1)
}

/// A = Σ cᵢ · 2^(8·w·i), for the first `len` coefficients.
fn pack(coeffs: &[Integer], len: usize, w: usize) -> Integer {
    let mut bytes = vec![0u8; len * w];
    for (c, slot) in coeffs[..len].iter().zip(bytes.chunks_exact_mut(w)) {
        // c < n ≤ 2^(8w), so it always fits its slot.
        c.write_digits(slot, Order::Lsf);
    }
    Integer::from_digits(&bytes, Order::Lsf)
}

/// Splits C back into `len` slot values of `w` bytes each.
fn unpack(c: &Integer, len: usize, w: usize) -> Vec<Integer> {
    let mut bytes = c.to_digits::<u8>(Order::Lsf);
    debug_assert!(bytes.len() <= len * w, "slot overflow: bound too small");
    bytes.resize(len * w, 0);
    bytes
        .chunks_exact(w)
        .map(|s| Integer::from_digits(s, Order::Lsf))
        .collect()
}

impl PolyMul for Kronecker {
    fn mul(&self, ring: &Ring, a: &CyclicPoly, b: &CyclicPoly) -> CyclicPoly {
        ring.check(a);
        ring.check(b);
        let (u, v) = (terms(a), terms(b));
        if u == 0 || v == 0 {
            return ring.zero(); // the zero-operand case that crashed the 2021 code
        }
        let w = slot_bytes(ring.n(), u, v);
        let product = pack(&a.coeffs, u, w) * pack(&b.coeffs, v, w);
        ring.fold_reduce(unpack(&product, u + v - 1, w))
    }

    fn square(&self, ring: &Ring, a: &CyclicPoly) -> CyclicPoly {
        ring.check(a);
        let u = terms(a);
        if u == 0 {
            return ring.zero();
        }
        let w = slot_bytes(ring.n(), u, u);
        let packed = pack(&a.coeffs, u, w);
        let product = Integer::from(packed.square_ref()); // mpz_mul with one operand: GMP's squaring path
        ring.fold_reduce(unpack(&product, 2 * u - 1, w))
    }

    fn name(&self) -> &'static str {
        "kronecker"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slot_width_bound() {
        // n = 256: (n - 1)^2 = 65025 < 2^16, times 1 term -> 2 bytes.
        assert_eq!(slot_bytes(&Integer::from(256), 1, 5), 2);
        // times 2 terms -> 130050 needs 17 bits -> 3 bytes.
        assert_eq!(slot_bytes(&Integer::from(256), 2, 2), 3);
        assert_eq!(slot_bytes(&Integer::from(2), 1, 1), 1);
    }

    #[test]
    fn coefficients_at_the_slot_boundary() {
        // All-(n - 1) operands make every product coefficient hit the bound exactly.
        for n in [2u64, 3, 255, 256, 257, 65_535, 65_536, 1 << 32, u64::MAX] {
            for r in [1usize, 2, 3, 8, 17] {
                let ring = Ring::new(Integer::from(n), r);
                let a = ring.from_coeffs(vec![Integer::from(n - 1); r]);
                let want = super::super::Schoolbook.mul(&ring, &a, &a);
                assert_eq!(Kronecker.mul(&ring, &a, &a), want, "n = {n}, r = {r}");
                assert_eq!(Kronecker.square(&ring, &a), want, "n = {n}, r = {r}");
            }
        }
    }

    #[test]
    fn zero_and_constants() {
        let ring = Ring::new(Integer::from(10007), 7);
        let three = ring.from_coeffs([3]);
        assert_eq!(Kronecker.mul(&ring, &ring.zero(), &three), ring.zero());
        assert_eq!(Kronecker.square(&ring, &ring.zero()), ring.zero());
        assert_eq!(Kronecker.mul(&ring, &three, &three), ring.from_coeffs([9]));
    }
}
