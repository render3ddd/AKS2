//! Kronecker substitution ("binary segmentation"): the dissertation's Z version,
//! with the 2021 bugs fixed by construction (REFACTOR.md §1.4, CORRECTIONS.md),
//! tuned in M5.
//!
//! Each polynomial is packed into one big integer by giving every coefficient
//! a fixed-width slot of `b` bits, so A = Σ aᵢ·2^(b·i). One GMP multiply (or
//! square) gives C = A·B, whose slots are exactly the coefficients of the
//! linear product, provided no slot overflows. The slots are then unpacked,
//! folded mod xʳ − 1 and reduced mod n.
//!
//! **Slot width.** Every coefficient is in [0, n), so each coefficient of the
//! linear product is a sum of at most min(U, V) products of at most (n − 1)²,
//! where U and V are the operands' term counts. A slot of b bits with
//! 2ᵇ > min(U, V)·(n − 1)² therefore never overflows, and because nothing is
//! negative there is no sign bit and no borrow: plain unsigned digits.
//!
//! **M5 tuning** (profiled with callgrind; byte-level `mpz_import`/`export`
//! and per-coefficient allocation were over 60% of the run time, the GMP
//! squaring itself under 10%):
//!
//! - slots are bit-granular, not rounded up to bytes, and packing and
//!   unpacking work on 64-bit limbs, which GMP imports and exports with a copy;
//! - slot i and slot i + r are added on limbs before any big-integer work, so
//!   each output coefficient is reduced once;
//! - for n < 2⁶⁴ that reduction is word arithmetic; otherwise it reuses the
//!   output coefficient's allocation;
//! - [`PolyMul::square_into`] writes into an existing polynomial, so a
//!   square-and-multiply loop allocates nothing per coefficient.

use rug::integer::Order;
use rug::{Assign, Integer};

use super::{CyclicPoly, PolyMul, Ring};

/// Kronecker substitution with bit-packed limb slots and a dedicated square.
#[derive(Debug, Clone, Copy, Default)]
pub struct Kronecker;

/// Number of terms up to and including the highest non-zero coefficient.
fn terms(p: &CyclicPoly) -> usize {
    p.coeffs
        .iter()
        .rposition(|c| !c.is_zero())
        .map_or(0, |i| i + 1)
}

/// Bits per slot for a product of polynomials with `u` and `v` terms over ℤₙ:
/// the bit length of min(U, V) · (n − 1)², so 2^bits exceeds it.
fn slot_bits(n: &Integer, u: usize, v: usize) -> usize {
    let nm1 = Integer::from(n - 1u32);
    let bound = Integer::from(nm1.square_ref()) * u.min(v) as u64;
    (bound.significant_bits() as usize).max(1)
}

/// ORs `value` (little-endian limbs) into `buf` starting at bit `off`.
fn put(buf: &mut [u64], off: usize, value: &[u64]) {
    let (w, s) = (off / 64, off % 64);
    for (k, &limb) in value.iter().enumerate() {
        buf[w + k] |= limb << s;
        if s != 0 {
            buf[w + k + 1] |= limb >> (64 - s);
        }
    }
}

/// Reads `width` bits of `src` starting at bit `off` into `out`
/// (`width.div_ceil(64)` limbs; bits past the end of `src` read as 0).
fn get(src: &[u64], off: usize, width: usize, out: &mut [u64]) {
    let (w, s) = (off / 64, off % 64);
    let limbs = width.div_ceil(64);
    let at = |i: usize| src.get(i).copied().unwrap_or(0);
    for (k, o) in out[..limbs].iter_mut().enumerate() {
        *o = if s == 0 {
            at(w + k)
        } else {
            (at(w + k) >> s) | (at(w + k + 1) << (64 - s))
        };
    }
    if width % 64 != 0 {
        out[limbs - 1] &= (1u64 << (width % 64)) - 1;
    }
}

/// a += b on equal-length limb slices (a has room for the carry).
fn add_limbs(a: &mut [u64], b: &[u64]) {
    let mut carry = false;
    for (x, &y) in a.iter_mut().zip(b) {
        let (s1, c1) = x.overflowing_add(y);
        let (s2, c2) = s1.overflowing_add(u64::from(carry));
        *x = s2;
        carry = c1 || c2;
    }
    for x in &mut a[b.len()..] {
        if !carry {
            break;
        }
        let (s, c) = x.overflowing_add(1);
        *x = s;
        carry = c;
    }
}

/// The packed polynomial: Σ cᵢ·2^(b·i) for the first `len` coefficients.
fn pack(ring: &Ring, coeffs: &[Integer], len: usize, b: usize) -> Integer {
    let n_limbs = ring.n().significant_digits::<u64>();
    let mut buf = vec![0u64; (len * b).div_ceil(64) + n_limbs + 1];
    let mut limbs = vec![0u64; n_limbs];
    for (i, c) in coeffs[..len].iter().enumerate() {
        if n_limbs == 1 {
            limbs[0] = c.to_u64_wrapping(); // c < n < 2^64
        } else {
            c.write_digits(&mut limbs, Order::Lsf); // zero-pads to n_limbs
        }
        put(&mut buf, i * b, &limbs);
    }
    Integer::from_digits(&buf, Order::Lsf)
}

/// Unpacks the `len` slots of `product`, folds slot i + r onto slot i, and
/// writes each sum mod n into `out`.
fn unpack_fold_reduce(ring: &Ring, product: &Integer, len: usize, b: usize, out: &mut CyclicPoly) {
    let r = ring.r();
    let n = ring.n();
    let digits = product.to_digits::<u64>(Order::Lsf);
    let slot = b.div_ceil(64);
    let mut lo = vec![0u64; slot + 1];
    let mut hi = vec![0u64; slot];
    let small = n.to_u64();
    out.coeffs.resize(r, Integer::new());
    for (i, c) in out.coeffs.iter_mut().enumerate() {
        if i >= len {
            c.assign(0);
            continue;
        }
        lo[slot] = 0;
        get(&digits, i * b, b, &mut lo);
        let mut j = i + r;
        while j < len {
            get(&digits, j * b, b, &mut hi);
            add_limbs(&mut lo, &hi);
            j += r;
        }
        match small {
            Some(m) => {
                // Horner from the top limb: rem < m < 2^64 throughout.
                let m = u128::from(m);
                let rem = lo
                    .iter()
                    .rev()
                    .fold(0u128, |rem, &l| ((rem << 64) | u128::from(l)) % m);
                c.assign(rem as u64);
            }
            None => {
                c.assign_digits(&lo, Order::Lsf);
                if *c >= *n {
                    *c %= n;
                }
            }
        }
    }
}

impl Kronecker {
    fn product_into(
        &self,
        ring: &Ring,
        a: &CyclicPoly,
        b: Option<&CyclicPoly>,
        out: &mut CyclicPoly,
    ) {
        ring.check(a);
        let u = terms(a);
        let v = b.map_or(u, terms);
        if u == 0 || v == 0 {
            *out = ring.zero(); // the zero-operand case that crashed the 2021 code
            return;
        }
        let bits = slot_bits(ring.n(), u, v);
        let pa = pack(ring, &a.coeffs, u, bits);
        let product = match b {
            Some(b) => {
                ring.check(b);
                pa * pack(ring, &b.coeffs, v, bits)
            }
            None => pa.square(), // mpz_mul with one operand: GMP's squaring path
        };
        unpack_fold_reduce(ring, &product, u + v - 1, bits, out);
    }
}

impl PolyMul for Kronecker {
    fn mul(&self, ring: &Ring, a: &CyclicPoly, b: &CyclicPoly) -> CyclicPoly {
        let mut out = CyclicPoly { coeffs: Vec::new() };
        self.product_into(ring, a, Some(b), &mut out);
        out
    }

    fn square(&self, ring: &Ring, a: &CyclicPoly) -> CyclicPoly {
        let mut out = CyclicPoly { coeffs: Vec::new() };
        self.product_into(ring, a, None, &mut out);
        out
    }

    fn square_into(&self, ring: &Ring, a: &CyclicPoly, out: &mut CyclicPoly) {
        self.product_into(ring, a, None, out);
    }

    fn name(&self) -> &'static str {
        "kronecker"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::Schoolbook;

    #[test]
    fn slot_width_bound() {
        assert_eq!(slot_bits(&Integer::from(256), 1, 5), 16); // 255² = 65025 < 2^16
        assert_eq!(slot_bits(&Integer::from(256), 2, 2), 17); // 130050
        assert_eq!(slot_bits(&Integer::from(2), 1, 1), 1);
    }

    #[test]
    fn bit_put_get_round_trip() {
        for width in [1usize, 7, 63, 64, 65, 127, 128, 129, 200] {
            let count = 9;
            let mut buf = vec![0u64; (count * width).div_ceil(64) + 5];
            let vals: Vec<Vec<u64>> = (0..count)
                .map(|i| {
                    let mut v = vec![u64::MAX - i as u64 * 977; width.div_ceil(64)];
                    if width % 64 != 0 {
                        let l = v.len() - 1;
                        v[l] &= (1u64 << (width % 64)) - 1;
                    }
                    v
                })
                .collect();
            for (i, v) in vals.iter().enumerate() {
                put(&mut buf, i * width, v);
            }
            let mut out = vec![0u64; width.div_ceil(64)];
            for (i, v) in vals.iter().enumerate() {
                get(&buf, i * width, width, &mut out);
                assert_eq!(&out, v, "width {width}, slot {i}");
            }
        }
    }

    #[test]
    fn coefficients_at_the_slot_boundary() {
        // All-(n - 1) operands make every product coefficient hit the bound exactly.
        let ns: Vec<Integer> = [2u64, 3, 255, 256, 257, 65_535, 65_536, 1 << 32, u64::MAX]
            .into_iter()
            .map(Integer::from)
            .chain([
                Integer::from(u64::MAX) + 1u32,
                (Integer::from(1) << 200) - 1u32,
            ])
            .collect();
        for n in ns {
            for r in [1usize, 2, 3, 8, 17] {
                let ring = Ring::new(n.clone(), r);
                let a = ring.from_coeffs(vec![Integer::from(&n - 1u32); r]);
                let want = Schoolbook.mul(&ring, &a, &a);
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
        let mut out = ring.one();
        Kronecker.square_into(&ring, &ring.zero(), &mut out);
        assert_eq!(out, ring.zero());
    }
}
