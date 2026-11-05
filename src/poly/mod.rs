//! The ring ℤₙ[x]/(xʳ − 1), and interchangeable ways of multiplying in it.
//!
//! AKS needs exactly one polynomial ring, so rather than a general polynomial
//! library this module owns that ring outright (RUST_PLAN.md §3.2). Reducing
//! mod xʳ − 1 is just folding coefficient i + r onto i.
//!
//! - [`Ring`] holds n and r and does everything except multiplication.
//! - [`PolyMul`] is the multiplier interface, with back-ends [`Schoolbook`]
//!   (the oracle for tests) and [`Kronecker`] (binary segmentation: the
//!   dissertation's "Z version").

use std::fmt;

use rug::Integer;

mod kronecker;
mod schoolbook;

pub use kronecker::Kronecker;
pub use schoolbook::Schoolbook;

/// An element of ℤₙ[x]/(xʳ − 1): exactly r coefficients, each in [0, n).
///
/// Only a [`Ring`] can build one, so the invariant always holds for the ring
/// that made it. Mixing polynomials from different rings is a logic error,
/// caught by debug assertions.
#[derive(Clone, PartialEq, Eq)]
pub struct CyclicPoly {
    coeffs: Vec<Integer>,
}

impl CyclicPoly {
    /// The coefficients, constant term first; always r of them.
    pub fn coeffs(&self) -> &[Integer] {
        &self.coeffs
    }

    /// Whether this is the zero polynomial.
    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|c| c.is_zero())
    }
}

impl fmt::Debug for CyclicPoly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Sparse form, highest power first: e.g. [x^5 + 3x + 2].
        let terms: Vec<String> = self
            .coeffs
            .iter()
            .enumerate()
            .rev()
            .filter(|(_, c)| !c.is_zero())
            .map(|(i, c)| match i {
                0 => c.to_string(),
                1 => format!("{c}x"),
                _ => format!("{c}x^{i}"),
            })
            .collect();
        if terms.is_empty() {
            write!(f, "[0]")
        } else {
            write!(f, "[{}]", terms.join(" + "))
        }
    }
}

/// The ring ℤₙ[x]/(xʳ − 1), for n ≥ 2 and r ≥ 1.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ring {
    n: Integer,
    r: usize,
}

impl Ring {
    /// # Panics
    /// If n < 2 or r = 0.
    pub fn new(n: Integer, r: usize) -> Self {
        assert!(n >= 2, "modulus n must be at least 2");
        assert!(r >= 1, "r must be at least 1");
        Ring { n, r }
    }

    /// The coefficient modulus n.
    pub fn n(&self) -> &Integer {
        &self.n
    }

    /// The degree r of xʳ − 1, which is also the length of every element.
    pub fn r(&self) -> usize {
        self.r
    }

    /// 0.
    pub fn zero(&self) -> CyclicPoly {
        CyclicPoly {
            coeffs: vec![Integer::new(); self.r],
        }
    }

    /// 1.
    pub fn one(&self) -> CyclicPoly {
        self.monomial(0, Integer::from(1))
    }

    /// c·xᵏ, reducing k mod r and c mod n.
    pub fn monomial(&self, k: usize, c: Integer) -> CyclicPoly {
        let mut p = self.zero();
        p.coeffs[k % self.r] = self.reduce(c);
        p
    }

    /// xᵏ, for any k ≥ 0 (as a big integer, since AKS needs xⁿ = x^(n mod r)).
    pub fn x_pow(&self, k: &Integer) -> CyclicPoly {
        let r = u32::try_from(self.r).expect("r fits u32");
        self.monomial(k.mod_u(r) as usize, Integer::from(1))
    }

    /// Builds an element from any integer coefficients, constant term first:
    /// folds index i onto i mod r and reduces each coefficient into [0, n).
    pub fn from_coeffs<I>(&self, coeffs: I) -> CyclicPoly
    where
        I: IntoIterator,
        I::Item: Into<Integer>,
    {
        let mut p = self.zero();
        for (i, c) in coeffs.into_iter().enumerate() {
            p.coeffs[i % self.r] += c.into();
        }
        for c in &mut p.coeffs {
            *c = self.reduce(std::mem::take(c));
        }
        p
    }

    /// c mod n, in [0, n).
    fn reduce(&self, mut c: Integer) -> Integer {
        if c.is_negative() || c >= self.n {
            c.modulo_mut(&self.n);
        }
        c
    }

    fn check(&self, p: &CyclicPoly) {
        debug_assert_eq!(p.coeffs.len(), self.r, "polynomial from a different ring");
        debug_assert!(p.coeffs.iter().all(|c| !c.is_negative() && *c < self.n));
    }

    /// a + b.
    pub fn add(&self, a: &CyclicPoly, b: &CyclicPoly) -> CyclicPoly {
        self.check(a);
        self.check(b);
        let coeffs = a
            .coeffs
            .iter()
            .zip(&b.coeffs)
            .map(|(x, y)| {
                let mut s = Integer::from(x + y);
                if s >= self.n {
                    s -= &self.n;
                }
                s
            })
            .collect();
        CyclicPoly { coeffs }
    }

    /// a − b.
    pub fn sub(&self, a: &CyclicPoly, b: &CyclicPoly) -> CyclicPoly {
        self.check(a);
        self.check(b);
        let coeffs = a
            .coeffs
            .iter()
            .zip(&b.coeffs)
            .map(|(x, y)| {
                let mut d = Integer::from(x - y);
                if d.is_negative() {
                    d += &self.n;
                }
                d
            })
            .collect();
        CyclicPoly { coeffs }
    }

    /// Folds a raw (unreduced, non-negative) product of any length onto r
    /// coefficients and reduces each mod n once, after folding: the "lazy
    /// reduction" of RUST_PLAN.md §3.2, which does r reductions instead of
    /// 2r − 1.
    fn fold_reduce(&self, raw: Vec<Integer>) -> CyclicPoly {
        let mut coeffs = raw;
        if coeffs.len() > self.r {
            let (low, high) = coeffs.split_at_mut(self.r);
            for (i, c) in high.iter_mut().enumerate() {
                low[i % self.r] += std::mem::take(c);
            }
            coeffs.truncate(self.r);
        } else {
            coeffs.resize(self.r, Integer::new());
        }
        for c in &mut coeffs {
            debug_assert!(!c.is_negative());
            if *c >= self.n {
                *c %= &self.n;
            }
        }
        CyclicPoly { coeffs }
    }

    /// baseᵉ by left-to-right square-and-multiply, for e ≥ 0.
    pub fn pow<M: PolyMul + ?Sized>(&self, m: &M, base: &CyclicPoly, e: &Integer) -> CyclicPoly {
        assert!(!e.is_negative(), "negative exponent");
        self.check(base);
        if e.is_zero() {
            return self.one();
        }
        let mut acc = base.clone();
        for i in (0..e.significant_bits() - 1).rev() {
            acc = m.square(self, &acc);
            if e.get_bit(i) {
                acc = m.mul(self, &acc, base);
            }
        }
        acc
    }
}

/// A way of multiplying two elements of a [`Ring`].
///
/// Implementations must return exactly the ring product. Every back-end is
/// property-tested against [`Schoolbook`] (RUST_PLAN.md §7 rule 7).
pub trait PolyMul: Sync {
    /// a · b in the ring.
    fn mul(&self, ring: &Ring, a: &CyclicPoly, b: &CyclicPoly) -> CyclicPoly;

    /// a² in the ring. Back-ends override this when squaring is cheaper.
    fn square(&self, ring: &Ring, a: &CyclicPoly) -> CyclicPoly {
        self.mul(ring, a, a)
    }

    /// Short name for logs and the command line.
    fn name(&self) -> &'static str;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring(n: u64, r: usize) -> Ring {
        Ring::new(Integer::from(n), r)
    }

    #[test]
    fn construction_folds_and_reduces() {
        let z = ring(7, 3);
        // 1 + 2x + 3x² + 4x³ + 10x⁴ - 1x⁵  ->  (1+4) + (2+10)x + (3-1)x²  mod 7
        let p = z.from_coeffs([1i32, 2, 3, 4, 10, -1]);
        assert_eq!(p.coeffs(), [5, 5, 2]);
        assert_eq!(z.x_pow(&Integer::from(10)), z.monomial(1, Integer::from(1)));
        assert_eq!(z.monomial(4, Integer::from(-1)).coeffs(), [0, 6, 0]);
        assert_eq!(format!("{p:?}"), "[2x^2 + 5x + 5]");
        assert_eq!(format!("{:?}", z.zero()), "[0]");
    }

    #[test]
    fn add_sub_wrap() {
        let z = ring(5, 2);
        let a = z.from_coeffs([4, 3]);
        let b = z.from_coeffs([3, 4]);
        assert_eq!(z.add(&a, &b).coeffs(), [2, 2]);
        assert_eq!(z.sub(&b, &a).coeffs(), [4, 1]);
        assert_eq!(z.sub(&a, &a), z.zero());
    }

    #[test]
    fn pow_small_exponents() {
        let z = ring(1_000_003, 5);
        let x_minus_2 = z.from_coeffs([-2, 1]);
        assert_eq!(z.pow(&Schoolbook, &x_minus_2, &Integer::from(0)), z.one());
        assert_eq!(z.pow(&Schoolbook, &x_minus_2, &Integer::from(1)), x_minus_2);
        // (x - 2)^3 = x^3 - 6x^2 + 12x - 8
        let want = z.from_coeffs([-8, 12, -6, 1]);
        assert_eq!(z.pow(&Schoolbook, &x_minus_2, &Integer::from(3)), want);
    }

    #[test]
    fn frobenius_holds_for_a_prime() {
        // (x - a)^p = x^(p mod r) - a in Z_p[x]/(x^r - 1): the AKS identity.
        let p = Integer::from(10007);
        let z = Ring::new(p.clone(), 13);
        for a in 1..5 {
            let lhs = z.pow(&Kronecker, &z.from_coeffs([-a, 1]), &p);
            let rhs = z.sub(&z.x_pow(&p), &z.monomial(0, Integer::from(a)));
            assert_eq!(lhs, rhs, "a = {a}");
        }
    }
}
