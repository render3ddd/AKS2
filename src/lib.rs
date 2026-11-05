//! `aks2`: the Lenstra (2002) variant of the AKS deterministic primality test.
//!
//! This crate is a Rust rewrite of the corrected C++ implementation
//! (<https://github.com/render3d/AKS>, tag `v1.1`). See `RUST_PLAN.md` for the
//! milestones; modules are added as each milestone lands.
//!
//! - [`num`]: ⌊log₂ n⌋, perfect powers, multiplicative order, φ and λ (M1).
//! - [`aks`]: choosing r and the range of *a* (M1); the test itself (M3).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod aks;
pub mod num;

/// Crate version, as set in `Cargo.toml`.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
