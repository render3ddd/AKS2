//! The Lenstra (2002) variant of AKS. Milestone M1 provides the parameter
//! choice; the congruence tests and the driver follow in M3.

pub mod params;

pub use params::{ARange, Params, RSearch, a_bound, find_r};
