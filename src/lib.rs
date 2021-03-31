//! # Bell Number
//!
//! `dahl_bellnumber` is a collection of functions related to the [Bell number](https://en.wikipedia.org/wiki/Bell_number),
//! which gives the number of partitions of a set.
//!
//!

#![allow(dead_code)]

#[cfg(test)]
#[macro_use]
extern crate approx;
extern crate num_bigint;

use num_bigint::BigUint;
use num_rational::Ratio;
use num_traits::cast::ToPrimitive;
use num_traits::{One, Zero};
use std::convert::TryFrom;

/// Compute the [Bell number](https://en.wikipedia.org/wiki/Bell_number).
///
/// # Examples
///
/// ```
/// let answer = dahl_bellnumber::bell(5);
///
/// use std::convert::TryFrom;
/// use num_traits::cast::ToPrimitive;
///
/// assert_eq!(answer, num_bigint::BigUint::try_from(52_u32).unwrap());
/// assert_eq!(answer.to_f64().unwrap(), 52.0);
/// ```
pub fn bell(n: usize) -> BigUint {
    let mut r1: Vec<BigUint> = vec![Zero::zero(); n];
    let mut r2: Vec<BigUint> = vec![Zero::zero(); n];
    r1[0] = One::one();
    for k in 1..n {
        r2[0] = r1[k - 1].clone();
        for i in 1..(k + 1) {
            r2[i] = r1[i - 1].clone() + &r2[i - 1];
        }
        std::mem::swap(&mut r1, &mut r2);
    }
    r1[n - 1].clone()
}

/// Compute the natural logarithm of the [Bell number](https://en.wikipedia.org/wiki/Bell_number).
///
/// # Examples
///
/// ```
/// let answer = dahl_bellnumber::lbell(5);
///
/// assert!( (answer - 52.0_f64.ln()).abs() < 0.00000001 );
/// ```
pub fn lbell(n: usize) -> f64 {
    let value = bell(n);
    let n_bits = value.bits();
    let threshold = 1022_u64;
    let log2 = if n_bits > threshold {
        let n_shifted_bits = value.bits() - threshold;
        let shifted_value = value >> n_shifted_bits;
        if shifted_value.bits() > threshold {
            return f64::INFINITY;
        }
        let y: f64 = shifted_value.to_f64().unwrap();
        (n_shifted_bits as f64) + y.log2()
    } else {
        value.to_f64().unwrap().log2()
    };
    log2 / std::f64::consts::LOG2_E
}

struct UniformDistributionCache(Vec<Vec<BigUint>>);

impl UniformDistributionCache {
    pub fn new(n: usize) -> Self {
        let mut x: Vec<Vec<BigUint>> = (0..(n+1)).map(|k| vec![One::one(); n - k + 1]).collect();
        for k in (0..(n-1)).rev() {
            for r in 1..(n - k + 1) {
                x[k][r] = &x[k][r - 1] * (k + 1) + &x[k + 1][r - 1];
            }
        }
        // println!("{:?}",x);
        Self(x)
    }

    pub fn bell(&self, n: usize) -> BigUint {
        self.partition_counter(n-1, 1)
    }

    pub fn partition_counter(
        &self,
        n_remaining_items_after_allocation: usize,
        n_clusters_after_allocation: usize,
    ) -> BigUint {
        self.0[n_clusters_after_allocation-1][n_remaining_items_after_allocation].clone()
    }

    pub fn probs_for_uniform(&self, n_remaining_items_after_allocation: usize, n_clusters_after_allocation: usize) -> (f64, f64) {
        let a = self.partition_counter(n_remaining_items_after_allocation - 1, n_clusters_after_allocation);
        let b = self.partition_counter(n_remaining_items_after_allocation - 1, n_clusters_after_allocation + 1);
        let denominator = self.partition_counter(n_remaining_items_after_allocation, n_clusters_after_allocation);
        let left = Ratio::new(a, denominator.clone()).to_f64().unwrap();
        let right = Ratio::new(b, denominator).to_f64().unwrap();
        (left, right)
    }
}

/// C-friendly wrapper over the `bell` function.
#[doc(hidden)]
#[no_mangle]
pub extern "C" fn dahl_bellnumber__bell(n: i32) -> f64 {
    if n < 0 {
        return 0.0;
    }
    match usize::try_from(n) {
        Ok(n) => match bell(n).to_f64() {
            Some(x) => x,
            None => f64::INFINITY,
        },
        Err(_) => f64::INFINITY,
    }
}

/// C-friendly wrapper over the `lbell` function.
#[doc(hidden)]
#[no_mangle]
pub extern "C" fn dahl_bellnumber__lbell(n: i32) -> f64 {
    if n < 0 {
        return 0.0;
    }
    match usize::try_from(n) {
        Ok(n) => lbell(n),
        Err(_) => f64::INFINITY,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lbell() {
        assert_relative_eq!(lbell(220), 714.4032630589774);
        assert_relative_eq!(bell(5).to_f64().unwrap(), 52.0);
    }

    #[test]
    fn test_bell_generalized_cache() {
        let cache = UniformDistributionCache::new(100);
        //let cache = UniformDistributionCache::new(10);
        assert_eq!(cache.bell(10), bell(10));
        assert_eq!(cache.partition_counter(4, 1), BigUint::from(52_u8));
        assert_eq!(cache.partition_counter(3, 1), BigUint::from(15_u8));
        assert_eq!(cache.partition_counter(3, 2), BigUint::from(37_u8));
        assert_eq!(cache.partition_counter(2, 1), BigUint::from(5_u8));
        assert_eq!(cache.partition_counter(2, 2), BigUint::from(10_u8));
        assert_eq!(cache.partition_counter(2, 3), BigUint::from(17_u8));
        assert_eq!(cache.partition_counter(1, 1), BigUint::from(2_u8));
        assert_eq!(cache.partition_counter(1, 2), BigUint::from(3_u8));
        assert_eq!(cache.partition_counter(1, 3), BigUint::from(4_u8));
        assert_eq!(cache.partition_counter(1, 4), BigUint::from(5_u8));
        assert_eq!(cache.partition_counter(0, 1), BigUint::from(1_u8));
        assert_eq!(cache.partition_counter(0, 2), BigUint::from(1_u8));
        assert_eq!(cache.partition_counter(0, 3), BigUint::from(1_u8));
        assert_eq!(cache.partition_counter(0, 4), BigUint::from(1_u8));
        assert_eq!(cache.partition_counter(0, 5), BigUint::from(1_u8));
        let c = 20;
        let (a, b) = cache.probs_for_uniform(80, c);
        assert_eq!((c as f64) * a + b, 1.0);
        let c = 4;
        let (a, b) = cache.probs_for_uniform(50, c);
        assert_eq!((c as f64) * a + b, 1.0);
    }
}
