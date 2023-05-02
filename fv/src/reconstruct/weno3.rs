//! WENO3 interpolation
//! Reference: https://en.wikipedia.org/wiki/WENO
//! \todo use precomputed coefficients

/// | x_{-1} | x_0 | x_1 |
///          ^
///          |
///          return value
pub fn interp_weno3(phim1: f64, phi: f64, phip1: f64) -> f64 {
    let p0 = (1.0 / 2.0) * phi + (1.0 / 2.0) * phim1;
    let p1 = (-1.0 / 2.0) * phip1 + (3.0 / 2.0) * phi;

    let beta0 = (phim1 - phi).powi(2);
    let beta1 = (phi - phip1).powi(2);

    let alpha0 = (1.0 / 3.0) / ((beta0 + 1e-10) * (beta0 + 1.0e-10));
    let alpha1 = (2.0 / 3.0) / ((beta1 + 1e-10) * (beta1 + 1.0e-10));

    let alpha_sum_inv = 1.0 / (alpha0 + alpha1);

    let w0 = alpha0 * alpha_sum_inv;
    let w1 = alpha1 * alpha_sum_inv;

    w0 * p0 + w1 * p1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interp_weno3() {
        let phim1 = 1.0;
        let phi = 2.0;
        let phip1 = 3.0;
        let result = interp_weno3(phim1, phi, phip1);
        let expected_result = 1.5;
        approx::assert_abs_diff_eq!(
            result,
            expected_result,
            epsilon = f64::EPSILON
        );
    }
}
