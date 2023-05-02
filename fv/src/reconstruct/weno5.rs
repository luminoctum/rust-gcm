//! WENO5 interpolation
//! Reference: https://en.wikipedia.org/wiki/WENO
//! \todo use precomputed coefficients

/// | x_{-2} | x_{-1} | x_0 | x_1 | x_2 |
///                   ^
///                   |
///                   return value
pub fn interp_weno5(
    phim2: f64,
    phim1: f64,
    phi: f64,
    phip1: f64,
    phip2: f64,
) -> f64 {
    let p0 = (1.0 / 3.0) * phi + (5.0 / 6.0) * phim1 - (1.0 / 6.0) * phim2;
    let p1 = (-1.0 / 6.0) * phip1 + (5.0 / 6.0) * phi + (1.0 / 3.0) * phim1;
    let p2 = (1.0 / 3.0) * phip2 - (7.0 / 6.0) * phip1 + (11.0 / 6.0) * phi;

    let beta0 = 13.0 / 12.0 * (phi - 2.0 * phim1 + phim2).powi(2)
        + 0.25 * (3.0 * phi - 4.0 * phim1 + phim2).powi(2);
    let beta1 = 13.0 / 12.0 * (phip1 - 2.0 * phi + phim1).powi(2)
        + 0.25 * (phip1 - phim1).powi(2);
    let beta2 = 13.0 / 12.0 * (phip2 - 2.0 * phip1 + phi).powi(2)
        + 0.25 * (phip2 - 4.0 * phip1 + 3.0 * phi).powi(2);

    let alpha0 = 0.3 / (beta0 + 1e-10).powi(2);
    let alpha1 = 0.6 / (beta1 + 1e-10).powi(2);
    let alpha2 = 0.1 / (beta2 + 1e-10).powi(2);

    (alpha0 * p0 + alpha1 * p1 + alpha2 * p2) / (alpha0 + alpha1 + alpha2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interp_weno5() {
        let phim2 = 1.0;
        let phim1 = 2.0;
        let phi = 3.0;
        let phip1 = 4.0;
        let phip2 = 5.0;
        let result = interp_weno5(phim2, phim1, phi, phip1, phip2);
        let expected_result = 2.5000000000000004;
        approx::assert_abs_diff_eq!(
            result,
            expected_result,
            epsilon = f64::EPSILON
        );
    }
}
