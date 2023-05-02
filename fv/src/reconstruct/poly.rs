//! Polynomial Interpolation
//! Input arguments are cell averaged values
//! \todo use precomputed coefficients

/// | x_{-1} | x_0 |
///          ^
///          |
///          return value
pub fn interp_cp2(phim1: f64, phi: f64) -> f64 {
    0.5 * (phim1 + phi)
}

/// | x_{-1} | x_0 | x_1 |
///          ^
///          |
///          return value
pub fn interp_cp3(phim1: f64, phi: f64, phip1: f64) -> f64 {
    1.0 / 6.0 * (2. * phim1 + 5.0 * phi - 1.0 * phip1)
}

/// | x_{-2} | x_{-1} | x_0 | x_1 |
///                   ^
///                   |
///                   return value
pub fn interp_cp4(phim2: f64, phim1: f64, phi: f64, phip1: f64) -> f64 {
    -1.0 / 12.0 * (phim2 - 7.0 * phim1 - 7.0 * phi + phip1)
}

#[cfg(test)]
mod tests {
    use super::*;
    extern crate approx;

    #[test]
    fn test_interp_cp2() {
        let phim1 = 1.0;
        let phi = 2.0;
        let result = interp_cp2(phim1, phi);
        let expected_result = 1.5;
        approx::assert_abs_diff_eq!(
            result,
            expected_result,
            epsilon = f64::EPSILON
        );
    }

    #[test]
    fn test_interp_cp3() {
        let phim1 = 1.0;
        let phi = 2.0;
        let phip1 = 3.0;
        let result = interp_cp3(phim1, phi, phip1);
        let expected_result = 1.5;
        approx::assert_abs_diff_eq!(
            result,
            expected_result,
            epsilon = f64::EPSILON
        );
    }

    #[test]
    fn test_interp_cp4() {
        let phim2 = 0.0;
        let phim1 = 1.0;
        let phi = 2.0;
        let phip1 = 3.0;
        let result = interp_cp4(phim2, phim1, phi, phip1);
        let expected_result = 1.5;
        approx::assert_abs_diff_eq!(
            result,
            expected_result,
            epsilon = f64::EPSILON
        );
    }
}
