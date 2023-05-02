//! Equation of state module
//!
//! Equation of state provides the following functions:
//! 1. Convert conserved variables to primitive variables
//! 2. Convert primitive variables to conserved variables

use crate::block::block2d::Block2D;
use crate::utils::{
    common::set_comp,
    defs::{Real, NGHOST, NHYDRO},
};

pub struct EquationOfState {
    // components
    pub comps: [usize; NHYDRO],

    // primary variables
    pub w: Block2D<Real>,

    // conserved variables
    pub u: Block2D<Real>,

    // conserved variable register
    _u1: Block2D<Real>,
    _u2: Block2D<Real>,
}

impl EquationOfState {
    pub fn new(dim2: usize, dim1: usize) -> Self {
        let w = Block2D::new(NHYDRO, dim2, dim1, NGHOST);
        let u = Block2D::new(NHYDRO, dim2, dim1, NGHOST);

        let u1 = Block2D::new(NHYDRO, dim2, dim1, NGHOST);
        let u2 = Block2D::new(NHYDRO, dim2, dim1, NGHOST);

        let mut comps = [0; NHYDRO];
        for (i, comp) in comps.iter_mut().enumerate().take(NHYDRO) {
            *comp = w.icomp(i);
        }

        Self {
            comps,
            w,
            u,
            _u1: u1,
            _u2: u2,
        }
    }

    /// # Safety
    ///
    /// This function may be unsafe because it dereferences a raw pointer.
    pub fn conserved_to_primitive(&mut self) {
        let w = self.w.all_mut();
        let u = self.u.all();
        let gm1 = 0.4;

        let [idn, iv1, iv2, iv3, ipr] = self.comps;

        for (w, u) in w.zip(u) {
            let rho = u[idn];
            let vx = u[iv1] / rho;
            let vy = u[iv2] / rho;
            let vz = u[iv3] / rho;
            let ke = 0.5 * rho * (vx * vx + vy * vy + vz * vz);

            unsafe {
                set_comp(w, idn, rho);
                set_comp(w, iv1, vx);
                set_comp(w, iv2, vy);
                set_comp(w, iv3, vz);
                set_comp(w, ipr, gm1 * (u[ipr] - ke));
            }
        }
    }

    pub fn cost_conserved_to_primitive(&self) -> Real {
        let overhead = 1.0;
        (self.w.size() as Real) * (self.u.size() as Real) * overhead
    }

    /// # Safety
    ///
    /// This function may be unsafe because it dereferences a raw pointer.
    pub fn primitive_to_conserved(&mut self) {
        let w = self.w.all();
        let u = self.u.all_mut();
        let gm1 = 0.4;

        let idn = self.w.icomp(0);
        let iv1 = self.w.icomp(1);
        let iv2 = self.w.icomp(2);
        let iv3 = self.w.icomp(3);
        let ipr = self.w.icomp(4);

        for (w, u) in w.zip(u) {
            let rho = w[idn];
            let vx = w[iv1];
            let vy = w[iv2];
            let vz = w[iv3];
            let ke = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
            let ie = w[ipr] / gm1;

            unsafe {
                set_comp(u, idn, rho);
                set_comp(u, iv1, rho * vx);
                set_comp(u, iv2, rho * vy);
                set_comp(u, iv3, rho * vz);
                set_comp(u, ipr, ke + ie);
            }
        }
    }

    pub fn cost_primitive_to_conserved(&self) -> Real {
        let overhead = 1.0;
        (self.w.size() as Real) * (self.u.size() as Real) * overhead
    }
}

impl Default for EquationOfState {
    fn default() -> Self {
        let w = Block2D::new(0, 0, 0, 0);
        let u = Block2D::new(0, 0, 0, 0);

        let u1 = Block2D::new(0, 0, 0, 0);
        let u2 = Block2D::new(0, 0, 0, 0);

        let comps = [0; NHYDRO];
        Self {
            comps,
            w,
            u,
            _u1: u1,
            _u2: u2,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let eos = EquationOfState::new(20, 10);

        assert_eq!(eos.w.shape(), (NHYDRO, 20 + 2 * NGHOST, 10 + 2 * NGHOST));
        assert_eq!(eos.u.shape(), (NHYDRO, 20 + 2 * NGHOST, 10 + 2 * NGHOST));
    }

    #[test]
    fn test_conserved_to_primitive() {
        let mut eos = EquationOfState::new(3, 3);
        eos.u.data.fill(1.0);

        eos.conserved_to_primitive();

        // variable indices
        let idn = eos.w.icomp(0);
        let iv1 = eos.w.icomp(1);
        let iv2 = eos.w.icomp(2);
        let iv3 = eos.w.icomp(3);
        let ipr = eos.w.icomp(4);

        for (w, u) in eos.w.all().zip(eos.u.all()) {
            let rho = w[idn];
            let vx = w[iv1];
            let vy = w[iv2];
            let vz = w[iv3];
            let ke = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
            let pr = gm1() * (u[ipr] - ke);

            approx::assert_abs_diff_eq!(w[idn], rho);
            approx::assert_abs_diff_eq!(w[iv1], vx);
            approx::assert_abs_diff_eq!(w[iv2], vy);
            approx::assert_abs_diff_eq!(w[iv3], vz);
            approx::assert_abs_diff_eq!(w[ipr], pr);
        }
    }

    #[test]
    fn test_primitive_to_conserved() {
        let mut eos = EquationOfState::new(3, 3);
        eos.w.data.fill(1.0);

        eos.primitive_to_conserved();

        // variable indices
        let idn = eos.w.icomp(0);
        let iv1 = eos.w.icomp(1);
        let iv2 = eos.w.icomp(2);
        let iv3 = eos.w.icomp(3);
        let ipr = eos.w.icomp(4);

        for (w, u) in eos.w.all().zip(eos.u.all()) {
            let rho = w[idn];
            let vx = w[iv1];
            let vy = w[iv2];
            let vz = w[iv3];
            let ke = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
            let ie = w[ipr] / gm1();

            approx::assert_abs_diff_eq!(u[idn], rho);
            approx::assert_abs_diff_eq!(u[iv1], rho * vx);
            approx::assert_abs_diff_eq!(u[iv2], rho * vy);
            approx::assert_abs_diff_eq!(u[iv3], rho * vz);
            approx::assert_abs_diff_eq!(u[ipr], ke + ie);
        }
    }

    fn gm1() -> Real {
        0.4
    }
}
