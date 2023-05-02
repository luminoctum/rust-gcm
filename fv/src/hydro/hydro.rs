//! Hydrodynamics module
//!
//! Hydrodynamics provides the following functions:
//! 1. Reconstruct the left and right states (hydro_reconstruct.rs)
//! 2. Solve the Riemann problem (hydro_riemann.rs)
//! 3. Add flux divergence to the conserved variables (hydro_flux_divergence.rs)
//!
//! Hydrodynamics will modify the conserved variables in the EquationOfState

use itertools::izip;
use crate::block::block2d::Block2D;
use crate::eos::eos::EquationOfState;
use crate::riemann::roe_shallow_water::roe_shallow_water;
use crate::utils::defs::{Real, DIMENSION, NGHOST, NHYDRO};

pub struct Hydro {
    // components
    pub comps: [usize; NHYDRO],

    // left and right states
    pub wls: [Block2D<Real>; DIMENSION],
    pub wrs: [Block2D<Real>; DIMENSION],

    // fluxes
    pub flx: [Block2D<Real>; DIMENSION],
}

impl Hydro {
    pub const X1DIR: usize = 0;
    pub const X2DIR: usize = 1;
    pub const X3DIR: usize = 2;

    pub fn new(dim2: usize, dim1: usize) -> Self {
        let wls = [
            Block2D::new(NHYDRO, dim2, dim1, NGHOST),
            Block2D::new(NHYDRO, dim2, dim1, NGHOST),
        ];

        let wrs = [
            Block2D::new(NHYDRO, dim2, dim1, NGHOST),
            Block2D::new(NHYDRO, dim2, dim1, NGHOST),
        ];

        let flx = [
            Block2D::new(NHYDRO, dim2, dim1, NGHOST),
            Block2D::new(NHYDRO, dim2, dim1, NGHOST),
        ];

        // let comps: Vec<_> = (0..NHYDRO).map(|x| wls[0].icomp(x)).collect();
        let mut comps = [0; NHYDRO];
        for (i, comp) in comps.iter_mut().enumerate().take(NHYDRO) {
            *comp = wls[0].icomp(i);
        }

        Self {
            comps,
            wls,
            wrs,
            flx,
        }
    }

    /// # Safety
    ///
    /// Riemann solver for the x1 direction
    /// |  w_{-1} *|*  w_{0}  |
    ///           ^ ^
    ///           | |
    ///       wl(i) wr(i)
    pub unsafe fn riemann_solver_x1(&mut self) {
        let wl = self.wls[Hydro::X1DIR].interior_f1();
        let wr = self.wrs[Hydro::X1DIR].interior_f1();
        let flx = self.flx[Hydro::X1DIR].interior_f1_mut();
        let pos = [0];

        for (flx, wl, wr) in izip!(flx, wl, wr) {
            unsafe {
                roe_shallow_water(flx, wl, wr, Hydro::X1DIR, &self.comps, &pos);
            }
        }
    }

    /// Riemann solver for the x2 direction
    /// ------
    /// w_{-1}
    ///   * <----- wl(i)
    /// ------
    ///   * <----- wr(i)
    /// w_{0}
    /// ------
    pub fn riemann_solver_x2(&mut self) {
        unimplemented!();
        // WIP pseudocode
        // let wl = wls[X2DIR].interior_f2();
        // let wr = wrs[X2DIR].interior_f2();
        // let flx = flx[X2DIR].interior_f2_mut();
        //
        // for (flx, wl, wr) in izip!(flx, wl, wr)
        // {
        // roe_shallow_water(flx, wl, wr, X2DIR);
        // }
    }

    /// Add flux divergence to the conserved variables
    /// |  u_{0}  |
    /// ^         ^
    /// |         |
    /// flx(i)    flx(i+1)
    pub fn add_flux_divergence(&self, _eos: &mut EquationOfState, _dt: Real) {
        unimplemented!();
        // WIP pseudocode
        // let mut u = eos.u.interior_mut();
        // let flx = flx[ivx].interior();
        // let flx_p1 = flx[ivx].interior().offset_x1(1);
        // let area = coord.area1(0);
        // let area_p1 = coord.area(1);
        // let vol = coord.vol();
        //
        // for (mut u, flx, flx_p1, area, area_p1, vol) in
        // izip!(w, flx, flx_p1, area, area_p1, vol)
        // {
        // for n in ns..ne {
        // w[n] -= dt * (flx_p1[n]*area_p1 - flx[n]*area) / vol;
        // }
        // }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::defs::{NGHOST, NHYDRO};

    #[test]
    fn hydro_new() {
        let dim2 = 5;
        let dim1 = 7;
        let hydro = Hydro::new(dim2, dim1);

        for wl in hydro.wls.iter() {
            assert_eq!(wl.nvar, NHYDRO);
            assert_eq!(wl.len1, dim1 + 2 * NGHOST);
            assert_eq!(wl.len2, dim2 + 2 * NGHOST);
            assert_eq!(wl.nghost, NGHOST);
        }

        for wr in hydro.wrs.iter() {
            assert_eq!(wr.nvar, NHYDRO);
            assert_eq!(wr.len1, dim1 + 2 * NGHOST);
            assert_eq!(wr.len2, dim2 + 2 * NGHOST);
            assert_eq!(wr.nghost, NGHOST);
        }
    }
}
