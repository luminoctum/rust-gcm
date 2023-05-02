//! Implements the hydrodynamics module of the code.
//!
//! This code provides the reconstruction of the left and right states
//! from the cell-centered primitive variables.

use itertools::izip;
use crate::eos::eos::EquationOfState;
use crate::reconstruct::{weno3::interp_weno3, weno5::interp_weno5};
use crate::utils::{
    common::set_comp,
};

use crate::hydro::hydro::Hydro;

impl Hydro {
    /// Reconstruct the left and right states in the x1 direction
    /// | w_{-2} | w_{-1} |* w_{0} *| w_{1} | w_{2} |
    ///                    ^       ^
    ///                    |       |
    ///                    wr(i)   wl(i+1)
    pub fn reconstruct_x1(&mut self, eos: &EquationOfState, order: usize) {
        let wl = self.wls[Self::X1DIR].interior_x1_mut(1);
        let wr = self.wrs[Self::X1DIR].interior_x1_mut(0);

        let wm2 = eos.w.interior_x1(-2);
        let wm1 = eos.w.interior_x1(-1);
        let w = eos.w.interior_x1(0);
        let wp1 = eos.w.interior_x1(1);
        let wp2 = eos.w.interior_x1(2);

        match order {
            1 => {
                for (wl, wr, w) in izip!(wl, wr, w) {
                    for n in self.comps {
                        unsafe {
                            set_comp(wl, n, w[n]);
                            set_comp(wr, n, w[n]);
                        }
                    }
                }
            },

            3 => {
                for (wl, wr, wm1, w, wp1) in izip!(wl, wr, wm1, w, wp1) {
                    for n in self.comps {
                        unsafe {
                            set_comp(wl, n, interp_weno3(wp1[n], w[n], wm1[n]));
                            set_comp(wr, n, interp_weno3(wm1[n], w[n], wp1[n]));
                        }
                    }
                }
            },

            5 => {
                for (wl, wr, wm2, wm1, w, wp1, wp2) in
                    izip!(wl, wr, wm2, wm1, w, wp1, wp2)
                {
                    for n in self.comps {
                        unsafe {
                            set_comp(
                                wl,
                                n,
                                interp_weno5(wp2[n], wp1[n], w[n], wm1[n], wm2[n]),
                            );
                            set_comp(
                                wr,
                                n,
                                interp_weno5(wm2[n], wm1[n], w[n], wp1[n], wp2[n]),
                            );
                        }
                    }
                }
            },

            _ => panic!("Invalid order"),
        }
    }

    /// Reconstruct the left and right states in the x2 direction
    /// ------
    /// w_{-2}
    /// ------
    /// w_{-1}
    /// ------ <- wr(i)
    /// w_{0}
    /// ------ <- wl(i+1)
    /// w_{1}
    /// ------
    /// w_{2}
    pub fn reconstruct_x2(&mut self, eos: &EquationOfState, order: usize) {
        let wl = self.wls[Self::X2DIR].interior_x2_mut(1);
        let wr = self.wrs[Self::X2DIR].interior_x2_mut(0);

        let wm2 = eos.w.interior_x2(-2);
        let wm1 = eos.w.interior_x2(-1);
        let w = eos.w.interior_x2(0);
        let wp1 = eos.w.interior_x2(1);
        let wp2 = eos.w.interior_x2(2);

        match order {
            1 => {
                for (wl, wr, w) in izip!(wl, wr, w) {
                    for n in self.comps {
                        unsafe {
                            set_comp(wl, n, w[n]);
                            set_comp(wr, n, w[n]);
                        }
                    }
                }
            },

            3 => {
                for (wl, wr, wm1, w, wp1) in izip!(wl, wr, wm1, w, wp1) {
                    for n in self.comps {
                        unsafe {
                            set_comp(wl, n, interp_weno3(wp1[n], w[n], wm1[n]));
                            set_comp(wr, n, interp_weno3(wm1[n], w[n], wp1[n]));
                        }
                    }
                }
            },

            5 => {
                for (wl, wr, wm2, wm1, w, wp1, wp2) in
                    izip!(wl, wr, wm2, wm1, w, wp1, wp2)
                {
                    for n in self.comps {
                        unsafe {
                            set_comp(
                                wl,
                                n,
                                interp_weno5(wp2[n], wp1[n], w[n], wm1[n], wm2[n]),
                            );
                            set_comp(
                                wr,
                                n,
                                interp_weno5(wm2[n], wm1[n], w[n], wp1[n], wp2[n]),
                            );
                        }
                    }
                }
            },

            _ => panic!("Invalid order"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eos::eos::EquationOfState;
    use crate::reconstruct::weno5::interp_weno5;
    use crate::utils::defs::{Real, NHYDRO};
    use crate::hydro::hydro::Hydro;

    #[test]
    fn reconstruct_x1() {
        let dim2 = 5;
        let dim1 = 7;
        let mut hydro = Hydro::new(dim2, dim1);
        let mut eos = EquationOfState::new(dim2, dim1);

        // Fill the eos with some values.
        for i in 0..eos.w.size() {
            eos.w.data[i] = i as Real;
        }

        let eos = eos;
        hydro.reconstruct_x1(&eos, 5);

        for n in 0..NHYDRO {
            for j in 0..dim2 as i32 {
                for i in 0..=dim1 as i32 {
                    let wl = hydro.wls[Hydro::X1DIR].get(n, j, i);
                    let wr = hydro.wrs[Hydro::X1DIR].get(n, j, i);

                    let wm3 = eos.w.get(n, j, i - 3);
                    let wm2 = eos.w.get(n, j, i - 2);
                    let wm1 = eos.w.get(n, j, i - 1);
                    let w = eos.w.get(n, j, i);
                    let wp1 = eos.w.get(n, j, i + 1);
                    let wp2 = eos.w.get(n, j, i + 2);

                    let expected_wl = interp_weno5(wp1, w, wm1, wm2, wm3);
                    let expected_wr = interp_weno5(wm2, wm1, w, wp1, wp2);

                    approx::assert_abs_diff_eq!(
                        wl,
                        expected_wl,
                        epsilon = f64::EPSILON
                    );

                    approx::assert_abs_diff_eq!(
                        wr,
                        expected_wr,
                        epsilon = f64::EPSILON
                    );
                }
            }
        }

        hydro.reconstruct_x1(&eos, 3);

        for n in 0..NHYDRO {
            for j in 0..dim2 as i32 {
                for i in 0..=dim1 as i32 {
                    let wl = hydro.wls[Hydro::X1DIR].get(n, j, i);
                    let wr = hydro.wrs[Hydro::X1DIR].get(n, j, i);

                    let wm2 = eos.w.get(n, j, i - 2);
                    let wm1 = eos.w.get(n, j, i - 1);
                    let w = eos.w.get(n, j, i);
                    let wp1 = eos.w.get(n, j, i + 1);

                    let expected_wl = interp_weno3(w, wm1, wm2);
                    let expected_wr = interp_weno3(wm1, w, wp1);

                    approx::assert_abs_diff_eq!(
                        wl,
                        expected_wl,
                        epsilon = f64::EPSILON
                    );

                    approx::assert_abs_diff_eq!(
                        wr,
                        expected_wr,
                        epsilon = f64::EPSILON
                    );
                }
            }
        }
    }

    #[test]
    fn reconstruct_x2() {
        let dim2 = 5;
        let dim1 = 7;
        let mut hydro = Hydro::new(dim2, dim1);
        let mut eos = EquationOfState::new(dim2, dim1);

        // Fill the eos with some values.
        for i in 0..eos.w.size() {
            eos.w.data[i] = i as Real;
        }

        let eos = eos;
        hydro.reconstruct_x2(&eos, 5);

        for n in 0..NHYDRO {
            for j in 0..=dim2 as i32 {
                for i in 0..dim1 as i32 {
                    let wl = hydro.wls[Hydro::X2DIR].get(n, j, i);
                    let wr = hydro.wrs[Hydro::X2DIR].get(n, j, i);

                    let wm3 = eos.w.get(n, j - 3, i);
                    let wm2 = eos.w.get(n, j - 2, i);
                    let wm1 = eos.w.get(n, j - 1, i);
                    let w = eos.w.get(n, j, i);
                    let wp1 = eos.w.get(n, j + 1, i);
                    let wp2 = eos.w.get(n, j + 2, i);

                    let expected_wl = interp_weno5(wp1, w, wm1, wm2, wm3);
                    let expected_wr = interp_weno5(wm2, wm1, w, wp1, wp2);

                    approx::assert_abs_diff_eq!(
                        wl,
                        expected_wl,
                        epsilon = f64::EPSILON
                    );

                    approx::assert_abs_diff_eq!(
                        wr,
                        expected_wr,
                        epsilon = f64::EPSILON
                    );
                }
            }
        }

        hydro.reconstruct_x2(&eos, 3);

        for n in 0..NHYDRO {
            for j in 0..=dim2 as i32 {
                for i in 0..dim1 as i32 {
                    let wl = hydro.wls[Hydro::X2DIR].get(n, j, i);
                    let wr = hydro.wrs[Hydro::X2DIR].get(n, j, i);

                    let wm2 = eos.w.get(n, j - 2, i);
                    let wm1 = eos.w.get(n, j - 1, i);
                    let w = eos.w.get(n, j, i);
                    let wp1 = eos.w.get(n, j + 1, i);

                    let expected_wl = interp_weno3(w, wm1, wm2);
                    let expected_wr = interp_weno3(wm1, w, wp1);

                    approx::assert_abs_diff_eq!(
                        wl,
                        expected_wl,
                        epsilon = f64::EPSILON
                    );

                    approx::assert_abs_diff_eq!(
                        wr,
                        expected_wr,
                        epsilon = f64::EPSILON
                    );
                }
            }
        }
    }
}
