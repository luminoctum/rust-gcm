//! LMARS Riemann solver
//! Reference: Chen et al. (2017)
//! \todo Does not work now. Needs fix

use crate::{defs, variable::Variable};

fn riemann_solver_lmars(
    flx: &mut Variable,
    wli: &Variable,
    wri: &Variable,
    ivx: usize,
) {
    let ivy = IVX + ((ivx - IVX) + 1) % 3;
    let ivz = IVX + ((ivx - IVX) + 2) % 3;

    // Assuming Thermodynamics struct and pmy_block are defined elsewhere
    // let pthermo = pmy_block.pimpl.pthermo;
    // let gamma = pmy_block.peos.get_gamma();

    // Assuming gamma is defined as a constant
    let gamma: Real = 1.4;

    // Correction for gamma
    // left
    let mut fsig = 1.0;
    let mut feps = 1.0;
    for n in 1..=NVAPOR {
        // Assuming pthermo.get_cv_ratio(n) and pthermo.get_mass_ratio(n) are
        // defined as constants
        let cv_ratio_n: Real = 1.0; // Replace with the correct value
        let mass_ratio_n: Real = 1.0; // Replace with the correct value
        fsig += wli[n] * (cv_ratio_n - 1.0);
        feps += wli[n] * (1.0 / mass_ratio_n - 1.0);
    }
    let kappal = 1.0 / (gamma - 1.0) * fsig / feps;

    // right
    fsig = 1.0;
    feps = 1.0;
    for n in 1..=NVAPOR {
        // Assuming pthermo.get_cv_ratio(n) and pthermo.get_mass_ratio(n) are
        // defined as constants
        let cv_ratio_n: Real = 1.0; // Replace with the correct value
        let mass_ratio_n: Real = 1.0; // Replace with the correct value
        fsig += wri[n] * (cv_ratio_n - 1.0);
        feps += wri[n] * (1.0 / mass_ratio_n - 1.0);
    }
    let kappar = 1.0 / (gamma - 1.0) * fsig / feps;

    // enthalpy
    let hl = wli[IPR] / wli[IDN] * (kappal + 1.0)
        + 0.5 * (wli[ivx].powi(2) + wli[ivy].powi(2) + wli[ivz].powi(2));
    let hr = wri[IPR] / wri[IDN] * (kappar + 1.0)
        + 0.5 * (wri[ivx].powi(2) + wri[ivy].powi(2) + wri[ivz].powi(2));

    let rhobar = 0.5 * (wli[IDN] + wri[IDN]);
    let cbar = (0.5
        * (1.0 + (1.0 / kappar + 1.0 / kappal) / 2.0)
        * (wli[IPR] + wri[IPR])
        / rhobar)
        .sqrt();
    let pbar = 0.5 * (wli[IPR] + wri[IPR])
        + 0.5 * (rhobar * cbar) * (wli[ivx] - wri[ivx]);
    let ubar = 0.5 * (wli[ivx] + wri[ivx])
        + 0.5 / (rhobar * cbar) * (wli[IPR] - wri[IPR]);

    let mut rd = 1.0;
    if ubar > 0.0 {
        // volume mixing ratio to mass mixing ratio
        for n in 1..=NVAPOR {
            rd -= wli[n];
        }

        flx[IDN] = ubar * wli[IDN] * rd;
        for n in 1..=NVAPOR {
            flx[n] = ubar * wli[IDN] * wli[n];
        }
        flx[ivx] = ubar * wli[IDN] * wli[ivx] + pbar;
        flx[ivy] = ubar * wli[IDN] * wli[ivy];
        flx[ivz] = ubar * wli[IDN] * wli[ivz];
        flx[IPR] = ubar * wli[IDN] * hl;
    } else {
        // volume mixing ratio to mass mixing ratio
        for n in 1..=NVAPOR {
            rd -= wri[n];
        }

        flx[IDN] = ubar * wri[IDN] * rd;
        for n in 1..=NVAPOR {
            flx[n] = ubar * wri[IDN] * wri[n];
        }
        flx[ivx] = ubar * wri[IDN] * wri[ivx] + pbar;
        flx[ivy] = ubar * wri[IDN] * wri[ivy];
        flx[ivz] = ubar * wri[IDN] * wri[ivz];
        flx[IPR] = ubar * wri[IDN] * hr;
    }
}
