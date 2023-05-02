/// ! Roe shallow water Riemann Solver
use crate::utils::{
    common::{add_comp, set_comp},
    defs::{Real, NHYDRO, X1DIR, X2DIR},
};

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers.
pub unsafe fn roe_shallow_water(
    flx: *mut Real,
    wli: &[Real],
    wri: &[Real],
    dir: usize,
    comps: &[usize; NHYDRO],
    _pos: &[usize],
) {
    let idn = comps[0];
    let iv1 = comps[1];
    let iv2 = comps[2];

    let (ivx, ivy) = match dir {
        X1DIR => (iv1, iv2),
        X2DIR => (iv2, iv1),
        _ => panic!("Invalid direction"),
    };

    let mut wave: [[Real; 3]; 3] = [[0.0; 3]; 3];
    let mut speed: [Real; 3] = [0.0; 3];

    // These are syntax for running on GPU, KEEP IT COMMENTED
    // let idx = thread::index_1d();
    // flx, wli, wri = shift_position!(pos, i, flx, wli, wri)

    let ubar = (wli[ivx] * wli[idn].sqrt() + wri[ivx] * wri[idn].sqrt())
        / (wli[idn].sqrt() + wri[idn].sqrt());
    let vbar = (wli[ivy] * wli[idn].sqrt() + wri[ivy] * wri[idn].sqrt())
        / (wli[idn].sqrt() + wri[idn].sqrt());
    let cbar = (0.5 * (wli[idn] + wri[idn])).sqrt();

    let delh = wri[idn] - wli[idn];
    let delu = wri[ivx] - wli[ivx];
    let delv = wri[ivy] - wli[ivy];
    let hbar = (wli[idn] * wri[idn]).sqrt();

    let a1 = 0.5 * (cbar * delh - hbar * delu) / cbar;
    let a2 = hbar * delv;
    let a3 = 0.5 * (cbar * delh + hbar * delu) / cbar;

    wave[0][0] = a1;
    wave[0][1] = a1 * (ubar - cbar);
    wave[0][2] = a1 * vbar;
    wave[1][0] = 0.0;
    wave[1][1] = 0.0;
    wave[1][2] = a2;
    wave[2][0] = a3;
    wave[2][1] = a3 * (ubar + cbar);
    wave[2][2] = a3 * vbar;

    speed[0] = (ubar - cbar).abs();
    speed[1] = ubar.abs();
    speed[2] = (ubar + cbar).abs();

    set_comp(flx, idn, 0.5 * (wli[idn] * wli[ivx] + wri[idn] * wri[ivx]));
    set_comp(
        flx,
        ivx,
        0.5 * (wli[idn] * wli[ivx].powi(2)
            + 0.5 * wli[idn].powi(2)
            + wri[idn] * wri[ivx].powi(2)
            + 0.5 * wri[idn].powi(2)),
    );
    set_comp(
        flx,
        ivy,
        0.5 * (wli[idn] * wli[ivx] * wli[ivy]
            + wri[idn] * wri[ivx] * wri[ivy]),
    );

    for r in 0..3 {
        add_comp(flx, idn, -0.5 * speed[r] * wave[r][0]);
        add_comp(flx, ivx, -0.5 * speed[r] * wave[r][1]);
        add_comp(flx, ivy, -0.5 * speed[r] * wave[r][2]);
    }
}
