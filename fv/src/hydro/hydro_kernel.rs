#![cfg_attr(
    target_os = "cuda",
    no_std,
    feature(register_attr),
    register_attr(nvvm_internal)
)]

use cuda_std::prelude::*;

#[kernel]
pub unsafe fn interp_weno3_gpu(
    wl: *mut Real,
    wr: *mut Real,
    wm1: &[Real],
    w: &[Real],
    wp1: &[Real],
    pos: usize,
    comps: &[usize],
) {
    let idx = thread::index_1d();

    // wl, wr, wm2, wm1, w, wp1, wp2 = shift_position!(pos, idx,
    // wl, wr, wm2, wm1, w, wp1, wp2)

    // let j = shift[i] + idx.y * chunk;
    for n in comps {
        set_comp(wl, n, interp_weno3(wp1, w, wm1));
        set_comp(wr, n, interp_weno3(wm1, w, wp1));
    }
}
