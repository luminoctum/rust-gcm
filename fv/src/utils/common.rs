//! Common functions for the hydrodynamics module

/// # Safety
///
/// Not recommended for use. Used in test environment only.
pub unsafe fn set_comp<T: Copy>(var: *mut T, n: usize, val: T) {
    unsafe {
        *var.add(n) = val;
    }
}

/// # Safety
///
/// Not recommended for use. Used in test environment only.
pub unsafe fn add_comp<T: Copy + std::ops::AddAssign>(
    var: *mut T,
    n: usize,
    val: T,
) {
    unsafe {
        *var.add(n) += val;
    }
}
