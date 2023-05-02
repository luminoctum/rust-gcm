// Precision
pub type Real = f64;

// Spatial dimension of the problem
pub const DIMENSION: usize = 2;

// Number of ghost zones surrounding the active domain
pub const NGHOST: usize = 3;

// Number of hydrodynamic variables
pub const NHYDRO: usize = 5;
pub const NVAPOR: usize = 0;

// Direction constants
pub const X1DIR: usize = 0;
pub const X2DIR: usize = 1;
pub const X3DIR: usize = 2;

// \todo(CLI) KEEP THIS COMMENT, WE NEED A WAY TO SHOW THIS SYSYTEM
// system = "shallow_water"
