//! Mesh module
//!
//! This module contains the MeshBlock and Mesh structs.
//! The MeshBlock struct contains all physics moldules for a single meshblock.
//! The MeshBlock struct is the main struct that manages all physics modules.
//! The Mesh struct contains a vector of MeshBlocks.
//! The Mesh struct is the main struct that manages communication between
//! MeshBlocks.
use crate::eos::eos::EquationOfState;

pub struct MeshBlock {
    nx3: usize,
    nx2: usize,
    nx1: usize,

    // TODO: Add other modules
    pub peos: Box<EquationOfState>,
    pub hydro: Box<Hydro>,
    // field: Field,
    // rad: Radiation,
    // tracer: Tracer,
    // chem: Chemistry,
}

impl Default for MeshBlock {
    fn default() -> Self {
        Self {
            nx3: 0,
            nx2: 0,
            nx1: 0,
            peos: Box::new(EquationOfState::default()),
        }
    }
}

impl MeshBlock {
    pub fn new(nx2: usize, nx1: usize) -> Self {
        let peos = Box::new(EquationOfState::new(nx2, nx1));
        Self {
            nx1,
            nx2,
            nx3,
            peos,
        }
    }

    pub fn eos(&mut self, eos: Box<EquationOfState>) -> &mut Self {
        self.peos = eos;
        self
    }
}

pub struct Mesh {
    meshblock: Vec<MeshBlock>,
}
