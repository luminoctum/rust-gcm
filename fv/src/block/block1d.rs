//! Block1D module
//! Block1D is a 1D data block with nvar variables and nghost ghost cells
//! The data is stored in a 1D vector with size nvar * (dim1 + 2 * nghost)
//! The data is stored in row-major order

use crate::block::iterator1d::{Iterator1D, Iterator1DMut};

pub struct Block1D<T> {
    pub nvar: usize,
    pub len1: usize,
    pub nghost: usize,
    pub data: Vec<T>,
}

impl<T: Default + Copy> Block1D<T> {
    pub fn new(nvar: usize, dim1: usize, nghost: usize) -> Self {
        let len1 = dim1 + 2 * nghost;
        let data = vec![T::default(); nvar * len1];
        Self {
            nvar,
            len1,
            nghost,
            data,
        }
    }

    pub fn size(&self) -> usize {
        self.data.len()
    }

    pub fn shape(&self) -> (usize, usize) {
        (self.nvar, self.len1)
    }

    pub fn icomp(&self, n: usize) -> usize {
        n * self.len1
    }

    pub fn interior(&self) -> Iterator1D<T> {
        Iterator1D {
            data: &self.data,
            current: self.nghost,
            end: self.len1 - self.nghost,
        }
    }

    pub fn interior_mut(&mut self) -> Iterator1DMut<T> {
        Iterator1DMut {
            data: &mut self.data,
            current: self.nghost,
            end: self.len1 - self.nghost,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::common::set_comp;

    #[test]
    fn test_data_block_1d_new() {
        let block = Block1D::<f32>::new(3, 4, 1);
        assert_eq!(block.nvar, 3);
        assert_eq!(block.len1, 6);
        assert_eq!(block.nghost, 1);
        assert_eq!(block.data.len(), 18);
    }

    #[test]
    fn test_data_block_1d_size() {
        let block = Block1D::<f32>::new(3, 4, 1);
        assert_eq!(block.size(), 18);
    }

    #[test]
    fn test_data_block_1d_shape() {
        let block = Block1D::<f32>::new(3, 4, 1);
        assert_eq!(block.shape(), (3, 6));
    }

    #[test]
    fn test_data_block_1d_interior() {
        let mut block = Block1D::new(3, 4, 1);

        // fill data
        for i in 0..block.size() {
            block.data[i] = i;
        }

        let interior = block.interior();

        // variable indices (offsets)
        let iv0 = block.icomp(0);
        let iv1 = block.icomp(1);
        let iv2 = block.icomp(2);

        // compare results
        for (i, var) in interior.enumerate() {
            assert_eq!(var[iv0], i + 1);
            assert_eq!(var[iv1], i + 1 + 6);
            assert_eq!(var[iv2], i + 1 + 12);
        }
    }

    #[test]
    fn test_interior_mut() {
        let mut block = Block1D::new(2, 4, 1);

        // Fill the block with some test data
        block.data = vec![0, 1, 2, 3, 4, 0, 0, 7, 8, 9, 10, 0];

        let iv1 = block.icomp(1);

        // Test that the iterator produces the expected values
        for x in block.interior_mut() {
            unsafe {
                set_comp(x, 0, 1);
                set_comp(x, iv1, 2);
            }
        }

        assert_eq!(block.data, vec![0, 1, 1, 1, 1, 0, 0, 2, 2, 2, 2, 0]);
    }
}
