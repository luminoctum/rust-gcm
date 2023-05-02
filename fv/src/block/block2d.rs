//! Block2D module
//! Block2D is a 2D data block with nvar variables, dim2 rows, dim1 columns,
//! and nghost ghost cells. Block2D stores data in a 1D array of size
//! nvar*dim2*dim1. The data is stored in row-major order (C-style).

use crate::block::iterator2d::{Iterator2D, Iterator2DMut};

/// A 2D data block with nvar variables, dim2 rows, dim1 columns, and nghost
/// ghost cells.
pub struct Block2D<T> {
    pub nvar: usize,
    pub len2: usize,
    pub len1: usize,
    pub len12: usize,
    pub nghost: usize,
    pub data: Vec<T>,
}

impl<T: Default + Copy> Block2D<T> {
    pub fn new(nvar: usize, dim2: usize, dim1: usize, nghost: usize) -> Self {
        let len2 = dim2 + 2 * nghost;
        let len1 = dim1 + 2 * nghost;
        let len12 = len1 * len2;
        let data = vec![T::default(); nvar * len12];
        Self {
            nvar,
            len2,
            len1,
            len12,
            nghost,
            data,
        }
    }

    pub fn size(&self) -> usize {
        self.data.len()
    }

    pub fn shape(&self) -> (usize, usize, usize) {
        (self.nvar, self.len2, self.len1)
    }

    pub fn icomp(&self, n: usize) -> usize {
        n * self.len12
    }

    /// |o|o|o|o|o|o|o|o|
    /// |o|x|x|x|x|x|x|o|
    /// |o|x|x|x|x|x|x|o|
    /// |o|x|x|x|x|x|x|o|
    /// |o|x|x|x|x|x|x|o|
    /// |o|o|o|o|o|o|o|o|
    /// This iterator loops over x (interior cells)
    pub fn interior(&self) -> Iterator2D<T> {
        Iterator2D {
            data: &self.data,
            current1: self.nghost,
            current2: self.nghost,
            start1: self.nghost,
            len1: self.len1,
            end1: self.len1 - self.nghost,
            end2: self.len2 - self.nghost,
        }
    }

    pub fn interior_mut(&mut self) -> Iterator2DMut<T> {
        Iterator2DMut {
            data: &mut self.data,
            len1: &self.len1,
            current1: self.nghost,
            current2: self.nghost,
            start1: self.nghost,
            end1: self.len1 - self.nghost,
            end2: self.len2 - self.nghost,
        }
    }

    /// Expand the interior region in the x1 direction by 1
    /// and shift by an offset. Used in reconstruction.
    pub fn interior_x1(&self, offset: i32) -> Iterator2D<T> {
        let start = (self.nghost as i32 + offset - 1) as usize;
        let end = (self.len1 as i32 - self.nghost as i32 + offset + 1) as usize;

        Iterator2D {
            data: &self.data,
            current1: start,
            current2: self.nghost,
            start1: start,
            len1: self.len1,
            end1: end,
            end2: self.len2 - self.nghost,
        }
    }

    pub fn interior_x1_mut(&mut self, offset: i32) -> Iterator2DMut<T> {
        let start = (self.nghost as i32 + offset - 1) as usize;
        let end = (self.len1 as i32 - self.nghost as i32 + offset + 1) as usize;

        Iterator2DMut {
            data: &mut self.data,
            len1: &self.len1,
            current1: start,
            current2: self.nghost,
            start1: start,
            end1: end,
            end2: self.len2 - self.nghost,
        }
    }

    /// |o|o|o|o|o|o|o|o|
    /// |o*x*x*x*x*x*x*o|
    /// |o*x*x*x*x*x*x*o|
    /// |o*x*x*x*x*x*x*o|
    /// |o*x*x*x*x*x*x*o|
    /// |o|o|o|o|o|o|o|o|
    /// This iterator loops over * (interior cell faces)
    pub fn interior_f1(&self) -> Iterator2D<T> {
        Iterator2D {
            data: &self.data,
            current1: self.nghost,
            current2: self.nghost,
            start1: self.nghost,
            len1: self.len1,
            end1: self.len1 - self.nghost + 1,
            end2: self.len2 - self.nghost,
        }
    }

    pub fn interior_f1_mut(&mut self) -> Iterator2DMut<T> {
        Iterator2DMut {
            data: &mut self.data,
            len1: &self.len1,
            current1: self.nghost,
            current2: self.nghost,
            start1: self.nghost,
            end1: self.len1 - self.nghost + 1,
            end2: self.len2 - self.nghost,
        }
    }

    pub fn interior_f2(&self) -> Iterator2D<T> {
        Iterator2D {
            data: &self.data,
            current1: self.nghost,
            current2: self.nghost,
            start1: self.nghost,
            len1: self.len1,
            end1: self.len1 - self.nghost,
            end2: self.len2 - self.nghost + 1,
        }
    }

    pub fn interior_f2_mut(&mut self) -> Iterator2DMut<T> {
        Iterator2DMut {
            data: &mut self.data,
            len1: &self.len1,
            current1: self.nghost,
            current2: self.nghost,
            start1: self.nghost,
            end1: self.len1 - self.nghost,
            end2: self.len2 - self.nghost + 1,
        }
    }

    /// Similar to interior_x1, but in the x2 direction
    /// Expand the interior region in the x2 direction by 1
    /// and shift by an offset. Used in reconstruction.
    pub fn interior_x2(&self, offset: i32) -> Iterator2D<T> {
        let start = (self.nghost as i32 + offset - 1) as usize;
        let end = (self.len2 as i32 - self.nghost as i32 + offset + 1) as usize;

        Iterator2D {
            data: &self.data,
            current1: self.nghost,
            current2: start,
            start1: self.nghost,
            len1: self.len1,
            end1: self.len1 - self.nghost,
            end2: end,
        }
    }

    pub fn interior_x2_mut(&mut self, offset: i32) -> Iterator2DMut<T> {
        let start = (self.nghost as i32 + offset - 1) as usize;
        let end = (self.len2 as i32 - self.nghost as i32 + offset + 1) as usize;

        Iterator2DMut {
            data: &mut self.data,
            len1: &self.len1,
            current1: self.nghost,
            current2: start,
            start1: self.nghost,
            end1: self.len1 - self.nghost,
            end2: end,
        }
    }

    pub fn all(&self) -> Iterator2D<T> {
        Iterator2D {
            data: &self.data,
            current1: 0,
            current2: 0,
            start1: 0,
            len1: self.len1,
            end1: self.len1,
            end2: self.len2,
        }
    }

    pub fn all_mut(&mut self) -> Iterator2DMut<T> {
        Iterator2DMut {
            data: &mut self.data,
            len1: &self.len1,
            current1: 0,
            current2: 0,
            start1: 0,
            end1: self.len1,
            end2: self.len2,
        }
    }

    pub fn at(&self, j: i32, i: i32) -> &[T] {
        let j1 = (self.nghost as i32 + j) as usize;
        let i1 = (self.nghost as i32 + i) as usize;

        let index = j1 * self.len1 + i1;
        &self.data[index..self.data.len()]
    }

    pub fn get(&self, n: usize, j: i32, i: i32) -> T {
        let j1 = (self.nghost as i32 + j) as usize;
        let i1 = (self.nghost as i32 + i) as usize;

        let index = n * self.len12 + j1 * self.len1 + i1;
        unsafe { *self.data.get_unchecked(index) }
    }

    pub fn set(&mut self, n: usize, j: i32, i: i32) -> &mut T {
        let j1: usize = (self.nghost as i32 + j).try_into().unwrap();
        let i1: usize = (self.nghost as i32 + i).try_into().unwrap();

        let index = n * self.len12 + j1 * self.len1 + i1;
        unsafe { self.data.get_unchecked_mut(index) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::common::set_comp;

    #[test]
    fn test_data_block_2d_new() {
        let nvar = 4;
        let dim2 = 10;
        let dim1 = 20;
        let nghost = 2;
        let w = Block2D::<f64>::new(nvar, dim2, dim1, nghost);

        assert_eq!(w.nvar, nvar);
        assert_eq!(w.len2, dim2 + 2 * nghost);
        assert_eq!(w.len1, dim1 + 2 * nghost);
        assert_eq!(w.nghost, nghost);
    }

    #[test]
    fn test_data_block_2d_get() {
        let nvar = 4;
        let dim2 = 10;
        let dim1 = 20;
        let nghost = 2;
        let mut w = Block2D::<f64>::new(nvar, dim2, dim1, nghost);

        // Fill the w with some values.
        for i in 0..w.size() {
            w.data[i] = i as f64;
        }

        approx::assert_abs_diff_eq!(w.get(0, 0, 0), 50.0);
    }

    #[test]
    fn test_data_block_2d_set() {
        let nvar = 4;
        let dim2 = 10;
        let dim1 = 20;
        let nghost = 2;
        let mut w = Block2D::<f64>::new(nvar, dim2, dim1, nghost);

        // Fill the w with some values.
        for i in 0..w.size() {
            w.data[i] = i as f64;
        }

        *w.set(0, 0, 0) = 100.0;

        approx::assert_abs_diff_eq!(w.get(0, 0, 0), 100.0);
        approx::assert_abs_diff_eq!(w.data[50], 100.0);
    }

    #[test]
    fn test_data_block_2d_interior() {
        let nvar = 2;
        let dim2 = 3;
        let dim1 = 3;
        let nghost = 1;
        let mut w = Block2D::new(nvar, dim2, dim1, nghost);

        // Fill the data_block with some values.
        for i in 0..w.size() {
            w.data[i] = i;
        }

        // Collect the interior points using the iterator.
        let interior_points: Vec<_> = w.interior().collect();

        // variable component index
        let iv0 = w.icomp(0);
        let iv1 = w.icomp(1);

        // Compare the actual and expected interior points.
        assert_eq!(interior_points.len(), 9);
        assert_eq!(interior_points[0][iv0], 6);
        assert_eq!(interior_points[0][iv1], 31);
        assert_eq!(interior_points[1][iv0], 7);
        assert_eq!(interior_points[1][iv1], 32);
    }

    #[test]
    fn test_data_block_2d_interior_mut() {
        let nvar = 2;
        let dim2 = 3;
        let dim1 = 3;
        let nghost = 1;
        let mut w = Block2D::new(nvar, dim2, dim1, nghost);

        // Fill the data_block with some values.
        for i in 0..w.size() {
            w.data[i] = i;
        }

        // Collect the interior points using the iterator.
        let interior_points: Vec<_> = w.interior_mut().collect();

        unsafe {
            set_comp(interior_points[0], 0, 100);
            assert_eq!(*interior_points[0], 100);
        }
    }

    #[test]
    fn test_data_block_2d_all() {
        let nvar = 2;
        let dim2 = 3;
        let dim1 = 3;
        let nghost = 1;
        let mut w = Block2D::new(nvar, dim2, dim1, nghost);

        // Fill the data_block with some values.
        for i in 0..w.size() {
            w.data[i] = i;
        }

        // Collect all points using the iterator.
        let all_points: Vec<_> = w.all().collect();

        // variable component index
        let iv0 = w.icomp(0);
        let iv1 = w.icomp(1);

        // Compare the actual and expected interior points.
        assert_eq!(all_points.len(), 25);
        assert_eq!(all_points[0][iv0], 0);
        assert_eq!(all_points[0][iv1], 25);
        assert_eq!(all_points[1][iv0], 1);
        assert_eq!(all_points[1][iv1], 26);
    }

    #[test]
    fn test_data_block_2d_all_mut() {
        let nvar = 2;
        let dim2 = 3;
        let dim1 = 3;
        let nghost = 1;
        let mut w = Block2D::new(nvar, dim2, dim1, nghost);

        // Fill the data_block with some values.
        for i in 0..w.size() {
            w.data[i] = i;
        }

        // Collect all points using the iterator.
        let all_points: Vec<_> = w.all_mut().collect();

        unsafe {
            set_comp(all_points[0], 0, 100);
            assert_eq!(*all_points[0], 100);
        }
    }

    #[test]
    fn test_interior_x1_iterator() {
        let nvar = 2;
        let dim2 = 3;
        let dim1 = 3;
        let nghost = 1;
        let mut data_block = Block2D::new(nvar, dim2, dim1, nghost);

        // Fill the data_block with some values.
        for i in 0..data_block.size() {
            data_block.data[i] = i;
        }

        // Collect the expanded interior points using the iterator.
        let interior_points_x1: Vec<_> =
            data_block.interior_x1(0).map(|x| x[0]).collect();

        // Define the expected interior points with an offset.
        let expected_interior_points_x1 =
            vec![5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19];

        // Compare the actual and expected interior points with an offset.
        assert_eq!(interior_points_x1, expected_interior_points_x1);
    }

    #[test]
    fn test_interior_x2_iterator() {
        let nvar = 2;
        let dim2 = 3;
        let dim1 = 3;
        let nghost = 1;
        let mut data_block = Block2D::new(nvar, dim2, dim1, nghost);

        // Fill the data_block with some values.
        for i in 0..data_block.size() {
            data_block.data[i] = i;
        }

        // Collect the expanded interior points using the iterator.
        let iv1 = data_block.icomp(1);
        let interior_points_x2: Vec<_> =
            data_block.interior_x2(0).map(|x| x[iv1]).collect();

        // Define the expected interior points with an offset.
        let expected_interior_points_x2 =
            vec![26, 27, 28, 31, 32, 33, 36, 37, 38, 41, 42, 43, 46, 47, 48];

        // Compare the actual and expected interior points with an offset.
        assert_eq!(interior_points_x2, expected_interior_points_x2);
    }
}
