/// ! Iterator over the interior of a Block2D.
/// ! Iterator replaces the for loops in the Block2D methods.

/// An iterator over the interior of a Block2D.
pub struct Iterator2D<'a, T> {
    pub data: &'a Vec<T>,
    pub current1: usize,
    pub current2: usize,
    pub start1: usize,
    pub len1: usize,
    pub end1: usize,
    pub end2: usize,
}

impl<'a, T> Iterator for Iterator2D<'a, T> {
    type Item = &'a [T];

    fn next(&mut self) -> Option<Self::Item> {
        if self.current1 >= self.end1 {
            self.current1 = self.start1;
            self.current2 += 1;
        }

        if self.current2 >= self.end2 {
            return None;
        }

        let index = self.len1 * self.current2 + self.current1;

        self.current1 += 1;

        Some(&self.data[index..self.data.len()])
    }
}


/// An iterator over the mutable interior of a Block2D.
pub struct Iterator2DMut<'a, T> {
    pub data: &'a mut Vec<T>,
    pub len1: &'a usize,
    pub current1: usize,
    pub current2: usize,
    pub start1: usize,
    pub end1: usize,
    pub end2: usize,
}

impl<'a, T> Iterator for Iterator2DMut<'a, T> {
    type Item = *mut T;

    fn next(&mut self) -> Option<Self::Item> {
        let len1 = *self.len1;

        if self.current1 >= self.end1 {
            self.current1 = self.start1;
            self.current2 += 1;
        }

        if self.current2 >= self.end2 {
            return None;
        }

        let index = len1 * self.current2 + self.current1;

        self.current1 += 1;

        Some(&mut self.data[index] as *mut T)
    }
}
