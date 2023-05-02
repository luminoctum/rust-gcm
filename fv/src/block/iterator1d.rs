/// ! Iterator over the interior of a Block1D.
/// ! Iterator replaces the for loops in the Block1D methods

pub struct Iterator1D<'a, T> {
    pub data: &'a Vec<T>,
    pub current: usize,
    pub end: usize,
}

impl<'a, T> Iterator for Iterator1D<'a, T> {
    type Item = &'a [T];

    fn next(&mut self) -> Option<Self::Item> {
        if self.current < self.end {
            let index = self.current;
            self.current += 1;
            Some(&self.data[index..self.data.len()])
        } else {
            None
        }
    }
}

pub struct Iterator1DMut<'a, T> {
    pub data: &'a mut Vec<T>,
    pub current: usize,
    pub end: usize,
}

impl<'a, T> Iterator for Iterator1DMut<'a, T> {
    type Item = *mut T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current < self.end {
            let index = self.current;
            self.current += 1;
            Some(&mut self.data[index] as *mut T)
        } else {
            None
        }
    }
}
