#[inline] pub fn idx(i: usize, j: usize, nx: usize) -> usize { j*nx + i }
#[inline] pub fn pid(i: isize, n: usize) -> usize { ((i).rem_euclid(n as isize)) as usize }
