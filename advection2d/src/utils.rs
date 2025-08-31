/// 2 次元インデックス `(i, j)` を 1 次元インデックスに変換する。
/// 配列は `nx` を幅とする行優先で格納されていると仮定する。
#[inline]
pub fn idx(i: usize, j: usize, nx: usize) -> usize {
    j * nx + i
}

/// 周期境界条件を考慮したインデックス変換。
/// 負の値 `i` も許容し，`0..n-1` の範囲へラップする。
#[inline]
pub fn pid(i: isize, n: usize) -> usize {
    (i.rem_euclid(n as isize)) as usize
}
