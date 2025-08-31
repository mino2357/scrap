//! Shared core for WENO5-family reconstruction (coefficients and beta).
//! 5 点サブステンシルの左右再構成で用いる候補多項式値 `cval` と
//! Jiang–Shu 型の滑らかさ指標 `beta` を返すユーティリティ。

pub const EPS5: f64 = 1e-6;
pub const D5: [f64; 3] = [1.0/10.0, 3.0/5.0, 3.0/10.0];

#[inline]
pub fn cvals_betas_5(arr: &[f64; 5]) -> ([f64; 3], [f64; 3]) {
    let im2 = arr[0];
    let im1 = arr[1];
    let i0 = arr[2];
    let ip1 = arr[3];
    let ip2 = arr[4];
    let c0 = (2.0 * im2 - 7.0 * im1 + 11.0 * i0) / 6.0;
    let c1 = (-im1 + 5.0 * i0 + 2.0 * ip1) / 6.0;
    let c2 = (2.0 * i0 + 5.0 * ip1 - ip2) / 6.0;
    let b0 = (13.0 / 12.0) * (im2 - 2.0 * im1 + i0).powi(2)
        + 0.25 * (im2 - 4.0 * im1 + 3.0 * i0).powi(2);
    let b1 = (13.0 / 12.0) * (im1 - 2.0 * i0 + ip1).powi(2) + 0.25 * (im1 - ip1).powi(2);
    let b2 = (13.0 / 12.0) * (i0 - 2.0 * ip1 + ip2).powi(2)
        + 0.25 * (3.0 * i0 - 4.0 * ip1 + ip2).powi(2);
    ([c0, c1, c2], [b0, b1, b2])
}
