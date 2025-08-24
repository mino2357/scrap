pub trait Scheme {
    fn rhs(
        &self,
        q: &[f64],
        u: &[f64],
        v: &[f64],
        dx: f64,
        dy: f64,
        nx: usize,
        ny: usize,
        out: &mut [f64],
    );
}

pub mod centered8;
pub mod cip;
pub mod mp5;
pub mod tvd;
pub mod upwind1;
pub mod weno5;
pub use centered8::Centered8;
pub use cip::{Cip, CipB, CipCsl, CipCsl2, CipCsl2Mh};
pub use mp5::Mp5;
pub use tvd::{TvdMinmod, TvdVanLeer};
pub use upwind1::Upwind1;
pub use weno5::Weno5Js;
