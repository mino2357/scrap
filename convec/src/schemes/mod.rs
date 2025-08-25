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

pub mod centered10;
pub mod centered12;
pub mod centered14;
pub mod centered6;
pub mod centered8;
pub mod cip;
pub mod mp5;
pub mod tvd;
pub mod upwind1;
pub mod upwind3x3;
pub mod weno5;
pub mod weno5z;
pub mod weno7z;
pub mod weno9z;
pub use centered6::Centered6;
pub use centered8::Centered8;
pub use centered10::Centered10;
pub use centered12::Centered12;
pub use centered14::Centered14;
pub use cip::{Cip, CipB, CipCsl, CipCsl2, CipCsl2Mh};
pub use mp5::Mp5;
pub use tvd::{TvdMinmod, TvdVanLeer};
pub use upwind1::Upwind1;
pub use upwind3x3::Upwind3x3;
pub use weno5::Weno5Js;
pub use weno5z::Weno5Z;
pub use weno7z::Weno7Z;
pub use weno9z::Weno9Z;
