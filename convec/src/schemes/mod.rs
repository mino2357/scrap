pub trait Scheme {
    fn rhs(&self, q:&[f64], u:&[f64], v:&[f64], dx:f64, dy:f64, nx:usize, ny:usize, out:&mut [f64]);
}

pub mod centered8;
pub mod weno5;
pub use centered8::Centered8;
pub use weno5::Weno5Js;
