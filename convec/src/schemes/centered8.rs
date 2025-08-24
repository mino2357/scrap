use crate::schemes::Scheme;
use crate::utils::{idx, pid};

pub struct Centered8;

const C8: [f64; 9] = [ 1.0/280.0, -4.0/105.0, 1.0/5.0, -4.0/5.0, 0.0,
                       4.0/5.0,   -1.0/5.0,   4.0/105.0, -1.0/280.0 ];

fn d1x(f:&[f64], nx:usize, ny:usize, dx:f64, out:&mut [f64]){
    for j in 0..ny {
        for i in 0..nx {
            let mut acc = 0.0;
            for s in -4..=4 { if s!=0 {
                let w = C8[(s+4) as usize];
                let ii = pid(i as isize + s as isize, nx);
                acc += w * f[idx(ii,j,nx)];
            }}
            out[idx(i,j,nx)] = acc / dx;
        }
    }
}
fn d1y(f:&[f64], nx:usize, ny:usize, dy:f64, out:&mut [f64]){
    for j in 0..ny {
        for i in 0..nx {
            let mut acc = 0.0;
            for s in -4..=4 { if s!=0 {
                let w = C8[(s+4) as usize];
                let jj = pid(j as isize + s as isize, ny);
                acc += w * f[idx(i,jj,nx)];
            }}
            out[idx(i,j,nx)] = acc / dy;
        }
    }
}

impl Scheme for Centered8 {
    fn rhs(&self, q:&[f64], u:&[f64], v:&[f64], dx:f64, dy:f64, nx:usize, ny:usize, out:&mut [f64]){
        let mut dqx = vec![0.0; nx*ny];
        let mut dqy = vec![0.0; nx*ny];
        d1x(q, nx, ny, dx, &mut dqx);
        d1y(q, nx, ny, dy, &mut dqy);

        let mut uq = vec![0.0; nx*ny];
        let mut vq = vec![0.0; nx*ny];
        for k in 0..nx*ny { uq[k]=u[k]*q[k]; vq[k]=v[k]*q[k]; }

        let mut dx_uq = vec![0.0; nx*ny];
        let mut dy_vq = vec![0.0; nx*ny];
        d1x(&uq, nx, ny, dx, &mut dx_uq);
        d1y(&vq, nx, ny, dy, &mut dy_vq);

        for k in 0..nx*ny {
            out[k] = -0.5 * ( u[k]*dqx[k] + dx_uq[k] + v[k]*dqy[k] + dy_vq[k] );
        }
    }
}
