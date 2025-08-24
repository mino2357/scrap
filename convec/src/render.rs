use crate::config::{Colormap, Interp, OutFmt, OutputCfg, ScaleCfg};
use crate::utils::{idx, pid};
use anyhow::{Context, Result};
use image::{ImageBuffer, Rgb};
use std::fs;
use std::io::{self, Write};
use std::path::PathBuf;

pub struct FrameWriter {
    pub cfg: OutputCfg,
    nx: usize,
    ny: usize,
    next_id: usize,
}

impl FrameWriter {
    pub fn new(mut cfg: OutputCfg, nx: usize, ny: usize) -> Result<Self> {
        if cfg.enable {
            fs::create_dir_all(&cfg.dir).with_context(|| format!("cannot create {}", cfg.dir))?;
        }
        if cfg.out_w.is_none() {
            cfg.out_w = Some(nx);
        }
        if cfg.out_h.is_none() {
            cfg.out_h = Some(ny);
        }
        Ok(Self {
            cfg,
            nx,
            ny,
            next_id: 0,
        })
    }
    pub fn maybe_write(&mut self, q: &Vec<f64>, step: usize, time: f64) -> Result<()> {
        if !self.cfg.enable {
            return Ok(());
        }
        if step % self.cfg.stride == 0 {
            self.write_frame(q, time)?;
        }
        Ok(())
    }
    pub fn write_frame(&mut self, q: &Vec<f64>, time: f64) -> Result<()> {
        let id = self.cfg.start_index + self.next_id;
        let fname = format!("{}_{:06}", self.cfg.prefix, id);
        let path = PathBuf::from(&self.cfg.dir).join(match self.cfg.format {
            OutFmt::Png => format!("{fname}.png"),
            OutFmt::Ppm => format!("{fname}.ppm"),
        });

        let (vmin, vmax) = match self.cfg.scale {
            ScaleCfg::Fixed { min, max } => (min, max),
            ScaleCfg::Auto => {
                let mut mn = q[0];
                let mut mx = q[0];
                for &v in q.iter() {
                    if v < mn {
                        mn = v
                    };
                    if v > mx {
                        mx = v
                    };
                }
                if (mx - mn).abs() < 1e-14 {
                    (mn - 0.5, mx + 0.5)
                } else {
                    (mn, mx)
                }
            }
        };

        let w = self.cfg.out_w.unwrap();
        let h = self.cfg.out_h.unwrap();
        let mut buf: Vec<u8> = vec![0; w * h * 3];

        for py in 0..h {
            let gy = (py as f64 + 0.5) * (self.ny as f64) / (h as f64);
            let (j0, tj) = (gy.floor() as isize, gy.fract());
            for px in 0..w {
                let gx = (px as f64 + 0.5) * (self.nx as f64) / (w as f64);
                let (i0, ti) = (gx.floor() as isize, gx.fract());
                let val = match self.cfg.interp {
                    Interp::Nearest => {
                        let ii = pid(gx.round() as isize, self.nx);
                        let jj = pid(gy.round() as isize, self.ny);
                        q[idx(ii, jj, self.nx)]
                    }
                    Interp::Bilinear => {
                        let i0u = pid(i0, self.nx);
                        let i1u = pid(i0 + 1, self.nx);
                        let j0u = pid(j0, self.ny);
                        let j1u = pid(j0 + 1, self.ny);
                        let q00 = q[idx(i0u, j0u, self.nx)];
                        let q10 = q[idx(i1u, j0u, self.nx)];
                        let q01 = q[idx(i0u, j1u, self.nx)];
                        let q11 = q[idx(i1u, j1u, self.nx)];
                        (1.0 - tj) * ((1.0 - ti) * q00 + ti * q10)
                            + tj * ((1.0 - ti) * q01 + ti * q11)
                    }
                };
                let norm = ((val - vmin) / (vmax - vmin)).clamp(0.0, 1.0);
                let (r, g, b) = match self.cfg.colormap {
                    Colormap::Gray => {
                        let c = (255.0 * norm) as u8;
                        (c, c, c)
                    }
                    Colormap::Turbo => turbo_rgb(norm),
                };
                let yy = if self.cfg.flip_y { h - 1 - py } else { py };
                let p = (yy * w + px) * 3;
                buf[p] = r;
                buf[p + 1] = g;
                buf[p + 2] = b;
            }
        }
        if self.cfg.grid && self.cfg.grid_step > 0 {
            draw_grid(
                &mut buf,
                w,
                h,
                self.nx,
                self.ny,
                self.cfg.grid_step,
                self.cfg.grid_thick,
                [200, 200, 200],
            );
        }
        if self.cfg.axes {
            draw_rect(&mut buf, w, h, 0, 0, w - 1, h - 1, 3, [0, 0, 0]);
        }
        // overlay current time in seconds at top-left
        let label = format!("t = {:.2}", time);
        draw_time_label(&mut buf, w, h, &label);
        if self.cfg.colorbar {
            draw_colorbar(&mut buf, w, h, self.cfg.colormap);
        }

        match self.cfg.format {
            OutFmt::Png => {
                let img: ImageBuffer<Rgb<u8>, _> =
                    ImageBuffer::from_raw(w as u32, h as u32, buf).unwrap();
                img.save(&path)
                    .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
            }
            OutFmt::Ppm => {
                let mut file = fs::File::create(&path)?;
                write!(file, "P6\n{} {}\n255\n", w, h)?;
                file.write_all(&buf)?;
            }
        }
        self.next_id += 1;
        Ok(())
    }
}

fn draw_rect(
    buf: &mut [u8],
    w: usize,
    h: usize,
    x0: usize,
    y0: usize,
    x1: usize,
    y1: usize,
    thick: usize,
    color: [u8; 3],
) {
    for t in 0..thick {
        for x in x0..=x1 {
            set_px(buf, w, h, x, y0.saturating_add(t).min(h - 1), color);
            set_px(buf, w, h, x, y1.saturating_sub(t), color);
        }
        for y in y0..=y1 {
            set_px(buf, w, h, x0.saturating_add(t).min(w - 1), y, color);
            set_px(buf, w, h, x1.saturating_sub(t), y, color);
        }
    }
}
fn draw_grid(
    buf: &mut [u8],
    w: usize,
    h: usize,
    nx: usize,
    ny: usize,
    step: usize,
    thick: usize,
    color: [u8; 3],
) {
    // 垂直線
    for k in (0..=nx).step_by(step) {
        let mut xpix = ((k as f64) * (w as f64) / (nx as f64)).round() as usize;
        if xpix >= w {
            xpix = w - 1;
        }
        for y in 0..h {
            for t in 0..thick {
                let x = xpix.saturating_add(t).min(w - 1);
                set_px(buf, w, h, x, y, color);
            }
        }
    }
    // 水平線
    for k in (0..=ny).step_by(step) {
        let mut ypix = ((k as f64) * (h as f64) / (ny as f64)).round() as usize;
        if ypix >= h {
            ypix = h - 1;
        }
        for x in 0..w {
            for t in 0..thick {
                let y = ypix.saturating_add(t).min(h - 1);
                set_px(buf, w, h, x, y, color);
            }
        }
    }
}
#[inline]
fn set_px(buf: &mut [u8], w: usize, h: usize, x: usize, y: usize, c: [u8; 3]) {
    if x < w && y < h {
        let p = (y * w + x) * 3;
        buf[p] = c[0];
        buf[p + 1] = c[1];
        buf[p + 2] = c[2];
    }
}
fn draw_colorbar(buf: &mut [u8], w: usize, h: usize, cmap: Colormap) {
    let cbw = (w as f64 * 0.05).round() as usize;
    let x0 = w.saturating_sub(cbw);
    for y in 0..h {
        let t = 1.0 - (y as f64 + 0.5) / (h as f64);
        let (r, g, b) = match cmap {
            Colormap::Gray => {
                let c = (255.0 * t) as u8;
                (c, c, c)
            }
            Colormap::Turbo => turbo_rgb(t),
        };
        for x in x0..w {
            let p = (y * w + x) * 3;
            buf[p] = r;
            buf[p + 1] = g;
            buf[p + 2] = b;
        }
    }
}
fn turbo_rgb(t: f64) -> (u8, u8, u8) {
    let x = t.clamp(0.0, 1.0);
    let r =
        34.61 + x * (1172.33 + x * (-10793.56 + x * (33300.12 + x * (-38394.49 + x * 14825.05))));
    let g = 23.31 + x * (557.33 + x * (1225.33 + x * (-3574.96 + x * (4520.31 + x * (-1974.13)))));
    let b = 27.2 + x * (321.15 + x * (1537.82 + x * (-4579.07 + x * (5496.05 + x * (-2163.56)))));
    let to8 = |v: f64| -> u8 { v.round().clamp(0.0, 255.0) as u8 };
    (to8(r), to8(g), to8(b))
}

fn draw_time_label(buf: &mut [u8], w: usize, h: usize, text: &str) {
    let char_w = 5;
    let char_h = 7;
    let spacing = 1;
    let width = text.chars().count() * (char_w + spacing) - spacing;
    let x0 = 5usize;
    let y0 = 5usize;
    fill_rect(
        buf,
        w,
        h,
        x0.saturating_sub(2),
        y0.saturating_sub(2),
        x0 + width + 1,
        y0 + char_h + 1,
        [255, 255, 255],
    );
    draw_text(buf, w, h, x0, y0, text, [0, 0, 0]);
}

fn fill_rect(
    buf: &mut [u8],
    w: usize,
    h: usize,
    x0: usize,
    y0: usize,
    x1: usize,
    y1: usize,
    color: [u8; 3],
) {
    for y in y0..=y1 {
        for x in x0..=x1 {
            set_px(buf, w, h, x, y, color);
        }
    }
}

fn draw_text(
    buf: &mut [u8],
    w: usize,
    h: usize,
    mut x: usize,
    y: usize,
    text: &str,
    color: [u8; 3],
) {
    for ch in text.chars() {
        x += draw_char(buf, w, h, x, y, ch, color) + 1;
    }
}

fn draw_char(
    buf: &mut [u8],
    w: usize,
    h: usize,
    x: usize,
    y: usize,
    ch: char,
    color: [u8; 3],
) -> usize {
    let glyph: [u8; 7] = match ch {
        '0' => [
            0b01110, 0b10001, 0b10011, 0b10101, 0b11001, 0b10001, 0b01110,
        ],
        '1' => [
            0b00100, 0b01100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110,
        ],
        '2' => [
            0b01110, 0b10001, 0b00001, 0b00010, 0b00100, 0b01000, 0b11111,
        ],
        '3' => [
            0b11110, 0b00001, 0b00001, 0b01110, 0b00001, 0b00001, 0b11110,
        ],
        '4' => [
            0b10010, 0b10010, 0b10010, 0b11111, 0b00010, 0b00010, 0b00010,
        ],
        '5' => [
            0b11111, 0b10000, 0b10000, 0b11110, 0b00001, 0b00001, 0b11110,
        ],
        '6' => [
            0b01110, 0b10000, 0b10000, 0b11110, 0b10001, 0b10001, 0b01110,
        ],
        '7' => [
            0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b01000, 0b01000,
        ],
        '8' => [
            0b01110, 0b10001, 0b10001, 0b01110, 0b10001, 0b10001, 0b01110,
        ],
        '9' => [
            0b01110, 0b10001, 0b10001, 0b01111, 0b00001, 0b00001, 0b01110,
        ],
        '.' => [
            0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b01100, 0b01100,
        ],
        '=' => [
            0b00000, 0b00000, 0b11111, 0b00000, 0b11111, 0b00000, 0b00000,
        ],
        't' | 'T' => [
            0b01110, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100,
        ],
        _ => [0; 7],
    };
    for (j, row) in glyph.iter().enumerate() {
        for i in 0..5 {
            if (row >> (4 - i)) & 1 == 1 {
                set_px(buf, w, h, x + i, y + j, color);
            }
        }
    }
    5
}
