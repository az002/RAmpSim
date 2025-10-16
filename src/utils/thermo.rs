use bam::Record as BamRecord;

const R_GAS_CAL: f64 = 1.987_204_258_64083;

#[inline]
fn comp_code(x: u8) -> u8 { x ^ 0b11 }

#[derive(Clone, Copy, Debug)]
pub struct HS {
    pub dh: f64, //kcal/mol
    pub ds: f64, //cal/(K mol)
}

static NT2BITS: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0; t[b'C' as usize] = 1; t[b'G' as usize] = 2; t[b'T' as usize] = 3;
    t[b'a' as usize] = 0; t[b'c' as usize] = 1; t[b'g' as usize] = 2; t[b't' as usize] = 3;
    t
};

#[inline]
fn is_wc_codes(top: u8, bot: u8) -> bool {
    top < 4 && bot < 4 && bot == comp_code(top)
}

fn hs_matrix() -> [[Option<HS>; 4]; 4] {
    let mut m: [[Option<HS>; 4]; 4] = [[None; 4]; 4];

    let set = |m: &mut [[Option<HS>;4];4], x: u8, y: u8, dh: f64, ds: f64| {
        m[x as usize][y as usize] = Some(HS { dh, ds });
    };

    set(&mut m, 0, 0, -7.9,  -22.2);
    set(&mut m, 0, 3, -7.2,  -20.4);
    set(&mut m, 3, 0, -7.2,  -21.3);
    set(&mut m, 1, 0, -8.5,  -22.7);
    set(&mut m, 2, 3, -8.4,  -22.4);
    set(&mut m, 1, 3, -7.8,  -21.0);
    set(&mut m, 2, 0, -8.2,  -22.2);
    set(&mut m, 1, 2, -10.6, -27.2);
    set(&mut m, 2, 1, -9.8,  -24.4);
    set(&mut m, 2, 2, -8.0,  -19.9);

    for x in 0..4u8 {
        for y in 0..4u8 {
            if m[x as usize][y as usize].is_none() {
                let xr = comp_code(y);
                let yr = comp_code(x);
                m[x as usize][y as usize] = m[xr as usize][yr as usize];
            }
        }
    }

    m
}

pub enum MismatchStrategy {
    SkipStacking,
}

pub fn score_bam_alignment_multiT(
    record: &BamRecord,
    temp_c: f64,
    strategy: MismatchStrategy,
) -> f64 {
    let t_k = temp_c + 273.15;
    let hs_tab = hs_matrix();

    let mut last_top: Option<u8> = None;
    let mut last_bot: Option<u8> = None; 

    let mut sum_dh = 0.0f64;   // kcal/mol
    let mut sum_ds = 0.0f64;   // cal/(K mol)

    for e in record.alignment_entries().expect("iter entries") {
        if !e.is_aln_match() {
            // ignore indels
            continue;
        }
        let top_code = NT2BITS[e.record_nt().unwrap_or(b'N') as usize];
        let bot_code = NT2BITS[e.ref_nt().unwrap_or(b'N') as usize];

        if let (Some(lp), Some(lb)) = (last_top, last_bot) {
            let prev_wc = is_wc_codes(lp, lb);
            let curr_wc = is_wc_codes(top_code, bot_code);

            match (prev_wc, curr_wc) {
                (true, true) => {
                    if let Some(hs) = hs_tab[lp as usize][top_code as usize] {
                        sum_dh += hs.dh;
                        sum_ds += hs.ds;
                    }
                }
                _ => {
                    // Step touches a mismatch
                    match strategy {
                        MismatchStrategy::SkipStacking => {
                        }
                    }
                }
            }
        }

        last_top = Some(top_code);
        last_bot = Some(bot_code);
    }

    let dg_stack = sum_dh - (t_k * (sum_ds / 1000.0)); // cal to kcal
    let r_kcal = R_GAS_CAL / 1000.0;
    (-dg_stack / (r_kcal * t_k)).exp()
}
