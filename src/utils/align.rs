use std::iter::repeat;
use serde::{Serialize, Deserialize};
pub type TextSlice<'a> = &'a [u8];

#[allow(non_snake_case)]
pub struct aligner{
    xo: i32,
    xe: i32,
    yo: i32,
    ye: i32,
    ma: i32,
    mm: i32,
}

pub const MIN_SCORE: i32 = -858_993_459;

impl aligner {
    pub fn new(xo: i32, xe: i32, yo: i32, ye: i32, ma: i32, mm: i32) -> Self {
        aligner {
            xo,
            xe,
            yo,
            ye,
            ma,
            mm,
        }
    }

    pub fn semiglobal(&self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
        let (m, n) = (x.len(), y.len());
        let mut Y = vec![vec![0; m + 1]; 2];
        let mut S = vec![vec![0; m + 1]; 2];
        let mut traceback = Traceback::with_capacity(m, n);
        traceback.init(m,n);

        Y[0][0] = MIN_SCORE;
        S[0][0] = 0;
        let mut tb = TracebackCell::new();
        tb.set_all(TB_START);
        traceback.set(0, 0, tb);     

        for i in 1..=m {
            let mut tb = TracebackCell::new();
            Y[0][i] = self.yo + self.ye * ((i - 1) as i32);
            S[0][i] = Y[0][i];
            tb.set_all(TB_INS);
            traceback.set(i, 0, tb);
        }

        let mut best_pos = 0;
        let mut best_score = S[0][m];

        for j in 1..=n {
            let curr = j % 2;
            let prev = 1 - curr;
            {
                S[curr][0] = 0;
                Y[curr][0] = MIN_SCORE;
                let mut tb = TracebackCell::new();
                tb.set_all(TB_START);
                traceback.set(0, j, tb);
            }
            for i in 1..=m {
                S[curr][i] = MIN_SCORE;
                Y[curr][i] = MIN_SCORE;
                let p = x[i - 1];
                let q = y[j - 1];
                let m_score = S[prev][i - 1] + match p == q {
                    true => self.ma,
                    false => self.mm,
                };
                let i_score = Y[curr][i - 1] + self.ye;
                let s_score = S[curr][i - 1] + self.yo;

                let mut best_s_score = S[curr][i];
                if m_score > best_s_score {
                    best_s_score = m_score;
                    tb.set_s_bits(if p == q { TB_MATCH } else { TB_SUBST });
                }
                if i_score > best_s_score {
                    best_s_score = i_score;
                    tb.set_s_bits(TB_INS);
                }
                if s_score > best_s_score {
                    best_s_score = s_score;
                    tb.set_s_bits(traceback.get(i - 1, j).get_s_bits());
                }

                let best_i_score;
                if i_score > s_score {
                    best_i_score = i_score;
                    tb.set_i_bits(TB_INS);
                } else {
                    best_i_score = s_score;
                    tb.set_i_bits(traceback.get(i - 1, j).get_s_bits());
                }
                
                let mut tb = TracebackCell::new();
                if best_s_score == m_score {
                    tb.set_s_bits(if p == q { TB_MATCH } else { TB_SUBST });
                } else if best_s_score == best_i_score {
                    tb.set_s_bits(TB_INS);
                } else {
                    tb.set_s_bits(TB_DEL);
                }
                S[curr][i] = best_s_score;
                Y[curr][i] = best_i_score;
                traceback.set(i, j, tb);
                if i == m && S[curr][i] > best_score {
                    best_score = S[curr][i];
                    best_pos = j;
                }
            }
        }

        let mut i = m;
        let mut j = best_pos;
        let mut last_layer = traceback.get(i, j).get_s_bits();
        loop {
            let next_layer: u16;
            match last_layer {
                TB_START => break,
                TB_INS => {
                    next_layer = traceback.get(i, j).get_i_bits();
                    i -= 1;
                }
                TB_MATCH => {
                    next_layer = traceback.get(i - 1, j - 1).get_s_bits();
                    i -= 1;
                    j -= 1;
                }
                TB_SUBST => {
                    next_layer = traceback.get(i - 1, j - 1).get_s_bits();
                    i -= 1;
                    j -= 1;
                }
                _ => panic!("Didn't expect this!"),
            }
            last_layer = next_layer;
        }

        Alignment{
            score: best_score,
            ystart: j,
            xstart: 0,
            yend: best_pos,
            xend: m,
            ylen: n,
            xlen: m,
        }
    }  
}

#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub struct TracebackCell {
    v: u16,
}

// Traceback bit positions (LSB)
const I_POS: u8 = 0; // Meaning bits 0,1,2,3 corresponds to I and so on
const D_POS: u8 = 4;
const S_POS: u8 = 8;

// Traceback moves
const TB_START: u16 = 0b0000;
const TB_INS: u16 = 0b0001;
const TB_DEL: u16 = 0b0010;
const TB_SUBST: u16 = 0b0011;
const TB_MATCH: u16 = 0b0100;

const TB_XCLIP_PREFIX: u16 = 0b0101; // prefix clip of x
const TB_XCLIP_SUFFIX: u16 = 0b0110; // suffix clip of x
const TB_YCLIP_PREFIX: u16 = 0b0111; // prefix clip of y
const TB_YCLIP_SUFFIX: u16 = 0b1000; // suffix clip of y

const TB_MAX: u16 = 0b1000; // Useful in checking that the
                            // TB value we got is a valid one
impl TracebackCell {
    /// Initialize a blank traceback cell
    #[inline(always)]
    pub fn new() -> TracebackCell {
        Default::default()
    }

    /// Sets 4 bits [pos, pos+4) with the 4 LSBs of value
    #[inline(always)]
    fn set_bits(&mut self, pos: u8, value: u16) {
        let bits: u16 = (0b1111) << pos;
        assert!(
            value <= TB_MAX,
            "Expected a value <= TB_MAX while setting traceback bits"
        );
        self.v = (self.v & !bits) // First clear the bits
            | (value << pos) // And set the bits
    }

    #[inline(always)]
    pub fn set_i_bits(&mut self, value: u16) {
        // Traceback corresponding to matrix I
        self.set_bits(I_POS, value);
    }

    #[inline(always)]
    pub fn set_d_bits(&mut self, value: u16) {
        // Traceback corresponding to matrix D
        self.set_bits(D_POS, value);
    }

    #[inline(always)]
    pub fn set_s_bits(&mut self, value: u16) {
        // Traceback corresponding to matrix S
        self.set_bits(S_POS, value);
    }

    // Gets 4 bits [pos, pos+4) of v
    #[inline(always)]
    fn get_bits(self, pos: u8) -> u16 {
        (self.v >> pos) & (0b1111)
    }

    #[inline(always)]
    pub fn get_i_bits(self) -> u16 {
        self.get_bits(I_POS)
    }

    #[inline(always)]
    pub fn get_d_bits(self) -> u16 {
        self.get_bits(D_POS)
    }

    #[inline(always)]
    pub fn get_s_bits(self) -> u16 {
        self.get_bits(S_POS)
    }

    /// Set all matrices to the same value.
    pub fn set_all(&mut self, value: u16) {
        self.set_i_bits(value);
        self.set_d_bits(value);
        self.set_s_bits(value);
    }
}

/// Internal traceback.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
struct Traceback {
    rows: usize,
    cols: usize,
    matrix: Vec<TracebackCell>,
}

impl Traceback {
    fn with_capacity(m: usize, n: usize) -> Self {
        let rows = m + 1;
        let cols = n + 1;
        Traceback {
            rows,
            cols,
            matrix: Vec::with_capacity(rows * cols),
        }
    }

    fn init(&mut self, m: usize, n: usize) {
        self.matrix.clear();
        let mut start = TracebackCell::new();
        start.set_all(TB_START);
        // set every cell to start
        self.resize(m, n, start);
    }

    #[inline(always)]
    fn set(&mut self, i: usize, j: usize, v: TracebackCell) {
        debug_assert!(i < self.rows);
        debug_assert!(j < self.cols);
        self.matrix[i * self.cols + j] = v;
    }

    #[inline(always)]
    fn get(&self, i: usize, j: usize) -> &TracebackCell {
        debug_assert!(i < self.rows);
        debug_assert!(j < self.cols);
        &self.matrix[i * self.cols + j]
    }

    fn get_mut(&mut self, i: usize, j: usize) -> &mut TracebackCell {
        debug_assert!(i < self.rows);
        debug_assert!(j < self.cols);
        &mut self.matrix[i * self.cols + j]
    }

    fn resize(&mut self, m: usize, n: usize, v: TracebackCell) {
        self.rows = m + 1;
        self.cols = n + 1;
        self.matrix.resize(self.rows * self.cols, v);
    }
}

pub struct Alignment {
    pub score: i32,
    pub ystart: usize,
    pub xstart: usize,
    pub yend: usize,
    pub xend: usize,
    pub ylen: usize,
    pub xlen: usize,
}

// const DEFAULT_ALIGNER_CAPACITY: usize = 200;

// pub struct Scoring {
//     pub x_open: i32,
//     pub x_extend: i32,
//     pub y_open: i32,
//     pub y_extend: i32,
//     pub ma: i32,
//     pub mm: i32,
//     pub xclip_prefix: i32,
//     pub xclip_suffix: i32,
//     pub yclip_prefix: i32,
//     pub yclip_suffix: i32,
// }

// impl Scoring {
//     pub fn new(x_open: i32, x_extend: i32, y_open: i32, y_extend: i32, ma: i32, mm: i32) -> Self {
//         assert!(x_open <= 0, "gap_open can't be positive");
//         assert!(x_extend <= 0, "gap_extend can't be positive");
//         assert!(y_open <= 0, "gap_open can't be positive");
//         assert!(y_extend <= 0, "gap_extend can't be positive");

//         Scoring {
//             x_open,
//             x_extend,
//             y_open,
//             y_extend,
//             ma,
//             mm,
//             xclip_prefix: MIN_SCORE,
//             xclip_suffix: MIN_SCORE,
//             yclip_prefix: MIN_SCORE,
//             yclip_suffix: MIN_SCORE,
//         }
//     }

//     pub fn score(&self, a: u8, b: u8) -> i32 {
//         if a == b {
//             self.ma
//         } else {
//             self.mm
//         }
//     }
// }

// pub struct Aligner {
//     I: [Vec<i32>; 2],
//     D: [Vec<i32>; 2],
//     S: [Vec<i32>; 2],
//     Lx: Vec<usize>,
//     Ly: Vec<usize>,
//     Sn: Vec<i32>,
//     traceback: Traceback,
//     scoring: Scoring,
// }

// impl Aligner {
//     /// Create new aligner instance with given gap open and gap extend penalties
//     /// and the score function.
//     ///
//     /// # Arguments
//     ///
//     /// * `gap_open` - the score for opening a gap (should be negative)
//     /// * `gap_extend` - the score for extending a gap (should be negative)
//     /// * `match_fn` - function that returns the score for substitutions
//     ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
//     pub fn new(x_open: i32, x_extend: i32, y_open: i32, y_extend: i32, ma: i32, mm: i32) -> Self {
//         assert!(x_open <= 0, "gap_open can't be positive");
//         assert!(x_extend <= 0, "gap_extend can't be positive");
//         assert!(y_open <= 0, "gap_open can't be positive");
//         assert!(y_extend <= 0, "gap_extend can't be positive");

//         Aligner {
//             I: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
//             D: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
//             S: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
//             Lx: Vec::with_capacity(n + 1),
//             Ly: Vec::with_capacity(m + 1),
//             Sn: Vec::with_capacity(m + 1),
//             traceback: Traceback::with_capacity(m, n),
//             scoring: Scoring::new(x_open, x_extend, y_open, y_extend, ma, mm),
//         }
//     }


//     pub fn semiglobal(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment{
//         let (m, n) = (x.len(), y.len());
//         self.traceback.init(m, n);
//         self.scoring.yclip_prefix = 0;
//         self.scoring.yclip_suffix = 0;

//         // Set the initial conditions
//         // We are repeating some work, but that's okay!
//         for k in 0..2 {
//             self.I[k].clear();
//             self.D[k].clear();
//             self.S[k].clear();

//             self.D[k].extend(repeat(MIN_SCORE).take(m + 1));
//             self.I[k].extend(repeat(MIN_SCORE).take(m + 1));
//             self.S[k].extend(repeat(MIN_SCORE).take(m + 1));

//             self.S[k][0] = 0;

//             if k == 0 {
//                 let mut tb = TracebackCell::new();
//                 tb.set_all(TB_START);
//                 self.traceback.set(0, 0, tb);
//                 self.Lx.clear();
//                 self.Lx.extend(repeat(0usize).take(n + 1));
//                 self.Ly.clear();
//                 self.Ly.extend(repeat(0usize).take(m + 1));
//                 self.Sn.clear();
//                 self.Sn.extend(repeat(MIN_SCORE).take(m + 1));
//                 self.Sn[0] = 0;
//                 self.Ly[0] = n;
//             }

//             for i in 1..=m {
//                 let mut tb = TracebackCell::new();
//                 tb.set_all(TB_START);
//                 if i == 1 {
//                     self.I[k][i] = self.scoring.y_open;
//                     tb.set_i_bits(TB_START);
//                 } else {
//                     // Insert all i characters
//                     self.I[k][i] = self.scoring.y_open + self.scoring.y_extend * ((i-1) as i32);
//                     tb.set_i_bits(TB_INS);
//                 }

//                 self.S[k][i] = self.I[k][i];
//                 tb.set_s_bits(TB_INS);

//                 if k == 0 {
//                     self.traceback.set(i, 0, tb);
//                 }
//                 // Track the score if we do suffix clip (y) from here
//                 if self.S[k][i] + self.scoring.yclip_suffix > self.Sn[i] {
//                     self.Sn[i] = self.S[k][i] + self.scoring.yclip_suffix;
//                     self.Ly[i] = n;
//                 }
//             }
//         }

//         for j in 1..=n {
//             let curr = j % 2;
//             let prev = 1 - curr;

//             {
//                 // Handle i = 0 case
//                 let mut tb = TracebackCell::new();
//                 self.I[curr][0] = MIN_SCORE;

//                 if j == 1 {
//                     self.D[curr][0] = self.scoring.x_open;
//                     tb.set_d_bits(TB_START);
//                 } else {
//                     // Delete all j characters
//                     let d_score = self.scoring.x_open + self.scoring.x_extend * ((j-1) as i32);
//                     let c_score =
//                         self.scoring.yclip_prefix + self.scoring.x_open;
//                     if d_score > c_score {
//                         self.D[curr][0] = d_score;
//                         tb.set_d_bits(TB_DEL);
//                     } else {
//                         self.D[curr][0] = c_score;
//                         tb.set_d_bits(TB_YCLIP_PREFIX);
//                     }
//                 }
//                 if self.D[curr][0] > self.scoring.yclip_prefix {
//                     self.S[curr][0] = self.D[curr][0];
//                     tb.set_s_bits(TB_DEL);
//                 } else {
//                     self.S[curr][0] = self.scoring.yclip_prefix;
//                     tb.set_s_bits(TB_YCLIP_PREFIX);
//                 }

//                 if j == n && self.Sn[0] > self.S[curr][0] {
//                     // Check if the suffix clip score is better
//                     self.S[curr][0] = self.Sn[0];
//                     tb.set_s_bits(TB_YCLIP_SUFFIX);
//                 // Track the score if we do suffix clip (y) from here
//                 } else if self.S[curr][0] + self.scoring.yclip_suffix > self.Sn[0] {
//                     self.Sn[0] = self.S[curr][0] + self.scoring.yclip_suffix;
//                     self.Ly[0] = n - j;
//                 }

//                 self.traceback.set(0, j, tb);
//             }

//             for i in 1..=m {
//                 self.S[curr][i] = MIN_SCORE;
//             }

//             let q = y[j - 1];

//             for i in 1..m + 1 {
//                 let p = x[i - 1];
//                 let mut tb = TracebackCell::new();

//                 let m_score = self.S[prev][i - 1] + self.scoring.score(p, q);

//                 let i_score = self.I[curr][i - 1] + self.scoring.y_extend;
//                 let s_score = self.S[curr][i - 1] + self.scoring.y_open;
//                 let best_i_score;
//                 if i_score > s_score {
//                     best_i_score = i_score;
//                     tb.set_i_bits(TB_INS);
//                 } else {
//                     best_i_score = s_score;
//                     tb.set_i_bits(self.traceback.get(i - 1, j).get_s_bits());
//                 }

//                 let d_score = self.D[prev][i] + self.scoring.x_extend;
//                 let s_score = self.S[prev][i] + self.scoring.x_open;
//                 let best_d_score;
//                 if d_score > s_score {
//                     best_d_score = d_score;
//                     tb.set_d_bits(TB_DEL);
//                 } else {
//                     best_d_score = s_score;
//                     tb.set_d_bits(self.traceback.get(i, j - 1).get_s_bits());
//                 }

//                 tb.set_s_bits(TB_XCLIP_SUFFIX);
//                 let mut best_s_score = self.S[curr][i];

//                 if m_score > best_s_score {
//                     best_s_score = m_score;
//                     tb.set_s_bits(if p == q { TB_MATCH } else { TB_SUBST });
//                 }

//                 if best_i_score > best_s_score {
//                     best_s_score = best_i_score;
//                     tb.set_s_bits(TB_INS);
//                 }

//                 if best_d_score > best_s_score {
//                     best_s_score = best_d_score;
//                     tb.set_s_bits(TB_DEL);
//                 }

//                 let yclip_score = self.scoring.yclip_prefix
//                     + self.scoring.y_open
//                     + self.scoring.y_extend * ((i-1) as i32);
//                 if yclip_score > best_s_score {
//                     best_s_score = yclip_score;
//                     tb.set_s_bits(TB_YCLIP_PREFIX);
//                 }

//                 self.S[curr][i] = best_s_score;
//                 self.I[curr][i] = best_i_score;
//                 self.D[curr][i] = best_d_score;

//                 // Track the score if we do suffix clip (y) from here
//                 if self.S[curr][i] + self.scoring.yclip_suffix > self.Sn[i] {
//                     self.Sn[i] = self.S[curr][i] + self.scoring.yclip_suffix;
//                     self.Ly[i] = n - j;
//                 }

//                 self.traceback.set(i, j, tb);
//             }
//         }

//         // Handle suffix clipping in the j=n case
//         for i in 0..=m {
//             let j = n;
//             let curr = j % 2;
//             if self.Sn[i] > self.S[curr][i] {
//                 self.S[curr][i] = self.Sn[i];
//                 self.traceback.get_mut(i, j).set_s_bits(TB_YCLIP_SUFFIX);
//             }
//         }

//         // Since there could be a change in the last column of S,
//         // recompute the last column of I as this could also change
//         for i in 1..=m {
//             let j = n;
//             let curr = j % 2;
//             let s_score = self.S[curr][i - 1] + self.scoring.y_open;
//             if s_score > self.I[curr][i] {
//                 self.I[curr][i] = s_score;
//                 let s_bit = self.traceback.get(i - 1, j).get_s_bits();
//                 self.traceback.get_mut(i, j).set_i_bits(s_bit);
//             }
//             if s_score > self.S[curr][i] {
//                 self.S[curr][i] = s_score;
//                 self.traceback.get_mut(i, j).set_s_bits(TB_INS);
//             }
//         }

//         let mut i = m;
//         let mut j = n;
//         let mut xstart: usize = 0usize;
//         let mut ystart: usize = 0usize;
//         let mut xend = m;
//         let mut yend = n;

//         let mut last_layer = self.traceback.get(i, j).get_s_bits();

//         loop {
//             let next_layer: u16;
//             match last_layer {
//                 TB_START => break,
//                 TB_INS => {
//                     next_layer = self.traceback.get(i, j).get_i_bits();
//                     i -= 1;
//                 }
//                 TB_DEL => {
//                     next_layer = self.traceback.get(i, j).get_d_bits();
//                     j -= 1;
//                 }
//                 TB_MATCH => {
//                     next_layer = self.traceback.get(i - 1, j - 1).get_s_bits();
//                     i -= 1;
//                     j -= 1;
//                 }
//                 TB_SUBST => {
//                     next_layer = self.traceback.get(i - 1, j - 1).get_s_bits();
//                     i -= 1;
//                     j -= 1;
//                 }
//                 TB_XCLIP_PREFIX => {
//                     xstart = i;
//                     i = 0;
//                     next_layer = self.traceback.get(0, j).get_s_bits();
//                 }
//                 TB_XCLIP_SUFFIX => {
//                     i -= self.Lx[j];
//                     xend = i;
//                     next_layer = self.traceback.get(i, j).get_s_bits();
//                 }
//                 TB_YCLIP_PREFIX => {
//                     ystart = j;
//                     j = 0;
//                     next_layer = self.traceback.get(i, 0).get_s_bits();
//                 }
//                 TB_YCLIP_SUFFIX => {
//                     j -= self.Ly[i];
//                     yend = j;
//                     next_layer = self.traceback.get(i, j).get_s_bits();
//                 }
//                 _ => panic!("Didn't expect this!"),
//             }
//             last_layer = next_layer;
//         }

//         Alignment {
//             score: self.S[n % 2][m],
//             ystart,
//             xstart,
//             yend,
//             xend,
//             ylen: n,
//             xlen: m,
//         }
//     }
// }