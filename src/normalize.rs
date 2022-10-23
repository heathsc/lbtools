use std::collections::{HashMap, VecDeque};

use crate::{config::Config, coverage::*, gc::N_GC_BINS};

// const MIN_COUNTS: usize = 1;

fn collect_bin_data(cfg: &Config, rc: RawCounts) -> Vec<Vec<usize>> {
    let mut bin_counts: Vec<Vec<usize>> = vec![Vec::new(); N_GC_BINS as usize];
    for contig in cfg
        .ctg_hash()
        .values()
        .filter(|c| c.use_for_normalization())
    {
        let ctg = contig.name();
        if let Some(raw_cts) = rc.get(ctg) {
            let gc = cfg
                .gc_data()
                .ctg_data(ctg)
                .expect("Missing GC data for contig");
            for (ix, ct) in raw_cts.iter().enumerate() {
                if let Some(j) = gc.gc(ix) {
                    bin_counts[j as usize].push(*ct);
                }
            }
        }
    }
    bin_counts
}

struct Obs {
    n: usize,  // Number of observations
    ix: usize, // Original GC bin
    quartiles: [usize; 3],
}

impl Obs {
    fn new(ix: usize, v: &mut [usize]) -> Option<Self> {
        let n = v.len();
        if n == 0 {
            None
        } else {
            v.sort_unstable();
            let quartiles = [v[n >> 2], v[n >> 1], v[(n * 3) >> 2]];
            Some(Self { n, ix, quartiles })
        }
    }

    fn weight(&self) -> f64 {
        self.n as f64
        //        let sd = ((self.quartiles[2] - self.quartiles[0]).max(1) as f64) / 1.35;
        //        (self.n as f64) / (sd * sd)
    }
}

// Accumulate regression equations
// X is a n x 3 matrix; W and Y are n x 1 vectors
// X'WX is symmetric 3 x 3 matrix and X'WX[3,1] = X'WX[2,2] = X'WX[1,3]
// We calculate just the lower triangle of X'WX
struct Accum {
    xwx: [f64; 6],
    xwy: [f64; 3],
    x0: f64,
    win_size: f64,
}

impl Accum {
    fn new(x0: usize, win_size: usize) -> Self {
        Self {
            xwx: [0.0; 6],
            xwy: [0.0; 3],
            x0: x0 as f64,
            win_size: win_size as f64,
        }
    }

    // Accumulate contribution of observation to LS matrices
    fn accum(&mut self, o: &Obs) {
        let x = (o.ix as f64) - self.x0;
        let d = x.abs() / self.win_size;
        // sanity check
        assert!(d <= 1.0);
        // tricube kernel
        let z = (1.0 - d * d * d);
        // weight = (1/estimate variance) * tricubic kernel
        let w = o.weight() * z * z * z;
        let x2 = x * x;
        let x3 = x * x2;
        let x4 = x2 * x2;
        // Use the median as our y value
        let y = o.quartiles[1] as f64;
        // Accumulate lower triangle of XWX
        self.xwx[0] += w; // Sigma w
        self.xwx[1] += w * x; // w * x
        self.xwx[2] += w * x2; // w * x^2
        self.xwx[3] += w * x2; // w * x^2
        self.xwx[4] += w * x3; // w * x^3
        self.xwx[5] += w * x4; // w * x^4
        self.xwy[0] += w * y; // w * y
        self.xwy[1] += w * x * y; // w * x * y
        self.xwy[2] += w * x2 * y; // w * x^2 * y
    }

    // Calculate Choleshy decomposition of X'WX matrix
    // Decomposition stored in place
    fn chol(&mut self) {
        let x = &mut self.xwx;
        assert!(x[0] > 0.0, "Matrix not PD");
        x[0] = x[0].sqrt(); // L[0,0]
        x[1] /= x[0]; // L[1,0]
        let z = x[2] - x[1] * x[1];
        assert!(z > 0.0, "Matrix not PD");
        x[2] = z.sqrt(); // L[1,1]
        x[3] /= x[0]; // L[2,0]
        x[4] = (x[4] - x[1] * x[3]) / x[2]; // L[2,1]
        let z = x[5] - x[3] * x[3] - x[4] * x[4];
        assert!(z > 0.0, "Matrix not PD");
        x[5] = z.sqrt(); // L[2;2]
    }

    fn solve(&mut self, b: &mut [f64]) {
        self.chol();
        let c = &mut self.xwx;
        let y = &mut self.xwy;
        let a0 = y[0] / c[0];
        let a1 = (y[1] - c[1] * a0) / c[2];
        let a2 = (y[2] - c[3] * a0 - c[4] * a1) / c[5];
        b[2] = a2 / c[5];
        b[1] = (a1 - c[4] * b[2]) / c[2];
        b[0] = (a0 - c[3] * b[2] - c[1] * b[1]) / c[0];
    }
}

struct Fit {
    x: isize,       // centre point of regression
    beta: [f64; 3], // regression coefficients
}

impl Fit {
    // Fit local regression with the observations in obs at the position
    // given by the observation obs[i]
    fn fit_local_regression(obs: &[Obs], i: usize) -> Self {
        // x coordinate of location where we are performing the fit
        let x0 = obs[i].ix;
        // window size (max distance from index location)
        let d = (obs.last().unwrap().ix - x0).max(x0 - obs[0].ix);

        let mut ls = Accum::new(x0, d);

        // Accumulate Least square matrices
        for o in obs {
            ls.accum(o);
        }
        // Get Cholesky decomposition of XWX (in place)
        let mut beta = [0.0; 3];
        ls.solve(&mut beta);
        Fit {
            x: x0 as isize,
            beta,
        }
    }

    fn pred(&self, pos: isize) -> Option<f64> {
        let x = (pos - self.x) as f64;
        let y = self.beta[0] + x * self.beta[1] + x * x * self.beta[2];
        if y < 1.0 {
            None
        } else {
            Some(y)
        }
    }
}

fn smooth(mut bc: Vec<Vec<usize>>) -> Vec<Option<f64>> {
    let n = bc.len();

    // Get median and weights (from inverse of estimated samples variance / n)
    let obs: Vec<_> = bc
        .iter_mut()
        .enumerate()
        .flat_map(|(ix, v)| Obs::new(ix, v))
        .collect();

    // Perform smoothing using a local quadratic function and a tricubic kernel

    // Number of points in smoothing region
    const REGION_SIZE: usize = 31;

    let region_size = obs.len().min(REGION_SIZE);
    // Can't fit local quadratic with less than 3 points...
    assert!(region_size > 2);

    // Fit local model with region_size points in each window
    // with the fitted point being as close as possible to the
    // centre of the window
    let mut left = 0;
    let mut right = region_size - 1;
    let l = obs.len();
    let mut fit = Vec::with_capacity(l);
    for i in 0..l {
        fit.push(Fit::fit_local_regression(&obs[left..=right], i - left));
        // Update window for next point
        if right - i - 1 < i + 1 - left && right < l - 1 {
            left += 1;
            right += 1;
        }
    }

    // Storage for predictions
    let mut pred = vec![None; n];

    pred[fit[0].x as usize] = Some(fit[0].beta[0]);
    for f in fit.windows(2) {
        for x in f[0].x + 1..=f[1].x {
            let k = if f[1].x - x > x - f[0].x { 0 } else { 1 };
            pred[x as usize] = f[k].pred(x)
        }
    }
    pred
}

/// Normalize coverage data for a sample based on GC content
/// This is done by getting the median coverage per GC bin from
/// contigs (normally the autosomes)
pub fn normalize_sample(cfg: &Config, rc: RawCounts) -> NormCov {
    // First collect counts per GC bin
    let mut bin_counts = collect_bin_data(cfg, rc);

    let pred = smooth(bin_counts);

    for (i, p) in pred.iter().enumerate() {
        println!("{}\t{:?}", i, p);
    }

    let mut nc = HashMap::new();
    nc
}
