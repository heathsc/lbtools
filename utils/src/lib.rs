#[macro_use]
extern crate anyhow;

use std::{fmt, io::BufRead, str::FromStr};

use clap::ArgMatches;
use special::Beta;

/// LogLevel
///
/// Represents minimum level of messages that will be logged
///
#[derive(Debug, Clone, Copy)]
pub struct LogLevel {
    pub level: usize,
}

impl FromStr for LogLevel {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "error" => Ok(LogLevel { level: 0 }),
            "warn" => Ok(LogLevel { level: 1 }),
            "info" => Ok(LogLevel { level: 2 }),
            "debug" => Ok(LogLevel { level: 3 }),
            "trace" => Ok(LogLevel { level: 4 }),
            "none" => Ok(LogLevel { level: 5 }),
            _ => Err("no match"),
        }
    }
}

impl LogLevel {
    pub fn is_none(&self) -> bool {
        self.level > 4
    }
    pub fn get_level(&self) -> usize {
        if self.level > 4 {
            0
        } else {
            self.level
        }
    }
}

impl fmt::Display for LogLevel {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let level_str = ["error", "warn", "info", "debug", "trace", "none"];
        if self.level < 6 {
            write!(f, "{}", level_str[self.level])
        } else {
            write!(f, "unknown")
        }
    }
}

/// Initialize logging from command line arguments
pub fn init_log(m: &ArgMatches) {
    let verbose = m
        .get_one::<LogLevel>("loglevel")
        .copied()
        .unwrap_or_else(|| LogLevel::from_str("info").expect("Could not set loglevel info"));
    let quiet = verbose.is_none() || m.get_flag("quiet");
    let ts = m
        .get_one::<stderrlog::Timestamp>("timestamp")
        .copied()
        .unwrap_or(stderrlog::Timestamp::Off);

    stderrlog::new()
        .quiet(quiet)
        .verbosity(verbose.get_level())
        .timestamp(ts)
        .init()
        .unwrap();
}

/// Read in next line and split on tabs after trimming white space
pub fn get_next_line<'a, R: BufRead>(
    rdr: &mut R,
    buf: &'a mut String,
) -> anyhow::Result<Option<Vec<&'a str>>> {
    buf.clear();
    if rdr.read_line(buf)? == 0 {
        Ok(None)
    } else {
        Ok(Some(buf.trim().split('\t').collect()))
    }
}

/// Robust estimation of sd from IQR following approach of
/// Wan et al. (2014) doi: 10.1186/1471-2288-14-135
///
/// Table 2 from Wan et al.

const SD_TAB: &[f64] = &[
    0.990, 1.144, 1.206, 1.239, 1.260, 1.274, 1.284, 1.292, 1.298, 1.303, 1.307, 1.311, 1.313,
    1.316, 1.318, 1.320, 1.322, 1.323, 1.324, 1.326, 1.327, 1.328, 1.329, 1.330, 1.330, 1.331,
    1.332, 1.332, 1.333, 1.333, 1.334, 1.334, 1.335, 1.335, 1.336, 1.336, 1.336, 1.337, 1.337,
    1.337, 1.338, 1.338, 1.338, 1.338, 1.339, 1.339, 1.339, 1.339, 1.339, 1.340,
];

pub fn robust_sd(iqr: f64, n: usize) -> Option<f64> {
    if n > 4 {
        let q = (n - 1) >> 2;
        let z = SD_TAB.get(q - 1).copied().unwrap_or_else(|| {
            let zmax = *SD_TAB.last().unwrap();
            let nn = n as f64;
            (2.0 * qnorm((0.75 * nn - 0.125) / (nn + 0.25)).unwrap()).max(zmax)
        });
        Some(iqr / z)
    } else {
        None
    }
}

/// Percentage points of normal distribution using a Rust translation of
/// Wichura, M. J. (1988) Algorithm AS 241: The percentage points of
/// the normal distribution.  _Applied Statistics_, *37*, 477-484.
/// This is the same source used by R.

/// Required constants
const QN_COEFF: [[f64; 15]; 3] = [
    [
        3.3871328727963666080E0,
        1.3314166789178437745E2,
        1.9715909503065514427E3,
        1.3731693765509461125E4,
        4.5921953931549871457E4,
        6.7265770927008700853E4,
        3.3430575583588128105E4,
        2.5090809287301226727E3,
        4.2313330701600911252E1,
        6.8718700749205790830E2,
        5.3941960214247511077E3,
        2.1213794301586595867E4,
        3.9307895800092710610E4,
        2.8729085735721942674E4,
        5.2264952788528545610E3,
    ],
    [
        1.42343711074968357734E0,
        4.63033784615654529590E0,
        5.76949722146069140550E0,
        3.64784832476320460504E0,
        1.27045825245236838258E0,
        2.41780725177450611770E-1,
        2.27238449892691845833E-2,
        7.74545014278341407640E-4,
        2.05319162663775882187E0,
        1.67638483018380384940E0,
        6.89767334985100004550E-1,
        1.48103976427480074590E-1,
        1.51986665636164571966E-2,
        5.47593808499534494600E-4,
        1.05075007164441684324E-9,
    ],
    [
        6.65790464350110377720E0,
        5.46378491116411436990E0,
        1.78482653991729133580E0,
        2.96560571828504891230E-1,
        2.65321895265761230930E-2,
        1.24266094738807843860E-3,
        2.71155556874348757815E-5,
        2.01033439929228813265E-7,
        5.99832206555887937690E-1,
        1.36929880922735805310E-1,
        1.48753612908506148525E-2,
        7.86869131145613259100E-4,
        1.84631831751005468180E-5,
        1.42151175831644588870E-7,
        2.04426310338993978564E-15,
    ],
];

const SPLIT1: f64 = 0.425E0;
const SPLIT2: f64 = 5.0E0;
const CONST1: f64 = 0.180625E0;
const CONST2: f64 = 1.6E0;

pub fn qnorm(p: f64) -> anyhow::Result<f64> {
    let f = |z: f64, v: &[f64; 15]| {
        (((((((v[7] * z + v[6]) * z + v[5]) * z + v[4]) * z + v[3]) * z + v[2]) * z + v[1]) * z
            + v[0])
            / (((((((v[14] * z + v[13]) * z + v[12]) * z + v[11]) * z + v[10]) * z + v[9]) * z
                + v[8])
                * z
                + 1.0)
    };

    let q = p - 0.5;
    if q.abs() <= SPLIT1 {
        let r = CONST1 - q * q;
        Ok(q * f(r, &QN_COEFF[0]))
    } else {
        let r = if q < 0.0 { p } else { 1.0 - p };
        if r <= 0.0 {
            return Err(anyhow!("qnorm(): Invalid parameter {}", p));
        }
        let r = (-(r.ln())).sqrt();
        let z = if r <= SPLIT2 {
            f(r - CONST2, &QN_COEFF[1])
        } else {
            f(r - SPLIT2, &QN_COEFF[2])
        };
        Ok(if q < 0.0 { -z } else { z })
    }
}

pub fn pt(t: f64, df: f64, lower_tail: bool) -> anyhow::Result<f64> {
    if df <= 0.0 {
        Err(anyhow!("pt(): Invalid df {}", df))
    } else {
        let a = df * 0.5;
        let lbeta = a.ln_beta(0.5);
        let x = df / (df + t * t);
        let z = 0.5 * x.inc_beta(a, 0.5, lbeta);
        let flip = lower_tail ^ t.is_sign_negative();
        Ok(if flip { 1.0 - z } else { z })
    }
}

pub fn qt(p: f64, df: f64) -> anyhow::Result<f64> {
    if df <= 0.0 {
        Err(anyhow!("qt(): Invalid df {}", df))
    } else if !(0.0..=1.0).contains(&p) {
        Err(anyhow!("qt(): p out of range: {}", p))
    } else if p == 0.0 {
        Ok(f64::NEG_INFINITY)
    } else if p == 1.0 {
        Ok(f64::INFINITY)
    } else {
        let a = df * 0.5;
        let lbeta = a.ln_beta(0.5);
        let (x, flip) = if p < 0.5 {
            (2.0 * p, false)
        } else {
            (2.0 * (1.0 - p), true)
        };
        let z = x.inv_inc_beta(a, 0.5, lbeta);
        let t = (df / z - df).sqrt();
        Ok(if flip { t } else { -t })
    }
}

/// Perform multiple test correction for a p value vector using the FDR method
/// of Benjamini & Hochberg (1995).
pub fn fdr(p: &[f64]) -> Vec<f64> {
    fdr_n(p, p.len())
}

pub fn fdr_n(p: &[f64], n: usize) -> Vec<f64> {
    let mut v: Vec<_> = p.iter().enumerate().collect();
    v.sort_unstable_by(|a, b| a.1.partial_cmp(b.1).unwrap());
    let n = n as f64;
    let mut min_p: f64 = 1.0;
    let mut q = vec![0.0; p.len()];
    for (i, (k, p)) in v.iter().enumerate().rev() {
        min_p = min_p.min((n / ((i + 1) as f64)) * *p);
        q[*k] = min_p;
    }
    q
}
