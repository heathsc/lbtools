use crate::{config::Config, io, region::Region, sample::Sample};
use anyhow::Context;
use std::{
    fmt::{self, Formatter},
    io::Write,
    path::Path,
};

use compress_io::compress::CompressIo;

struct RegData<'a> {
    region: &'a Region,
    n: usize,
    mean: f64,
    sd: f64,
    sd_ratio: f64, // For confidence intervals
}

impl<'a> RegData<'a> {
    fn new(region: &'a Region, n: usize, mean: f64, sd: f64) -> Self {
        let sd_ratio = if let Some(d) = region.delta_cn() {
            let r = utils::qt(0.975, (n - 1) as f64)
                .expect("Error getting confidence intervals for ctDNA");
            if d.is_negative() {
                -r
            } else {
                r
            }
        } else {
            1.0
        };
        Self {
            region,
            n,
            mean,
            sd,
            sd_ratio,
        }
    }
}

struct SampleDataItem<'a> {
    reg_data: &'a RegData<'a>,
    copy_num: f64,
    t: f64,
    p: f64,
    q: Option<f64>,
    ct_dna: String,
}

impl<'a> SampleDataItem<'a> {
    fn new(reg_data: &'a RegData, z: f64) -> Self {
        let reg = reg_data.region;
        let mean = reg_data.mean;
        let sd = reg_data.sd;

        let copy_num = z - mean + 2.0;
        let diff = z - mean;
        let mut t = diff / sd;

        let ct_dna = reg
            .delta_cn()
            .map(|delta| {
                let sdr = reg_data.sd_ratio;
                if delta.is_negative() {
                    t = -t;
                };
                let d = delta as f64;
                let get_ct = |x: f64| (x / d).min(1.0).max(0.0);
                let (a, b, c) = (
                    get_ct(diff - sdr * sd),
                    get_ct(diff),
                    get_ct(diff + sdr * sd),
                );
                format!("{} ({}-{})", b, a, c)
            })
            .unwrap_or(String::from("NA"));

        let p = utils::pt(t, (reg_data.n - 1) as f64, false).expect("Error getting p-value");

        Self {
            reg_data,
            copy_num,
            t,
            p,
            q: None,
            ct_dna,
        }
    }
    fn add_q(&mut self, q: f64) {
        self.q = Some(q)
    }
}

impl<'a> fmt::Display for SampleDataItem<'a> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{:.6}\t{:.6}\t{}\t{:.6e}\t{:.6e}\t",
            self.reg_data.region.desc(),
            self.reg_data.n,
            self.reg_data.sd,
            self.copy_num,
            self.ct_dna,
            self.t,
            self.p
        )?;
        if let Some(q) = self.q {
            write!(f, "{:.6e}", q)
        } else {
            write!(f, "NA")
        }
    }
}

struct SampleData<'a> {
    sample: &'a Sample,
    data: Vec<SampleDataItem<'a>>,
}

impl<'a> SampleData<'a> {
    fn new(sample: &'a Sample) -> Self {
        Self {
            sample,
            data: Vec::new(),
        }
    }

    fn add_item(&mut self, item: SampleDataItem<'a>) {
        self.data.push(item)
    }
}

/// Strategy
///
/// Read in control data for each region so that we can get robust estimates of mean and sd
/// Then go through each test sample, getting the average coverage for each region
/// Perform t-test on each test sample and region, comparing against the control data
pub fn process_data(cfg: &Config) -> anyhow::Result<()> {
    debug!("Starting processing");

    let mut reg_data = Vec::with_capacity(cfg.regions().len());
    for reg in cfg.regions().iter() {
        debug!("Reading data for {}", reg.desc());
        let mut v = Vec::new();
        for s in cfg.sample_list().iter().filter(|x| x.is_control()) {
            if let Some(p) = s.ctg_path(reg.ctg()) {
                if let Some(z) = io::read_region_data(p, reg)? {
                    v.push(z)
                }
            }
        }
        v.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        let l = v.len();
        let q1 = v[l >> 2];
        let q2 = v[l >> 1];
        let q3 = v[(3 * l) >> 2];
        if let Some(sd) = utils::robust_sd(q3 - q1, l) {
            let mean = ((q1 + q2 + q3) / 3.0);
            debug!("n: {}, mean: {}, sd: {}", l, mean, sd);
            reg_data.push(RegData::new(reg, l, mean, sd))
        } else {
            warn!(
                "Not enough control data to get robust estimates of sdev for coverage of {}",
                reg.desc()
            )
        }
    }

    // Check test samples
    let mut sample_data = Vec::with_capacity(cfg.sample_list().len());
    for s in cfg.sample_list().iter().filter(|x| !x.is_control()) {
        let mut sdata = SampleData::new(s);
        for rdata in reg_data.iter() {
            let reg = rdata.region;
            if let Some(p) = s.ctg_path(reg.ctg()) {
                if let Some(z) = io::read_region_data(p, reg)? {
                    sdata.add_item(SampleDataItem::new(rdata, z));
                }
            }
        }
        sample_data.push(sdata);
    }

    // Correct p-values for multiple testing by controlling FDR
    // First, collect all p-values;
    let mut p = Vec::new();
    for sd in sample_data.iter() {
        for item in sd.data.iter() {
            p.push(item.p)
        }
    }
    debug!("Total number of p-values - {}", p.len());
    // get corrected vector
    let mut q = utils::fdr(&p);
    // Add corrected p-value for sample results
    let mut q_it = q.drain(..);
    for sd in sample_data.iter_mut() {
        for item in sd.data.iter_mut() {
            item.add_q(q_it.next().unwrap())
        }
    }
    // Write out results
    let mut wrt = CompressIo::new()
        .opt_path(cfg.output_file())
        .bufwriter()
        .with_context(|| "Failed to open output file")?;

    writeln!(
        wrt,
        "sample\tregion\tn\tsd\tcopy number\tctDNA\tt\tp\tp(FDR corrected)\t"
    )?;
    for sd in sample_data.iter() {
        for item in sd.data.iter() {
            writeln!(wrt, "{}\t{}", sd.sample.name(), item)?
        }
    }
    Ok(())
}
