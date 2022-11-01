use crate::{config::Config, io};
use std::collections::BTreeMap;
use std::path::PathBuf;

/// Strategy
///
/// Read in control data contig by contig
/// When complete for a contig, generate median vector
/// We can then start the output of normalized data, first for the
/// controls samples that are also output samples, and then for
/// the rest of the samples
pub fn process_samples(cfg: &Config) -> anyhow::Result<()> {
    debug!("Starting processing");
    let in_dir = cfg
        .output_dir()
        .map(|p| p.to_owned())
        .unwrap_or_else(PathBuf::new);

    for ctg in cfg.ctg_list().iter() {
        debug!("Reading data from {}", ctg);
        let mut bt = BTreeMap::new();
        for s in cfg.sample_list().iter().filter(|x| x.is_control()) {
            if let Some(p) = s.ctg_path(ctg) {
                let mut v = io::read_sample_contig_data(p)?;
                for (i, x) in v.drain(..) {
                    let e = bt.entry(i).or_insert_with(|| Vec::new());
                    e.push(x)
                }
            }
        }

        debug!("Calculate median vector for {}", ctg);
        let mut iqr = Vec::with_capacity(bt.len());
        let med: BTreeMap<_, _> = bt
            .iter_mut()
            .map(|(i, v)| {
                v.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
                let l = v.len();
                let q1 = v[l >> 2];
                let q2 = v[l >> 1];
                let q3 = v[(3 * l) >> 2];
                iqr.push(q3 - q1);
                (*i, (q2, q3 - q1))
            })
            .collect();

        // We will exclude the top and bottom 0.5% of bins depending in IQR
        iqr.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        let l = iqr.len() as f64;
        let low = iqr[(l * 0.005) as usize];
        let high = iqr[(l * 0.995) as usize];
        debug!("Output normalized data for {}", ctg);

        for s in cfg.sample_list().iter().filter(|x| x.is_output()) {
            if let Some(p) = s.ctg_path(ctg) {
                let mut opath = in_dir.clone();
                opath.push(s.name());
                let oname = format!("{}_{}.txt", cfg.output_prefix(), ctg);
                opath.push(&oname);
                io::output_sample_contig_data(&p, &opath, &med, low, high)?;
            }
        }
    }
    Ok(())
}
