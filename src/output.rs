use anyhow::Context;
use std::{
    fs,
    io::{BufWriter, Write},
    path::PathBuf,
};

use crate::{config::Config, coverage::Coverage};

fn get_file_path(cfg: &Config, sample_idx: usize, ctg: &str) -> PathBuf {
    let mut p = if let Some(d) = cfg.output_dir() {
        d.to_owned()
    } else {
        PathBuf::new()
    };
    p.push(cfg.sample_list()[sample_idx].name());
    let name = format!("{}_{}.txt", cfg.output_prefix(), ctg);
    p.push(&name);
    p
}

pub fn setup_output(cfg: &Config) -> anyhow::Result<()> {
    // Create output directories
    let p = if let Some(d) = cfg.output_dir() {
        if !d.exists() {
            fs::create_dir_all(d)
                .with_context(|| format!("Error creating output directory {}", d.display()))?;
        }
        d.to_owned()
    } else {
        PathBuf::new()
    };
    for s in cfg.sample_list() {
        let p1 = p.join(s.name());
        if !p1.exists() {
            fs::create_dir(&p1)
                .with_context(|| format!("Error creating output directory {}", p1.display()))?;
        }
    }
    Ok(())
}

pub fn output_sample_cfg(
    cfg: &Config,
    sample_idx: usize,
    ctg: &str,
    mut cov: Coverage,
) -> anyhow::Result<()> {
    let opath = get_file_path(cfg, sample_idx, ctg);
    let mut wrt = BufWriter::new(
        fs::File::create(&opath)
            .with_context(|| format!("problem creating output file {}", opath.display()))?,
    );
    let bs = cfg.block_size() as f64;
    for (i, (rc, norm)) in cov.drain(..).enumerate() {
        let x = (((i as f64) + 0.5) * bs).round() as usize;
        writeln!(wrt, "{}\t{}\t{:.4}\t{:.4}", ctg, x, norm, (rc as f64) / bs)?
    }
    Ok(())
}
