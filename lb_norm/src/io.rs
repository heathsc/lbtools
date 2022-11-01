use std::{collections::BTreeMap, io::Write, path::Path};

use anyhow::Context;
use compress_io::compress::CompressIo;
use utils::get_next_line;

pub fn read_sample_contig_data(p: &Path) -> anyhow::Result<Vec<(usize, f64)>> {
    let mut v = Vec::new();
    trace!("Opening sample file {} for reading", p.display());
    let mut rdr = CompressIo::new().path(p).bufreader()?;
    trace!("Reading from {}", p.display());
    let mut buf = String::new();
    let mut line = 0;

    while let Some(fields) = get_next_line(&mut rdr, &mut buf)
        .with_context(|| format!("Error after reading {} lines from {}", line, p.display()))?
    {
        line += 1;
        // Parse input line and store to hash if valid
        if fields.len() >= 3 {
            // Skip short lines
            let x = fields[1]
                .parse::<usize>()
                .with_context(|| format!("{}:{} Error reading position", p.display(), line))?;
            let z = fields[2]
                .parse::<f64>()
                .with_context(|| format!("{}:{} Error reading copy number", p.display(), line))?;
            v.push((x, z))
        }
    }
    Ok(v)
}

pub fn output_sample_contig_data(
    in_path: &Path,
    out_path: &Path,
    med: &BTreeMap<usize, (f64, f64)>,
    low: f64,
    high: f64,
) -> anyhow::Result<()> {
    trace!(
        "Opening sample file {} for reading; corrected data will be written to {}",
        in_path.display(),
        out_path.display()
    );
    let mut rdr = CompressIo::new().path(in_path).bufreader()?;
    let mut wrt = CompressIo::new().path(out_path).bufwriter()?;
    trace!(
        "Reading from {} and writing to {}",
        in_path.display(),
        out_path.display()
    );
    let mut buf = String::new();
    let mut line = 0;

    while let Some(fields) = get_next_line(&mut rdr, &mut buf).with_context(|| {
        format!(
            "Error after reading {} lines from {}",
            line,
            in_path.display()
        )
    })? {
        line += 1;
        // Parse input line and store to hash if valid
        if fields.len() >= 3 {
            // Skip short lines
            let x = fields[1].parse::<usize>().with_context(|| {
                format!("{}:{} Error reading position", in_path.display(), line)
            })?;
            if let Some((m, iqr)) = med.get(&x) {
                if *iqr > low && *iqr < high {
                    let z = fields[2].parse::<f64>().with_context(|| {
                        format!("{}:{} Error reading copy number", in_path.display(), line)
                    })?;
                    let y = 2.0 + z - m;
                    writeln!(wrt, "{}\t{}\t{:.4}\t{}", fields[0], x, y, fields[2]).with_context(
                        || format!("Error writing corrected data to {}", out_path.display()),
                    )?;
                }
            }
        }
    }
    Ok(())
}
