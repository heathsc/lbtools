use std::path::Path;

use anyhow::Context;
use compress_io::compress::CompressIo;
use utils::get_next_line;

use crate::region::Region;

pub fn read_region_data(p: &Path, reg: &Region) -> anyhow::Result<Option<f64>> {
    trace!("Opening sample file {} for reading", p.display());
    let mut rdr = CompressIo::new().path(p).bufreader()?;
    trace!("Reading from {}", p.display());
    let mut buf = String::new();
    let mut line = 0;

    // List of sorted, non-overlapping ranges for this region
    let ranges = reg.ranges();
    // Get end coordinate of last range (so we know when we should stop reading)
    let last_x = ranges.last().unwrap().1;

    // Read in pos, cn pairs
    let mut v = Vec::new();
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
            v.push((x, z));
            if x > last_x {
                break;
            }
        }
    }

    if v.is_empty() {
        return Ok(None);
    }
    // Guess bin spacing
    let bin_size: usize = v.windows(2).map(|x| x[1].0 - x[0].0).min().unwrap_or(1);
    trace!("bin_size = {}", bin_size);
    let mut s = 0.0;
    let mut n: usize = 0;

    for (x, z) in v.drain(..) {
        // Get limits of bin
        let (lo, hi) = (x - (bin_size >> 1), x + (bin_size >> 1));

        // Check if x lies within a range
        if ranges.iter().any(|(a, b)| hi > *a && lo <= *b) {
            n += 1;
            s += z;
        }
    }
    Ok(if n > 0 { Some(s / (n as f64)) } else { None })
}
