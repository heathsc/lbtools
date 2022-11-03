use std::{collections::HashSet, path::Path, sync::Arc};

use crate::config::Contig;

use anyhow::Context;
use compress_io::compress::CompressIo;
use log::Level::Debug;
use regex::Regex;
use utils::get_next_line;

#[derive(Debug)]
pub struct Region {
    // Description for the region
    desc: String,
    // Expected change in copy-number, i.e., for a single copy deletion this should be Some(-1)
    // Use None where the expected change is not known
    delta_cn: Option<isize>,

    ctg: Contig,
    // Sorted, Non-overlapping vector of ranges.  Each range is a tuple (a,b) where b>=a
    ranges: Vec<(usize, usize)>,
}

impl Region {
    pub fn ranges(&self) -> &[(usize, usize)] {
        &self.ranges
    }
    pub fn desc(&self) -> &str {
        &self.desc
    }
    pub fn ctg(&self) -> &Contig {
        &self.ctg
    }
    pub fn delta_cn(&self) -> Option<isize> {
        self.delta_cn
    }
}

fn parse_usize_with_commas(s: &str) -> anyhow::Result<usize> {
    s.replace(',', "")
        .parse::<usize>()
        .with_context(|| "Error parsing coordinate")
}

fn parse_single_range(s1: &str, s2: &str) -> anyhow::Result<(usize, usize)> {
    let a = parse_usize_with_commas(s1)?;
    let b = parse_usize_with_commas(s2)?;
    if b >= a {
        Ok((a, b))
    } else {
        Err(anyhow!("Range error - {} > {}", s1, s2))
    }
}

/// Parse a comma delimited set of ranges (i.e., 10-50000, 60000-65000)
/// and make a sort list of non-overlapping tuples (a,b) where b >= a
fn parse_ranges(s: &str) -> anyhow::Result<Vec<(usize, usize)>> {
    let r = Regex::new(r"^\s*([0-9,]+)\s*[-:]\s*([0-9,]+)\s*$").unwrap();
    let mut v: Vec<(usize, usize)> = Vec::new();
    for s1 in s.split(',') {
        if let Some(c) = r.captures(s1) {
            v.push(
                parse_single_range(c.get(1).unwrap().as_str(), c.get(2).unwrap().as_str())
                    .with_context(|| format!("Illegal range: {}", s1))?,
            );
        } else {
            return Err(anyhow!("Illegal range: {}", s1));
        }
    }
    // Check that range(s) are valid
    if v.is_empty() {
        Err(anyhow!("No valid ranges found in {}", s))
    } else if v.len() > 1 {
        // Form sorted, non-overlapping vector from v
        v.sort_unstable_by_key(|(a, _)| *a);
        let mut v1 = Vec::new();
        let mut prev = v[0];
        for (a, b) in &v[1..] {
            // Check for overlap
            if *a <= prev.1 {
                prev.1 = prev.1.max(*b)
            } else {
                v1.push(prev);
                prev = (*a, *b)
            }
        }
        v1.push(prev);
        Ok(v1)
    } else {
        Ok(v)
    }
}

pub fn read_region_file<P: AsRef<Path>>(
    fname: P,
) -> anyhow::Result<(Vec<Region>, HashSet<Contig>)> {
    debug!("Reading in region list from {}", fname.as_ref().display());

    trace!("Opening region file for reading");
    let mut rdr = CompressIo::new().path(&fname).bufreader()?;

    trace!("Reading from file");
    let mut buf = String::new();
    let mut line = 0;
    let mut region_vec = Vec::new();
    let mut ctg_hash = HashSet::new();
    while let Some(fields) = get_next_line(&mut rdr, &mut buf).with_context(|| {
        format!(
            "Error after reading {} lines from {}",
            line,
            fname.as_ref().display()
        )
    })? {
        line += 1;
        // Parse input line and store if valid
        // Skip short lines
        if fields.len() >= 3 {
            let ranges = parse_ranges(fields[2]).with_context(|| {
                format!(
                    "{}:{} Parse error for sample {}",
                    fname.as_ref().display(),
                    line,
                    fields[0]
                )
            })?;
            let delta_cn = fields.get(3).and_then(|s| s.parse::<isize>().ok());
            if !ctg_hash.contains(fields[1]) {
                trace!("Adding contig {}", fields[1]);
                ctg_hash.insert(Arc::from(fields[1]));
            }
            let ctg = ctg_hash.get(fields[1]).unwrap().clone();
            let desc = fields[0].to_owned();
            region_vec.push(Region {
                desc,
                delta_cn,
                ctg,
                ranges,
            })
        }
    }
    if log_enabled!(Debug) {
        debug!("Regions:");
        for r in region_vec.iter() {
            debug!("\t{:?}", r)
        }
    }
    Ok((region_vec, ctg_hash))
}
