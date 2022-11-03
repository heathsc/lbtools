use std::{
    collections::{HashMap, HashSet},
    path::{Path, PathBuf},
    sync::Arc,
};

use anyhow::Context;
use compress_io::compress::CompressIo;
use regex::Regex;
use utils::get_next_line;

use crate::config::Contig;

pub struct Sample {
    name: String,
    files: HashMap<Contig, PathBuf>,
    control_flag: bool,
}

impl Sample {
    pub fn new(name: String, control_flag: bool) -> Self {
        Self {
            name,
            files: HashMap::new(),
            control_flag,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn is_control(&self) -> bool {
        self.control_flag
    }
    pub fn ctg_path(&self, ctg: &str) -> Option<&Path> {
        self.files.get(ctg).map(|x| x.as_path())
    }
}

/// Reading sample list file
/// Each line should contain the sample name and an indicator of whether this is
/// a test or control sample.  The indicator is the word "test" or "control"
/// (case insensitive), or any prefix of the two words (i.e., 'T' or 'C' works.
pub fn read_sample_list_from_file<P: AsRef<Path>>(fname: P) -> anyhow::Result<Vec<Sample>> {
    debug!("Reading in sample list from {}", fname.as_ref().display());

    trace!("Opening sample file for reading");
    let mut rdr = CompressIo::new().path(&fname).bufreader()?;

    trace!("Reading from file");
    let mut buf = String::new();
    let mut line = 0;
    let mut sample_vec = Vec::new();

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
        if fields.len() >= 2 {
            let stype = parse_test_control(fields[1]).with_context(|| {
                format!(
                    "{}:{} Parse error for sample {}",
                    fname.as_ref().display(),
                    line,
                    fields[0]
                )
            })?;

            sample_vec.push(Sample::new(fields[0].to_owned(), stype));
        }
    }

    debug!(
        "Finished reading in {} lines; found {} samples",
        line,
        sample_vec.len()
    );

    Ok(sample_vec)
}

/// Checks s to see if it is a prefix of 'Test' or 'Control' (case insensitive)
/// Returns Ok(true) if s is a prefix of 'Control', Ok(false) if a prefix of 'Test',
/// and Err otherwise.
fn parse_test_control(s: &str) -> anyhow::Result<bool> {
    let s = s.to_lowercase();
    if "test".starts_with(&s) {
        Ok(false)
    } else if "control".starts_with(&s) {
        Ok(true)
    } else {
        Err(anyhow!("Could not parse sample type {}", s))
    }
}

/// Collect input file paths for each sample in samples.  
/// Each file path is parsed to extract the contig name, and this
/// is checked to see if it exists in ctg_hash, and matching files are stored
pub fn get_input_files_and_contig_list(
    samples: &mut Vec<Sample>,
    dir: Option<&PathBuf>,
    prefix: &str,
    ctg_hash: &HashSet<Contig>,
) -> anyhow::Result<()> {
    let reg = Regex::new(format!("^{}_([^_]*)[.]txt$", prefix).as_str())?;
    for s in samples.iter_mut() {
        get_files_for_sample(s, dir, &reg, ctg_hash)?
    }

    Ok(())
}

fn get_files_for_sample(
    s: &mut Sample,
    dir: Option<&PathBuf>,
    reg: &Regex,
    ctg_hash: &HashSet<Contig>,
) -> anyhow::Result<()> {
    let mut in_dir = dir.map(|p| p.to_owned()).unwrap_or_else(PathBuf::new);
    in_dir.push(&s.name);

    for f in in_dir
        .read_dir()
        .with_context(|| format!("Error checking input directory {}", in_dir.display()))?
    {
        let entry =
            f.with_context(|| format!("Could not get directory entry from {}", in_dir.display()))?;
        let path = entry.path();
        if path.is_file() {
            let name = entry.file_name().into_string().expect("Illegal file name");
            if let Some(c) = reg.captures(name.as_str()) {
                let ctg = c.get(1).unwrap().as_str();
                if let Some(c) = ctg_hash.get(ctg) {
                    trace!(
                        "Adding file {} ({}) for sample {}",
                        path.display(),
                        ctg,
                        name
                    );
                    s.files.insert(c.clone(), path);
                }
            }
        }
    }

    if s.files.is_empty() {
        Err(anyhow!(
            "No input files found for sample {} in {}",
            s.name,
            in_dir.display()
        ))
    } else {
        debug!(
            "{} input files found for sample {} in {}",
            s.files.len(),
            s.name,
            in_dir.display()
        );
        Ok(())
    }
}
