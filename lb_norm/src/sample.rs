use std::{
    collections::{hash_map::Entry, HashMap, HashSet},
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
    output_flag: bool,
    control_flag: bool,
}

impl Sample {
    pub fn new(name: String) -> Self {
        Self {
            name,
            files: HashMap::new(),
            output_flag: false,
            control_flag: false,
        }
    }

    pub fn set_output_flag(&mut self) {
        self.output_flag = true;
    }

    pub fn set_control_flag(&mut self) {
        self.control_flag = true;
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn is_output(&self) -> bool {
        self.output_flag
    }
    pub fn is_control(&self) -> bool {
        self.control_flag
    }
    pub fn ctg_path(&self, ctg: &str) -> Option<&Path> {
        self.files.get(ctg).map(|x| x.as_path())
    }
}

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
        // Skip empty lines
        if !fields.is_empty() {
            sample_vec.push(Sample::new(fields[0].to_owned()));
        }
    }

    debug!(
        "Finished reading in {} lines; found {} samples",
        line,
        sample_vec.len()
    );

    Ok(sample_vec)
}

/// Add control samples to sample list, setting flags appropriately
pub fn merge_controls(samples: &mut Vec<Sample>, controls: &mut Vec<Sample>) {
    let mut h: HashMap<_, _> = samples
        .iter()
        .enumerate()
        .map(|(ix, s)| (s.name.clone(), ix))
        .collect();
    for mut c in controls.drain(..) {
        match h.entry(c.name.clone()) {
            Entry::Occupied(e) => {
                let i = *e.get();
                samples[i].set_control_flag();
            }
            Entry::Vacant(e) => {
                c.set_control_flag();
                let l = samples.len();
                samples.push(c);
                e.insert(l);
            }
        }
    }
}

/// Collect input file paths for each sample in samples.  
/// Each file path is parsed to extract the contig name.  
/// A vector of all contigs found is returned
pub fn get_input_files_and_contig_list(
    samples: &mut Vec<Sample>,
    dir: Option<&PathBuf>,
    prefix: &str,
) -> anyhow::Result<Vec<Contig>> {
    let mut ctg_hash = HashSet::new();
    let reg = Regex::new(format!("^{}_([^_]*)[.]txt$", prefix).as_str())?;
    for s in samples.iter_mut() {
        get_files_for_sample(s, dir, &reg, &mut ctg_hash)?
    }
    let v = ctg_hash.drain().collect();

    Ok(v)
}

fn get_files_for_sample(
    s: &mut Sample,
    dir: Option<&PathBuf>,
    reg: &Regex,
    ctg_hash: &mut HashSet<Contig>,
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
                if !ctg_hash.contains(ctg) {
                    trace!("Adding contig {}", ctg);
                    ctg_hash.insert(Arc::from(ctg));
                }
                let ctg = ctg_hash.get(ctg).unwrap().clone();
                trace!(
                    "Adding file {} ({}) for sample {}",
                    path.display(),
                    ctg,
                    name
                );
                s.files.insert(ctg, path);
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
