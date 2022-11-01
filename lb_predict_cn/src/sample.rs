use std::path::{Path, PathBuf};

use anyhow::Context;
use compress_io::compress::CompressIo;

use crate::utils::get_next_line;

/// Input sample
///
/// name - used to generate output file
/// input_path - path to input SAM/BAM/CRAM file
///
pub struct Sample {
    name: String,
    input_path: PathBuf,
}

impl Sample {
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn input_path(&self) -> &Path {
        self.input_path.as_ref()
    }
}

/// Read in sample list from file
/// Expects two tab separated columns.  
/// The first column has the sample name (used for the output files)
/// The second column has the path to the SAM/BAM/CRAM file for this sample
///
pub fn sample_vec_from_file<S: AsRef<Path>>(fname: S) -> anyhow::Result<Vec<Sample>> {
    debug!("Reading in sample list from {}", fname.as_ref().display());

    trace!("Opening sample file for reading");
    let mut rdr = CompressIo::new()
        .path(&fname)
        .bufreader()
        .with_context(|| format!("Error opening contig file {}", fname.as_ref().display()))?;

    trace!("Reading from sample file");
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

        // Parse input line and store to hash if valid
        if fields.len() >= 2 {
            // Skip short lines
            let sample = Sample {
                name: fields[0].to_owned(),
                input_path: PathBuf::from(fields[1]),
            };
            trace!(
                "Read in sample {} path {}",
                sample.name,
                sample.input_path.display()
            );
            sample_vec.push(sample)
        }
    }

    debug!(
        "Finished reading in {} lines; found {} samples",
        line,
        sample_vec.len()
    );
    Ok(sample_vec)
}
