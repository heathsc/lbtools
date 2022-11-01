use std::{collections::HashMap, path::Path, sync::Arc};

use anyhow::Context;
use compress_io::compress::CompressIo;

use utils::get_next_line;

/// Contig
///
/// name - this is shared across many data structures so we use Rc<str>
/// use_for_normalization - whether to use this contig for sample normalization (normally set for the autosomes)
///
pub struct Contig {
    name: Arc<str>,
    use_for_normalization: bool,
}

impl Contig {
    fn new(name: &str, flag: bool) -> Self {
        trace!("Creating new contig {}, normalize: {}", name, flag);
        Self {
            name: Arc::from(name.to_owned()),
            use_for_normalization: flag,
        }
    }

    fn add_to_hash(self, hash: &mut HashMap<Arc<str>, Contig>) {
        let name = Arc::clone(&self.name);
        hash.insert(name, self);
    }

    pub fn name(&self) -> &Arc<str> {
        &self.name
    }

    pub fn use_for_normalization(&self) -> bool {
        self.use_for_normalization
    }
}

fn parse_bool(s: &str) -> anyhow::Result<bool> {
    match s.to_ascii_lowercase().as_str() {
        "true" | "1" | "yes" => Ok(true),
        "false" | "0" | "no" => Ok(false),
        _ => Err(anyhow!("Could not parse {} as bool", s)),
    }
}

/// Read in contig list from file
/// Expects one or two tab separated columns.  
/// The first column has the contig name
/// The second column, if present, should be 0/no/false or 1/yes/true to indicate
/// whether or not the contig should be used for normalization.  If absent, ttue is assumed.
///
pub fn contig_hash_from_file<S: AsRef<Path>>(
    fname: S,
) -> anyhow::Result<HashMap<Arc<str>, Contig>> {
    debug!("Reading in contig list from {}", fname.as_ref().display());

    trace!("Opening contig file for reading");
    let mut rdr = CompressIo::new()
        .path(&fname)
        .bufreader()
        .with_context(|| format!("Error opening contig file {}", fname.as_ref().display()))?;

    trace!("Reading from contig file");
    let mut buf = String::new();
    let mut line = 0;
    let mut ctg_hash = HashMap::new();

    while let Some(fields) = get_next_line(&mut rdr, &mut buf).with_context(|| {
        format!(
            "Error after reading {} lines from {}",
            line,
            fname.as_ref().display()
        )
    })? {
        line += 1;

        // Parse input line and store to hash if valid
        if !fields.is_empty() {
            // Skip blank lines
            let flag = match fields.get(1) {
                Some(s) => parse_bool(*s)
                    .with_context(|| format!("Error at {}:{}", fname.as_ref().display(), line))?,
                None => true,
            };
            Contig::new(fields[0], flag).add_to_hash(&mut ctg_hash)
        }
    }

    debug!(
        "Finished reading in {} lines; found {} contigs",
        line,
        ctg_hash.len()
    );
    Ok(ctg_hash)
}
