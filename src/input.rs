use std::path::Path;

use anyhow::Context;
use r_htslib::*;

pub fn open_input<P: AsRef<Path>>(
    name: P,
    no_index: bool,
    reference: &Path,
    tpool: Option<&HtsThreadPool>,
) -> anyhow::Result<Hts> {
    let name = name.as_ref();
    debug!(
        "Try to open input file {} with reference {:?} and no_index option {}",
        name.display(),
        reference,
        no_index
    );

    // Set up htsFormat so we can add the reference
    let mut fmt = HtsFormat::default();
    fmt.opt_add(format!("reference={}", reference.display()))?;

    // Try opening input file
    let mut hts = Hts::open_format(name, "r", &fmt)
        .with_context(|| format!("Failed to open input file {}", name.display()))?;

    // Check that this is a SAM type file (SAM/BAM/CRAM)
    if !matches!(hts.rec_type(), Some(HtsRecType::Sam)) {
        Err(anyhow!(
            "Incorrect file format for input file {}",
            name.display()
        ))
    } else {
        // Try and load index
        if !no_index {
            hts.index_load()
                .with_context(|| format!("Index not found for input file {}", name.display()))?;
        }

        // Try and attach thread pool if presemt
        if let Some(tp) = tpool {
            debug!("Attach thread pool to file");
            hts.hts_file_mut().set_thread_pool(tp)?
        }
        Ok(hts)
    }
}
