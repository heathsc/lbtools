use std::path::PathBuf;

/// Input sample
///
/// name - used to generate output file
/// input_path - path to input SAM/BAM/CRAM file
///
pub struct Sample {
    name: String,
    input_path: PathBuf,
}
