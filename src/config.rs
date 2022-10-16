use std::{collections::HashMap, path::PathBuf, sync::Arc};

use crate::{contig::Contig, sample::Sample};

/// Config
///
/// Configuration info for the program
/// This is generated from the command line arguments
/// Once set it is read only
///
/// sample_list - list of input samples
/// ctg_list - look up hash for contigs
/// reference - path to reference (FASTA) file
/// output_dir - output directory
/// block_size - block size for coverage
/// min_template_len - minimum allowed template (fragment) length
/// max_template_len - maximum allowed template length
/// threads - number of threads
///
pub struct Config {
    sample_list: Vec<Sample>,
    ctg_list: HashMap<Arc<str>, Contig>,
    reference: PathBuf,
    output_dir: PathBuf,
    block_size: usize,
    min_template_len: usize,
    max_template_len: Option<usize>,
    output_prefix: String,
    threads: usize,
}

impl Config {
    pub fn new() -> Self {
        Self {
            sample_list: Vec::new(),
            ctg_list: HashMap::new(),
            reference: Default::default(),
            output_dir: Default::default(),
            block_size: 0,
            min_template_len: 0,
            max_template_len: None,
            output_prefix: "".to_string(),
            threads: 0,
        }
    }
}
