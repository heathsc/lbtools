use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::sample::Sample;

pub struct Config {
    sample_list: Vec<Sample>,
    ctg_list: Vec<Contig>,
    output_prefix: String,
    output_dir: Option<PathBuf>,
}

impl Config {
    pub fn new(output_prefix: String, sample_list: Vec<Sample>, ctg_list: Vec<Contig>) -> Self {
        Self {
            sample_list,
            ctg_list,
            output_prefix,
            output_dir: None,
        }
    }
    pub fn set_output_dir(&mut self, d: PathBuf) {
        self.output_dir = Some(d)
    }

    pub fn ctg_list(&self) -> &[Contig] {
        &self.ctg_list
    }

    pub fn sample_list(&self) -> &[Sample] {
        &self.sample_list
    }

    pub fn output_dir(&self) -> Option<&Path> {
        self.output_dir.as_deref()
    }

    pub fn output_prefix(&self) -> &str {
        &self.output_prefix
    }
}

pub type Contig = Arc<str>;
