use std::{
    path::{Path, PathBuf},
    sync::Arc,
};

use crate::{region::Region, sample::Sample};

pub struct Config {
    sample_list: Vec<Sample>,
    ctg_list: Vec<Contig>,
    regions: Vec<Region>,
    output_file: Option<PathBuf>,
}

impl Config {
    pub fn new(
        sample_list: Vec<Sample>,
        ctg_list: Vec<Contig>,
        regions: Vec<Region>,
        output_file: Option<PathBuf>,
    ) -> Self {
        Self {
            sample_list,
            ctg_list,
            regions,
            output_file,
        }
    }

    pub fn ctg_list(&self) -> &[Contig] {
        &self.ctg_list
    }

    pub fn sample_list(&self) -> &[Sample] {
        &self.sample_list
    }

    pub fn regions(&self) -> &[Region] {
        &self.regions
    }

    pub fn output_file(&self) -> Option<&Path> {
        self.output_file.as_deref()
    }
}

pub type Contig = Arc<str>;
