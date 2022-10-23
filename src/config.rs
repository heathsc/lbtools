use std::{
    collections::HashMap,
    path::{Path, PathBuf},
    sync::Arc,
};

use crate::{contig::Contig, gc::GcData, sample::Sample};

/// Config
///
/// Configuration info for the program
/// This is generated from the command line arguments
/// Once set it is read only
///
/// sample_list - list of input samples
/// ctg_hash - look up hash for contigs
/// gc_data - data on gc content per bin across the genome
/// reference - path to reference (FASTA) file
/// output_dir - output directory
/// block_size - block size for coverage
/// min_template_len - minimum allowed template (fragment) length
/// max_template_len - maximum allowed template length
/// threads - number of threads
///
pub struct Config {
    sample_list: Vec<Sample>,
    ctg_hash: HashMap<Arc<str>, Contig>,
    gc_data: GcData,
    reference: PathBuf,
    output_dir: Option<PathBuf>,
    block_size: u32,
    min_template_len: usize,
    max_template_len: Option<usize>,
    min_mapq: u8,
    min_qual: u8,
    keep_duplicates: bool,
    ignore_dup_flag: bool,
    output_prefix: String,
    hts_threads: usize,
    n_tasks: usize,
    n_readers: usize,
}

impl Config {
    pub fn new(
        sample_list: Vec<Sample>,
        ctg_hash: HashMap<Arc<str>, Contig>,
        gc_data: GcData,
        reference: PathBuf,
        output_prefix: String,
    ) -> Self {
        Self {
            sample_list,
            ctg_hash,
            gc_data,
            reference,
            output_prefix,
            output_dir: None,
            block_size: 1000,
            min_template_len: 0,
            max_template_len: None,
            keep_duplicates: false,
            ignore_dup_flag: false,
            min_mapq: 0,
            min_qual: 0,
            hts_threads: 1,
            n_tasks: 1,
            n_readers: 1,
        }
    }

    pub fn set_output_dir<P: AsRef<Path>>(&mut self, dir: P) {
        self.output_dir = Some(dir.as_ref().to_owned())
    }

    pub fn set_block_size(&mut self, bs: u32) {
        self.block_size = bs
    }

    pub fn set_min_template_len(&mut self, x: usize) -> anyhow::Result<()> {
        self.min_template_len = x;
        match self.max_template_len {
            Some(y) if y < x => Err(anyhow!("Invalid template lengths - maximum < minimum")),
            _ => Ok(()),
        }
    }

    pub fn set_min_mapq(&mut self, x: u8) {
        self.min_mapq = x;
    }

    pub fn set_min_qual(&mut self, x: u8) {
        self.min_qual = x;
    }

    pub fn set_keep_duplicates(&mut self) {
        self.keep_duplicates = true
    }

    pub fn set_ignore_dup_flag(&mut self) {
        self.ignore_dup_flag = true
    }

    pub fn set_max_template_len(&mut self, x: usize) -> anyhow::Result<()> {
        self.max_template_len = Some(x);
        if x < self.min_template_len {
            Err(anyhow!("Invalid template lengths - maximum < minimum"))
        } else {
            Ok(())
        }
    }

    pub fn set_hts_threads(&mut self, x: usize) {
        self.hts_threads = x
    }

    pub fn set_n_tasks(&mut self, x: usize) {
        self.n_tasks = x
    }

    pub fn set_n_readers(&mut self, x: usize) {
        self.n_readers = x
    }

    pub fn sample_list(&self) -> &[Sample] {
        &self.sample_list
    }

    pub fn ctg_hash(&self) -> &HashMap<Arc<str>, Contig> {
        &self.ctg_hash
    }

    pub fn gc_data(&self) -> &GcData {
        &self.gc_data
    }

    pub fn reference(&self) -> &Path {
        self.reference.as_ref()
    }

    pub fn output_prefix(&self) -> &str {
        &self.output_prefix
    }

    pub fn output_dir(&self) -> Option<&Path> {
        self.output_dir.as_deref()
    }

    pub fn block_size(&self) -> u32 {
        self.block_size
    }

    pub fn min_template_len(&self) -> usize {
        self.min_template_len
    }

    pub fn max_template_len(&self) -> Option<usize> {
        self.max_template_len
    }

    pub fn min_mapq(&self) -> u8 {
        self.min_mapq
    }

    pub fn min_qual(&self) -> u8 {
        self.min_qual
    }

    pub fn hts_threads(&self) -> usize {
        self.hts_threads
    }

    pub fn n_tasks(&self) -> usize {
        self.n_tasks
    }

    pub fn n_readers(&self) -> usize {
        self.n_readers
    }

    pub fn keep_duplicates(&self) -> bool {
        self.keep_duplicates
    }

    pub fn ignore_dup_flag(&self) -> bool {
        self.ignore_dup_flag
    }
}
