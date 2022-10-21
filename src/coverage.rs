use std::{collections::HashMap, sync::Arc};

pub struct Coverage {
    v: Vec<(usize, f64)>,
}

impl Coverage {
    pub fn new(v: Vec<(usize, f64)>) -> Self {
        Self { v }
    }
}

pub type RawCounts = HashMap<Arc<str>, Vec<usize>>;
pub type NormCov = HashMap<Arc<str>, Coverage>;
