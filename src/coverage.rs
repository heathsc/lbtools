use std::{collections::HashMap, sync::Arc};

pub type Coverage = Vec<(usize, Option<f64>)>;
pub type RawCounts = HashMap<Arc<str>, Vec<usize>>;
pub type NormCov = HashMap<Arc<str>, Coverage>;
