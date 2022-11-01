use std::{collections::HashMap, sync::Arc};

use r_htslib::*;

use crate::{config::Config, coverage::*};

#[derive(Debug)]
struct ReadFilter {
    min_mapq: u8,
    min_len: usize,
    max_len: Option<usize>,
    keep_duplicates: bool,
    forbid_flags_paired: u16,
    forbid_flags_unpaired: u16,
}

const FORBID_FLAGS: u16 = BAM_FUNMAP | BAM_FSUPPLEMENTARY | BAM_FSECONDARY | BAM_FQCFAIL;

impl ReadFilter {
    fn new(cfg: &Config) -> Self {
        let mut forbid_flags_unpaired = FORBID_FLAGS;
        if !(cfg.ignore_dup_flag() | cfg.keep_duplicates()) {
            forbid_flags_unpaired |= BAM_FDUP
        }
        let forbid_flags_paired = forbid_flags_unpaired | BAM_FUNMAP;
        Self {
            min_mapq: cfg.min_mapq(),
            min_len: cfg.min_template_len(),
            max_len: cfg.max_template_len(),
            keep_duplicates: cfg.keep_duplicates(),
            forbid_flags_paired,
            forbid_flags_unpaired,
        }
    }

    fn pass_filter(&self, brec: &BamRec, prev_pos: &Option<(usize, usize, Option<usize>)>) -> bool {
        let flag = brec.flag();
        let mapq = brec.qual();
        if (flag & BAM_FPAIRED) == 0 {
            // Unpaired reads
            if mapq >= self.min_mapq && (flag & self.forbid_flags_unpaired) == 0 {
                // Check for duplicate (same coordinates as previous read)
                if let Some((tid, x, None)) = prev_pos {
                    brec.tid().unwrap() != *tid || brec.pos().unwrap() != *x
                } else {
                    true
                }
            } else {
                false
            }
        } else {
            // Paired reads
            if mapq >= self.min_mapq
                && (flag & (self.forbid_flags_paired | BAM_FPROPER_PAIR)) == BAM_FPROPER_PAIR
            {
                if !self.keep_duplicates {
                    // Check for duplicate (same coordinates as previous read)
                    if let Some((tid, x, Some(y))) = prev_pos {
                        if brec.tid().unwrap() == *tid
                            && brec.pos().unwrap() == *x
                            && brec.mpos().unwrap() == *y
                        {
                            return false;
                        }
                    }
                }
                let m = flag & (BAM_FREVERSE | BAM_FMREVERSE);
                (m == BAM_FREVERSE || m == BAM_FMREVERSE)
                    && if let Some(x) = self.max_len {
                        let l = brec.template_len().unsigned_abs();
                        l >= self.min_len && l <= x
                    } else if self.min_len > 0 {
                        let l = brec.template_len().unsigned_abs();
                        l >= self.min_len
                    } else {
                        true
                    }
            } else {
                false
            }
        }
    }
}

struct RawCounter {
    ctg: Arc<str>,
    cov: Vec<usize>,
    block_size: usize,
    seq_len: usize,
}

impl RawCounter {
    fn new(ctg: &Arc<str>, seq_len: usize, block_size: usize) -> Self {
        let n_bins = (seq_len + block_size - 1) / block_size;
        Self {
            ctg: Arc::clone(ctg),
            cov: vec![0; n_bins],
            block_size,
            seq_len,
        }
    }

    fn add_raw_counts(&mut self, rec: &BamRec, min_qual: u8) {
        let read_start = rec.pos().unwrap();
        let mut x = read_start;
        let end = rec.endpos() + 1;
        let flag = rec.flag();
        let mut y = if (flag & BAM_FPAIRED) != 0 {
            // Paired
            if (flag & BAM_FREVERSE) == 0 {
                // Forward (this should not panic because we are only looking at correctly mapped pairs
                let mate_pos = rec.mpos().unwrap();
                if mate_pos < x {
                    // The read pair overlaps such that the start of the + read is after the end of the - read
                    // We will handle these with the - read as it is simpler, so we just return for now
                    return;
                } else {
                    // We don't count the overlapping part (it will be counted when the - read is processed)
                    end.min(mate_pos)
                }
            } else {
                let mate_pos = rec.mpos().unwrap();
                // These are overlapping pairs where the end of the - read is before the start of the + read.
                // We ignored these when we encountered the + read, so we must deal with them here
                if mate_pos > x {
                    // Find the start of the - read
                    let x1 = rec.endpos();

                    if x1 < mate_pos {
                        // The pair does not overlap; all of the - read is to the left of the + read
                        // Normally such read pairs should be filtered before getting to this stage
                        return;
                    }
                    // We will count only the bases between the start of the + read and the start of the - read
                    x = mate_pos;
                    x1 + 1
                } else {
                    end
                }
            }
        } else {
            end
        };

        // Make sure we are not over the end of the contig
        y = y.min(self.seq_len);

        if y == x {
            return;
        }
        if y > x {
            // We will count bases in the interval [x, y)
            // We need to get the base qualities to apply the base quality filter
            if let Some(qv) = rec.get_qual() {
                let mut x1 = read_start;
                assert!(x1 <= x);
                for q in qv.iter() {
                    if x1 >= x && *q >= min_qual {
                        self.cov[x1 / self.block_size] += 1;
                    }
                    x1 += 1;
                    if x1 > y {
                        break;
                    }
                }
            }
        } else {
            warn!(
                "Read {} skipped due to inconsistent flags",
                rec.qname().unwrap()
            )
        }
    }
}

/// Read SAM/BAM/CRAM data from input file and calculate binned coverage
pub fn read_coverage_data(
    cfg: &Config,
    hts: &mut Hts,
    ctg: Option<&Arc<str>>,
) -> anyhow::Result<RawCounts> {
    if let Some(c) = ctg {
        read_ctg_coverage_data(cfg, hts, c)
    } else {
        read_sample_coverage_data(cfg, hts)
    }
}

/// Read data from a particular contig (requires indexed file)
fn read_ctg_coverage_data(
    cfg: &Config,
    hts: &mut Hts,
    ctg: &Arc<str>,
) -> anyhow::Result<RawCounts> {
    let mut rc = HashMap::new();
    if let Some(seq_len) = hts.seq_length(ctg) {
        let tid = hts.name2tid(ctg);
        let block_size = cfg.block_size() as usize;
        let filter = ReadFilter::new(cfg);
        trace!("Filter set to: {:?}", filter);
        let mut raw_cov = RawCounter::new(ctg, seq_len, block_size);
        let rlist = hts.make_region_list(&[ctg]);
        let mut rdr: HtsItrReader<BamRec> = hts.itr_reader(&rlist);
        let mut rec = BamRec::new()?;

        // Keep track of previous read so that we can remove duplicates if required
        let mut prev_pos: Option<(usize, usize, Option<usize>)> = None;
        while rdr.read(&mut rec)? {
            assert_eq!(rec.tid(), tid);
            if filter.pass_filter(&rec, &prev_pos) {
                raw_cov.add_raw_counts(&rec, cfg.min_qual());
                prev_pos = Some((rec.tid().unwrap(), rec.pos().unwrap(), rec.mpos()));
            }
        }
        rc.insert(raw_cov.ctg, raw_cov.cov);
    } else {
        warn!("Contig {} not found in input file", ctg);
    }

    Ok(rc)
}

/// Read data for all requested contigs from file without index
fn read_sample_coverage_data(cfg: &Config, hts: &mut Hts) -> anyhow::Result<RawCounts> {
    let mut rc = HashMap::new();
    let mut rec = BamRec::new()?;
    let filter = ReadFilter::new(cfg);

    // Construct hash with keys being the tid of the required sequences and the
    // values being RawCounter structures
    let mut chash: HashMap<_, _> = cfg
        .ctg_hash()
        .keys()
        .filter_map(|ctg| {
            hts.seq_length(ctg).map(|l| {
                (
                    hts.name2tid(ctg).unwrap(),
                    RawCounter::new(ctg, l, cfg.block_size() as usize),
                )
            })
        })
        .collect();

    // Keep track of previous read so that we can remove duplicates if required
    // Only works if input is sorted on genomic order
    let mut prev_pos: Option<(usize, usize, Option<usize>)> = None;

    while rec.read(hts)? {
        if filter.pass_filter(&rec, &prev_pos) {
            if let Some(raw_cov) = chash.get_mut(&rec.tid().unwrap()) {
                raw_cov.add_raw_counts(&rec, cfg.min_qual());
                prev_pos = Some((rec.tid().unwrap(), rec.pos().unwrap(), rec.mpos()));
            }
        }
    }

    for (_, raw_cov) in chash.drain() {
        rc.insert(raw_cov.ctg, raw_cov.cov);
    }
    Ok(rc)
}
