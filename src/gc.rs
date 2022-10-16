use std::{collections::HashMap, io::BufRead, path::Path, sync::Arc, thread};

use anyhow::Context;
use compress_io::compress::CompressIo;
use crossbeam_channel::{unbounded, Receiver};
use r_htslib::Faidx;

use crate::contig::Contig;

pub const N_GC_BINS: u32 = 100;
const MIN_GC_COUNT: u32 = (0.9 * (N_GC_BINS as f64)) as u32;

const MTAB: [usize; 256] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

struct GcBuilder {
    ctg: Arc<str>,
    data: Vec<Option<u32>>,
    counts: [u32; 3],
    block_size: u32,
    current_pos: usize,
}

impl GcBuilder {
    // Returns the bin corresponding to a set of counts
    fn bin(&self) -> Option<u32> {
        let tot = self.counts[1] + self.counts[2];
        if tot >= MIN_GC_COUNT {
            Some(
                (((self.counts[2] as f64 / tot as f64) * (N_GC_BINS as f64)).floor() as u32)
                    .min(N_GC_BINS - 1),
            )
        } else {
            None
        }
    }

    fn new(ctg: &Arc<str>, block_size: u32) -> Self {
        Self {
            ctg: Arc::clone(ctg),
            data: Vec::new(),
            counts: [0; 3],
            current_pos: 0,
            block_size,
        }
    }

    fn add_u8(&mut self, c: u8) {
        let ix = self.current_pos / (self.block_size as usize);
        if ix != self.data.len() {
            self.update_vec();
        }
        self.counts[MTAB[c as usize]] += 1;
        self.current_pos += 1;
    }

    fn add_str(&mut self, s: &str) {
        for c in s.as_bytes() {
            self.add_u8(*c)
        }
    }

    fn update_vec(&mut self) {
        self.data.push(self.bin());
        self.counts = [0; 3];
    }
}

pub struct GcCtgData {
    name: Arc<str>,
    data: Vec<Option<u32>>,
    block_size: u32,
}

impl GcCtgData {
    pub fn gc_bin(&self, x: usize) -> Option<u32> {
        let ix = x / (self.block_size as usize);
        self.data.get(ix).and_then(|x| *x)
    }

    fn from_builder(mut gcb: GcBuilder) -> Self {
        gcb.update_vec();
        Self {
            name: gcb.ctg,
            data: gcb.data,
            block_size: gcb.block_size,
        }
    }
}

/// Read in next line.
fn get_next_line<R: BufRead>(rdr: &mut R, buf: &mut String) -> anyhow::Result<bool> {
    buf.clear();
    if rdr.read_line(buf)? == 0 {
        Ok(false)
    } else {
        Ok(true)
    }
}

pub struct GcData {
    chash: HashMap<Arc<str>, GcCtgData>,
}

impl GcData {
    pub fn ctg_data(&self, ctg: &str) -> Option<&GcCtgData> {
        self.chash.get(ctg)
    }

    pub fn from_reference<S: AsRef<Path>>(
        fname: S,
        block_size: u32,
        nt: usize,
        ctg_hash: &HashMap<Arc<str>, Contig>,
    ) -> anyhow::Result<Self> {
        debug!(
            "Reading reference sequence from {} and calculating gc bins with block size {}",
            fname.as_ref().display(),
            block_size
        );

        if nt == 1 {
            single_threaded_read(fname, block_size, ctg_hash)
        } else {
            // Check if the reference has an index
            trace!("Test for faidx index");
            match Faidx::load(&fname) {
                Ok(_) => {
                    trace!("Index found: use multithreaded reading");
                    multi_threaded_read(fname, block_size, nt, ctg_hash)
                }
                Err(e) => {
                    trace!("Couldn't open file for indexed reading: {}", e);
                    single_threaded_read(fname, block_size, ctg_hash)
                }
            }
        }
    }
}

fn multi_threaded_read<S: AsRef<Path>>(
    fname: S,
    block_size: u32,
    nt: usize,
    ctg_hash: &HashMap<Arc<str>, Contig>,
) -> anyhow::Result<GcData> {
    let fname = fname.as_ref();
    let mut v = Vec::with_capacity(nt);
    // Everything runs within a scope so that we can pass references to the threads
    thread::scope(|sc| {
        // Create channels to send jobs the threads
        trace!(
            "Spawning {} readers for reference file {}",
            nt,
            fname.display()
        );

        // Spawn reader threads
        let (snd, rcv) = unbounded();
        let jobs: Vec<_> = (0..nt)
            .map(|i| {
                let r = rcv.clone();
                sc.spawn(move || reader(fname, block_size, i + 1, r))
            })
            .collect();
        drop(rcv);

        // Send required contigs to child threads
        for ctg in ctg_hash.keys() {
            if snd.send(ctg).is_err() {
                error!("Error sending message to child readers");
                break;
            }
        }

        drop(snd);
        for jh in jobs {
            v.push(jh.join())
        }
    });

    trace!("Collecting results from child threads");
    let mut chash = HashMap::new();
    for (ix, ch) in v.drain(..).enumerate() {
        match ch {
            Ok(c) => {
                let mut h =
                    c.with_context(|| format!("Error returned from GC read thread {}", ix + 1))?;
                for (k, v) in h.drain() {
                    chash.insert(k, v);
                }
            }
            Err(_) => return Err(anyhow!("Error joining GC read thread {}", ix + 1)),
        }
    }

    debug!("Finished reading reference and calculating gc bins");
    let tst = chash.get("chr1").unwrap();
    for (i, k) in tst.data.iter().enumerate() {
        println!("{}\t{:?}", i * (block_size as usize), k);
    }
    Ok(GcData { chash })
}

fn reader(
    fname: &Path,
    block_size: u32,
    ix: usize,
    r: Receiver<&Arc<str>>,
) -> anyhow::Result<HashMap<Arc<str>, GcCtgData>> {
    trace!("Starting up GC reader thread {}", ix);
    let faidx =
        Faidx::load(fname).with_context(|| format!("Error opening file {}", fname.display()))?;
    let mut chash = HashMap::new();
    while let Ok(ctg) = r.recv() {
        trace!("GC reader {} processing contig {}", ix, ctg);
        let mut gcb = GcBuilder::new(ctg, block_size);
        let s = faidx
            .fetch_seq(ctg, 0, None)
            .with_context(|| format!("Error fetching sequence for contig {}", ctg))?;
        for c in s.seq().iter() {
            gcb.add_u8(*c)
        }
        store_ctg_data(gcb, &mut chash);
        trace!("GC reader {} finished processing contig {}", ix, ctg);
    }
    trace!("Closing down GC reader thread {}", ix);
    Ok(chash)
}

fn single_threaded_read<S: AsRef<Path>>(
    fname: S,
    block_size: u32,
    ctg_hash: &HashMap<Arc<str>, Contig>,
) -> anyhow::Result<GcData> {
    trace!("Opening reference file for reading");
    let mut rdr = CompressIo::new()
        .path(&fname)
        .bufreader()
        .with_context(|| format!("Error opening reference file {}", fname.as_ref().display()))?;

    trace!("Reading from reference file");
    let mut buf = String::new();
    let mut line = 0;
    let mut chash = HashMap::new();
    let mut gcb: Option<GcBuilder> = None;
    while get_next_line(&mut rdr, &mut buf).with_context(|| {
        format!(
            "Error after reading {} lines from {}",
            line,
            fname.as_ref().display()
        )
    })? {
        line += 1;
        if buf.starts_with('>') {
            // New contig
            if let Some(ctg) = buf.trim_start_matches('>').split_ascii_whitespace().next() {
                if let Some(b) = gcb.take() {
                    store_ctg_data(b, &mut chash)
                }
                if let Some((k, _)) = ctg_hash.get_key_value(ctg) {
                    gcb = Some(GcBuilder::new(k, block_size));
                    trace!("Processing ctg {}", ctg);
                }
            } else {
                return Err(anyhow!("Missing contig name at line {}", line));
            }
        } else if let Some(b) = gcb.as_mut() {
            b.add_str(buf.trim_end())
        }
    }
    if let Some(b) = gcb.take() {
        store_ctg_data(b, &mut chash)
    }
    debug!("Finished reading reference and calculating gc bins");
    let tst = chash.get("chr1").unwrap();
    for (i, k) in tst.data.iter().enumerate() {
        println!("{}\t{:?}", i * (block_size as usize), k);
    }
    Ok(GcData { chash })
}

fn store_ctg_data(b: GcBuilder, chash: &mut HashMap<Arc<str>, GcCtgData>) {
    let ctg_data = GcCtgData::from_builder(b);
    let c = Arc::clone(&ctg_data.name);
    trace!("Storing gc data for contig {}", c);
    chash.insert(c, ctg_data);
}
