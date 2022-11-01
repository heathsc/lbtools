/// Coordination of jobs between tasks
///
/// The pool of process tasks (threads) request jobs from and return processed
/// results to the controller.  The tasks return any results from a previous job (if existing)
/// and request a new job.  If no more jobs are available then None is returned, and the child
/// tasks will exit.
///
/// Tasks are loosely aligned to sample files to avoid excessive opening/closing of input files.
/// This is achieved by checking first if a read job is available for the same sample that the
/// child task was previously reading.  If not then an alternate sample will be assigned to the
/// task.
///
/// Possible job types are:
///
///   ReadData - read a contig (if file is indexed) or all contigs for a sample
///   NormalizeSample - Perform GC normalization on all contigs of a sample
///   OutputSampleCtg - Output a contig for a processed sample
///   Wait - No jobs are available, but more will be available in future
///
///   After processing a ReadData job the child tasks will return Completed::RawCounts.  When
///   all contigs have been processed for a sample it will be eligible for Normalization.
///
///   Processing of a NormalizeSample jobs will result in Completed::NormalizeCounts. As this
///   is for a complete sample this will be immediately eligible for Output
///
///   Processing of an output job has no results returned (just a request for a new job)
///
use std::{collections::hash_map, fmt, sync::Arc};

use anyhow::Context;
use crossbeam_channel::{Receiver, Sender};
use r_htslib::*;

use crate::{
    config::Config,
    coverage::{Coverage, NormCov, RawCounts},
    sample::Sample,
};

pub enum JobType {
    ReadData(Option<Arc<str>>),
    NormalizeSample(RawCounts),
    OutputSampleCtg(usize, Arc<str>, Coverage),
    Wait, // No jobs currently available, but there will be jobs in the future
}

impl fmt::Debug for JobType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ReadData(s) => write!(f, "JobType::ReadData({:?})", s),
            Self::NormalizeSample(_) => f.write_str("JobType::NormalizeSample"),
            Self::OutputSampleCtg(i, s, _) => {
                write!(f, "JobType::OutputSampleCtg({}, {:?})", *i, s)
            }
            Self::Wait => f.write_str("JobType::Wait"),
        }
    }
}

/// The child tasks send their results as Completed objects
pub enum Completed {
    RawCounts(usize, RawCounts), // (sample id, raw (un-normalized) counts)
    NormalizedCounts(usize, NormCov), // (sample id, normalized and raw counts
    None, // This is returned either initially or after a task receives a Wait or OutputSampleCtg job
}

impl fmt::Debug for Completed {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RawCounts(i, _) => write!(f, "Completed::RawCounts(Sample {})", *i),
            Self::NormalizedCounts(i, _) => write!(f, "Completed::NormalizedCounts(Sample {})", *i),
            Self::None => f.write_str("Completed::None"),
        }
    }
}

/// Sent from child tasks to request a new job
#[derive(Debug)]
pub struct JobRequest {
    pub prev_results: Completed, // Returned results from previous job by this task
    pub sample_idx: Option<usize>, // The current sample the task has been reading from (if any)
    pub task_idx: usize, // Id of task (used to select the channel to send the reply back to)
}

/// Sent to child task in response to a JobRequest
#[derive(Debug)]
pub struct Job {
    pub sample_idx: usize, // The sample that will be processed
    pub job_type: JobType, // The type of job
}

/// Keep track of pending Jobs (those hat have been sent out and the results have not yet come back)
#[derive(Default, Debug)]
struct Tracker {
    n_read_jobs_pending: usize,
    n_normalize_jobs_pending: usize,
}

impl Tracker {
    fn update_at_send(&mut self, job: &Job) {
        match job.job_type {
            JobType::ReadData(_) => self.n_read_jobs_pending += 1,
            JobType::NormalizeSample(_) => self.n_normalize_jobs_pending += 1,
            _ => (),
        }
    }

    fn update_at_recv(&mut self, jr: &JobRequest) {
        match jr.prev_results {
            Completed::RawCounts(_, _) => {
                assert!(self.n_read_jobs_pending > 0);
                self.n_read_jobs_pending -= 1;
            }
            Completed::NormalizedCounts(_, _) => {
                assert!(self.n_normalize_jobs_pending > 0);
                self.n_normalize_jobs_pending -= 1;
            }
            _ => (),
        }
    }

    fn pending(&self) -> bool {
        self.n_read_jobs_pending > 0 || self.n_normalize_jobs_pending > 0
    }
}

/// A sample that is currently being output
struct OnGoingOutput {
    sample_idx: usize,
    norm_cov: Vec<(Arc<str>, Coverage)>,
}

impl OnGoingOutput {
    fn new(sample_idx: usize, mut nc: NormCov) -> Self {
        trace!("OngoingOutput::new({})", sample_idx);
        let norm_cov: Vec<_> = nc.drain().collect();
        Self {
            sample_idx,
            norm_cov,
        }
    }

    fn next_job(&mut self) -> Option<Job> {
        trace!("OngoingOutput::next_job({})", self.sample_idx);
        self.norm_cov.pop().map(|(ctg, c)| Job {
            sample_idx: self.sample_idx,
            job_type: JobType::OutputSampleCtg(self.sample_idx, ctg, c),
        })
    }
}

/// An input file: keeps track of which contigs remain to be read
struct InputFile<'a, T> {
    sample: &'a Sample,
    sample_idx: usize,
    ctg_iter: hash_map::Keys<'a, Arc<str>, T>,
    indexed: Option<bool>,
    finished: bool,
}

impl<'a, T> InputFile<'a, T> {
    fn new(
        sample_idx: usize,
        sample: &'a Sample,
        ctg_iter: hash_map::Keys<'a, Arc<str>, T>,
    ) -> Self {
        Self {
            sample,
            sample_idx,
            ctg_iter,
            indexed: None,
            finished: false,
        }
    }

    fn check_finished(&mut self) -> anyhow::Result<bool> {
        if self.indexed.is_none() {
            let path = self.sample.input_path();
            trace!(
                "Collector checking input file {} for sample {}",
                self.sample.name(),
                path.display()
            );
            let mut hts = Hts::open(path, "r").with_context(|| {
                format!(
                    "Error opening input file {} for sample {}",
                    path.display(),
                    self.sample.name()
                )
            })?;
            // Check that this is a SAM type file (SAM/BAM/CRAM)
            if !matches!(hts.rec_type(), Some(HtsRecType::Sam)) {
                return Err(anyhow!(
                    "Incorrect file format for input file {}",
                    path.display()
                ));
            }
            self.indexed = Some(hts.index_load().is_ok());
        }
        Ok(self.finished)
    }

    fn next_job(&mut self) -> Option<Job> {
        if self.finished {
            None
        } else if self.indexed.unwrap() {
            match self.ctg_iter.next() {
                Some(c) => Some(Job {
                    sample_idx: self.sample_idx,
                    job_type: JobType::ReadData(Some(Arc::clone(c))),
                }),
                None => {
                    self.finished = true;
                    None
                }
            }
        } else {
            self.finished = true;
            Some(Job {
                sample_idx: self.sample_idx,
                job_type: JobType::ReadData(None),
            })
        }
    }
}

/// Selects an InputFile with pending contigs from sample_vec.  
/// Starts looking at index idx and processed through the whole vector,
/// wrapping around if required. On return idx will be set to the next
/// index after the selected sample (if the selection is made)
fn get_new_read_job<'a, T>(
    sample_vec: &mut [InputFile<'a, T>],
    idx: &mut usize,
) -> anyhow::Result<Option<Job>> {
    // Find first sample with available samples starting from *idx (and wrapping around)
    let l = sample_vec.len();
    for _ in 0..l {
        if !sample_vec[*idx].check_finished()? {
            break;
        }
        *idx = (*idx + 1) % l;
    }

    // Get job (if available) from selected sample
    let job = sample_vec[*idx].next_job();
    *idx = (*idx + 1) % l;
    Ok(job)
}

/// Main loop.  Recieves messages from child tasks and allocates jobs appropriately.  Will
/// end if channel r is closed (i.e., when all child tasks exit) or on error
pub fn controller(
    cfg: &Config,
    r: Receiver<JobRequest>,
    svec: Vec<Sender<Option<Job>>>,
) -> anyhow::Result<()> {
    debug!("Controller thread starting up");

    let ns = cfg.sample_list().len();
    let nc = cfg.ctg_hash().len();
    let mut track = Tracker::default();

    // Tracking for samples/ctgs to be read
    let mut sample_vec: Vec<_> = cfg
        .sample_list()
        .iter()
        .enumerate()
        .map(|(i, s)| InputFile::new(i, s, cfg.ctg_hash().keys()))
        .collect();
    assert!(!sample_vec.is_empty());
    let mut sample_idx = 0;

    let read_job_limit = cfg.n_readers();

    // Tracking for samples to be normalized
    let mut sample_data: Vec<Option<RawCounts>> = vec![None; ns];
    let mut pending_norm: Vec<(usize, RawCounts)> = Vec::new();

    // Tracking for samples/ctgs still to be output
    let mut pending_output: Vec<(usize, NormCov)> = Vec::new();
    let mut ongoing_output: Option<OnGoingOutput> = None;

    while let Ok(jr) = r.recv() {
        trace!("Controller received request {:?}; pending: {:?}", jr, track);

        track.update_at_recv(&jr);

        // Store data from previous results
        match jr.prev_results {
            // Returning raw counts.  Add to sample_data
            Completed::RawCounts(i, mut h) => {
                let cts = if let Some(mut d) = sample_data[i].take() {
                    for (k, v) in h.drain() {
                        d.insert(k, v);
                    }
                    d
                } else {
                    h
                };

                // If all contigs have been read then move to pending_norm else store in sample_data
                if cts.len() == nc {
                    pending_norm.push((i, cts))
                } else {
                    sample_data[i] = Some(cts)
                }
            }
            Completed::NormalizedCounts(i, v) => pending_output.push((i, v)),
            Completed::None => (),
        }

        // See if we can add new read jobs
        let new_reads = track.n_read_jobs_pending < read_job_limit;

        // First we check if we have more contigs to read from the requested sample
        let mut job = if new_reads {
            jr.sample_idx.and_then(|i| sample_vec[i].next_job())
        } else {
            None
        };

        // If not, we see if we are already in the process of outputting a sample
        job = job
            .or_else(|| ongoing_output.as_mut().and_then(|o| o.next_job()))
            // Otherwise, check if there is a complete sample ready to start outputting
            .or_else(|| {
                ongoing_output = pending_output
                    .pop()
                    .map(|(ix, nc)| OnGoingOutput::new(ix, nc));
                ongoing_output.as_mut().and_then(|o| o.next_job())
            })
            // If we have no available output jobs, check if there is a normalization jobs waiting
            .or_else(|| {
                pending_norm.pop().map(|(ix, c)| Job {
                    sample_idx: ix,
                    job_type: JobType::NormalizeSample(c),
                })
            });

        // If we still have no job, check for additional sample/ctgs for reading
        if job.is_none() {
            job = if new_reads {
                get_new_read_job(
                    // &mut ongoing_read,
                    // &mut pending_read,
                    //cfg.ctg_hash(),
                    &mut sample_vec,
                    &mut sample_idx,
                )?
            } else {
                None
            }
            // If we get here then we have no pending read, normalization or output jobs
            // Check if jobs have been sent for processing that have not returned.  If yes,
            // then return JobType::Wait otherwise processing is finished so we can return None.
            .or_else(|| {
                if track.pending() {
                    Some(Job {
                        sample_idx: 0,
                        job_type: JobType::Wait,
                    })
                } else {
                    None
                }
            });
        };

        if let Some(j) = job.as_ref() {
            track.update_at_send(j)
        }

        trace!(
            "Controller sending back job {:?} for task {}",
            job,
            jr.task_idx
        );
        svec[jr.task_idx - 1]
            .send(job)
            .expect("Error sending message to task");
    }
    debug!("Controller thread closing down");
    Ok(())
}
