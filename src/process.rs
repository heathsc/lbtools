use std::{collections::HashMap, sync::Arc, thread, time::Duration};

use crossbeam_channel::{bounded, Receiver, Sender};
use r_htslib::*;

use crate::{config::Config, controller::*, input::open_input};

fn process_task(
    cfg: &Config,
    ix: usize,
    tpool: Option<&HtsThreadPool>,
    snd: Sender<JobRequest>,
    recv: Receiver<Option<Job>>,
) -> anyhow::Result<()> {
    debug!("Process task {} starting up", ix);
    let mut sample_idx = None;
    let mut hts = None;
    snd.send(JobRequest {
        prev_results: Completed::None,
        sample_idx,
        task_idx: ix,
    })?;

    while let Some(job) = recv.recv()? {
        trace!("Task {} received job {:?}", ix, job);
        let i = job.sample_idx;

        let res = match job.job_type {
            JobType::ReadData(ctg) => {
                // If this is a new sample, open the file
                if sample_idx.map(|x| x != i).unwrap_or(true) {
                    let fname = cfg.sample_list()[i].input_path();
                    trace!("Task {} opening file {}", ix, fname.display());
                    let mut h = open_input(fname, ctg.is_some(), cfg.reference(), tpool)?;
                    sample_idx = Some(i);
                    hts = Some(h);
                }

                // TODO
                // Make dummy counts for testing
                let mut h = HashMap::new();
                if let Some(c) = ctg {
                    let v: Vec<usize> = Vec::new();
                    h.insert(c, v);
                }
                Completed::RawCounts(i, h)
            }
            JobType::NormalizeSample(rc) => {
                let mut h = HashMap::new();
                for ctg in cfg.ctg_hash().keys() {
                    let v = Coverage::new(Vec::new());
                    h.insert(Arc::clone(ctg), v);
                }
                Completed::NormalizedCounts(i, h)
            }
            JobType::OutputSampleCtg(_) => Completed::None,
            JobType::Wait => Completed::None,
        };
        let d = Duration::from_secs(1);
        thread::sleep(d);
        snd.send(JobRequest {
            prev_results: res,
            sample_idx,
            task_idx: ix,
        })?;
    }
    debug!("Process task {} closing down", ix);
    Ok(())
}

/// Create child threads to process samples
pub fn process_samples(cfg: &Config) -> anyhow::Result<()> {
    // Set up Hts thread pool
    let tpool = HtsThreadPool::new(cfg.hts_threads());
    let tpool_ref = tpool.as_ref();

    thread::scope(|sc| {
        let nt = cfg.n_tasks();

        // Channel for a task to request a new job
        let (send_ctrl, recv_ctrl) = bounded(nt * 8);

        // Storage for channels by which a task receives a new job
        let mut send_job = Vec::with_capacity(nt);

        // Spawn task processes
        let mut join_handles: Vec<_> = (0..nt)
            .map(|ix| {
                let (s, r) = bounded(1);
                send_job.push(s);
                let s = send_ctrl.clone();
                sc.spawn(move || process_task(cfg, ix + 1, tpool_ref, s, r))
            })
            .collect();

        // Spawn controller processes
        let control_jh = sc.spawn(|| controller(cfg, recv_ctrl, send_job));
    });

    Ok(())
}
