use anyhow::Context;
use std::{collections::HashMap, sync::Arc, thread, time::Duration};

use crossbeam_channel::{bounded, Receiver, Sender};
use r_htslib::*;

use crate::{
    config::Config,
    controller::*,
    coverage::{Coverage, NormCov, RawCounts},
    input::open_input,
    normalize::normalize_sample,
    reader::read_coverage_data,
};

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
                    let h = open_input(fname, ctg.is_none(), cfg.reference(), tpool)?;
                    sample_idx = Some(i);
                    hts = Some(h);
                }
                debug!(
                    "Task {} reading ctg {:?} for sample {}",
                    ix,
                    ctg,
                    cfg.sample_list()[sample_idx.unwrap()].name()
                );
                let h = read_coverage_data(cfg, hts.as_mut().unwrap(), ctg.as_ref())?;
                Completed::RawCounts(i, h)
            }
            JobType::NormalizeSample(rc) => {
                debug!(
                    "Task {} normalizing sample {}",
                    ix,
                    cfg.sample_list()[i].name()
                );
                let h = normalize_sample(cfg, rc);
                Completed::NormalizedCounts(i, h)
            }
            JobType::OutputSampleCtg(_, _, _) => Completed::None,
            JobType::Wait => {
                let d = Duration::from_secs(5);
                thread::sleep(d);
                Completed::None
            }
        };
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
    debug!(
        "Setting up hts thread pool with {} threads",
        cfg.hts_threads()
    );
    let tpool = HtsThreadPool::new(cfg.hts_threads());
    let tpool_ref = tpool.as_ref();

    let mut res = Vec::new();
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

        // Spawn controller process
        let control_jh = sc.spawn(|| controller(cfg, recv_ctrl, send_job));
        drop(send_ctrl);

        // Join task processes
        for jh in join_handles.drain(..) {
            res.push(jh.join())
        }

        // Join controller process
        res.push(control_jh.join())
    });

    // Handle any errors from controller process
    if let Some(x) = res.pop() {
        match x {
            Ok(y) => y.with_context(|| "Error returned from controller process")?,
            Err(_) => return Err(anyhow!("Error joining controller process")),
        }
    }

    // Handler any errors from tasks
    for (i, x) in res.drain(..).enumerate() {
        match x {
            Ok(y) => y.with_context(|| format!("Error returned from task {}", i + 1))?,
            Err(_) => return Err(anyhow!("Error joining task {}", i + 1)),
        }
    }

    Ok(())
}
