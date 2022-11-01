use std::{num::NonZeroUsize, path::PathBuf};

use clap::{
    crate_authors, crate_description, crate_name, crate_version, value_parser, Arg, ArgAction,
    Command,
};

use anyhow::Context;

use utils::{init_log, LogLevel};

use crate::{config::*, sample::*};

/// Set up definition of command options for clap
fn cli_model() -> Command {
    Command::new(crate_name!())
        .about(crate_description!())
        .version(crate_version!())
        .author(crate_authors!())
        .arg(
            Arg::new("timestamp")
                .short('X')
                .long("timestamp")
                .value_parser(value_parser!(stderrlog::Timestamp))
                .value_name("GRANULARITY")
                .default_value("none")
                .help("Prepend log entries with a timestamp"),
        )
        .arg(
            Arg::new("loglevel")
                .short('l')
                .long("loglevel")
                .value_name("LOGLEVEL")
                .value_parser(value_parser!(LogLevel))
                .ignore_case(true)
                .default_value("warn")
                .help("Set log level"),
        )
        .arg(
            Arg::new("quiet")
                .action(ArgAction::SetTrue)
                .long("quiet")
                .conflicts_with("loglevel")
                .help("Silence all output"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .value_parser(value_parser!(NonZeroUsize))
                .value_name("INT")
                .help("Set number of calculation threads [default: available cores]"),
        )
        .arg(
            Arg::new("input_prefix")
                .short('P')
                .long("input-prefix")
                .value_parser(value_parser!(String))
                .value_name("STRING")
                .default_value("cov")
                .help("Set prefix for input file names"),
        )
        .arg(
            Arg::new("output_prefix")
                .short('p')
                .long("output-prefix")
                .value_parser(value_parser!(String))
                .value_name("STRING")
                .default_value("ncov")
                .help("Set prefix for output file names"),
        )
        .arg(
            Arg::new("input_dir")
                .short('D')
                .long("input-dir")
                .value_parser(value_parser!(PathBuf))
                .value_name("PATH")
                .help("Set input directory [default: current directory]"),
        )
        .arg(
            Arg::new("output_dir")
                .short('d')
                .long("output-dir")
                .value_parser(value_parser!(PathBuf))
                .value_name("PATH")
                .help("Set output directory [default: current directory]"),
        )
        .arg(
            Arg::new("control_list")
                .short('c')
                .long("control-list")
                .value_parser(value_parser!(PathBuf))
                .value_name("PATH")
                .help("FIle with list of sample names to be used as controls [default: use samples from sample_list]"),
        )
        .arg(
            Arg::new("sample_list")
                .value_parser(value_parser!(PathBuf))
                .value_name("SAMPLE_FILE")
                .required(true)
                .help("Input file with list of sample names for normalization"),
        )
}

/// Handle command line options.  Set up Config structure
pub fn handle_cli() -> anyhow::Result<Config> {
    // Get matches from command line
    let m = cli_model().get_matches();

    // Setup logging
    init_log(&m);

    debug!("Processing command line options");

    let nt = m
        .get_one::<NonZeroUsize>("threads")
        .map(|x| usize::from(*x))
        .unwrap_or_else(num_cpus::get);

    let output_prefix = m
        .get_one::<String>("output_prefix")
        .expect("Missing default output prefix")
        .clone();

    let input_prefix = m
        .get_one::<String>("input_prefix")
        .expect("Missing default output prefix")
        .clone();

    // Read in sample list
    let mut samples = read_sample_list_from_file(
        m.get_one::<PathBuf>("sample_list")
            .expect("Missing sample list file"),
    )
    .with_context(|| "Could not open sample list file for input")?;
    for s in samples.iter_mut() {
        s.set_output_flag()
    }

    let input_dir = m.get_one::<PathBuf>("input_dir");

    // Read in control list if present and merge with sample list
    if let Some(cfile) = m.get_one::<PathBuf>("control_list") {
        let mut controls = read_sample_list_from_file(cfile)
            .with_context(|| "Could not open control list file for input")?;
        // Add controls to samples vec where they don't already exist
        merge_controls(&mut samples, &mut controls);
    } else {
        // If no control list specified, set all samples to be controls as well
        for s in samples.iter_mut() {
            s.set_control_flag();
        }
    }

    let contigs = get_input_files_and_contig_list(&mut samples, input_dir, &input_prefix)
        .with_context(|| "Error collecting input files")?;

    debug!("Number of contigs found: {}", contigs.len());

    let mut cfg = Config::new(output_prefix, samples, contigs);

    if let Some(p) = m.get_one::<PathBuf>("output_dir") {
        cfg.set_output_dir(p.to_owned())
    }

    cfg.set_threads(nt);

    // Make sure output does not overlap input
    if cfg.output_prefix() == input_prefix {
        let d1 = input_dir
            .unwrap_or(&PathBuf::from("./"))
            .canonicalize()
            .expect("Invalid input dir path");
        let d2 = cfg
            .output_dir()
            .unwrap_or(&PathBuf::from("./"))
            .canonicalize()
            .expect("Invalid output dir path");
        if d1 == d2 {
            return Err(anyhow!(
                "Both Input and output directories and input and output prefixes are the same"
            ));
        }
    }

    Ok(cfg)
}
