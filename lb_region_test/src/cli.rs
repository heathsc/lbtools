use std::path::PathBuf;

use clap::{
    crate_authors, crate_description, crate_name, crate_version, value_parser, Arg, ArgAction,
    Command,
};

use anyhow::Context;

use utils::{init_log, LogLevel};

use crate::{config::*, region::*, sample::*};

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
            Arg::new("input_prefix")
                .short('P')
                .long("input-prefix")
                .value_parser(value_parser!(String))
                .value_name("STRING")
                .default_value("cov")
                .help("Set prefix for input file names"),
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
            Arg::new("output")
                .short('o')
                .long("output-file")
                .value_parser(value_parser!(PathBuf))
                .value_name("PATH")
                .help("Set output file [default: <stdout>]"),
        )
        .arg(
            Arg::new("region_list")
                .short('r')
                .long("region-list")
                .value_parser(value_parser!(PathBuf))
                .value_name("PATH")
                .required(true)
                .help("FIle with list of regions to be tested"),
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

    let (regions, mut ctg_hash) = read_region_file(
        m.get_one::<PathBuf>("region_list")
            .expect("Missing region list file"),
    )
    .with_context(|| "Could not read from region list file")?;

    let input_prefix = m
        .get_one::<String>("input_prefix")
        .expect("Missing default output prefix")
        .clone();

    // Read in sample list
    let mut samples = read_sample_list_from_file(
        m.get_one::<PathBuf>("sample_list")
            .expect("Missing sample list file"),
    )
    .with_context(|| "Could not read from sample list file")?;

    let input_dir = m.get_one::<PathBuf>("input_dir");

    get_input_files_and_contig_list(&mut samples, input_dir, &input_prefix, &ctg_hash)
        .with_context(|| "Error collecting input files")?;

    let contigs: Vec<_> = ctg_hash.drain().collect();

    debug!("Number of contigs found: {}", contigs.len());

    let output = m.get_one::<PathBuf>("output").map(|s| s.to_owned());
    Ok(Config::new(samples, contigs, regions, output))
}
