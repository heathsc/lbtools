use std::{
    num::{NonZeroU32, NonZeroUsize},
    path::PathBuf,
};

use clap::{crate_version, value_parser, Arg, ArgAction, Command};

use crate::{
    config::Config,
    contig::contig_hash_from_file,
    gc::GcData,
    utils::{init_log, LogLevel},
};

/// Set up definition of command options for clap
fn cli_model() -> Command {
    Command::new("normalize")
        .about("Normalize binned coverage using gc content")
        .version(crate_version!())
        .author("Simon Heath")
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
                .short('q')
                .action(ArgAction::SetTrue)
                .long("quiet")
                .conflicts_with("loglevel")
                .help("Silence all output"),
        )
        .arg(
            Arg::new("block_size")
                .short('b')
                .long("block-size")
                .value_parser(value_parser!(NonZeroU32))
                .value_name("INT")
                .default_value("1000")
                .help("Set block size in base pairs"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .value_parser(value_parser!(NonZeroUsize))
                .value_name("INT")
                .help("Set number of threads [default: available cores]"),
        )
        .arg(
            Arg::new("max_template_len")
                .short('m')
                .long("max-template-len")
                .value_parser(value_parser!(usize))
                .value_name("INT")
                .help("Set maximum template length"),
        )
        .arg(
            Arg::new("prefix")
                .short('p')
                .long("prefix")
                .value_parser(value_parser!(String))
                .value_name("STRING")
                .default_value("cov")
                .help("Set maximum template length"),
        )
        .arg(
            Arg::new("dir")
                .short('d')
                .long("dir")
                .value_parser(value_parser!(PathBuf))
                .value_name("PATH")
                .help("Set maximum template length [default: current directory]"),
        )
        .arg(
            Arg::new("min_template_len")
                .short('M')
                .long("min-template-len")
                .value_parser(value_parser!(usize))
                .value_name("INT")
                .default_value("0")
                .help("Set minimum template length"),
        )
        .arg(
            Arg::new("sample_file")
                .value_parser(value_parser!(PathBuf))
                .value_name("SAMPLE_FILE")
                .required(true)
                .help("Input file with list of sample names and file paths"),
        )
        .arg(
            Arg::new("contig_file")
                .value_parser(value_parser!(PathBuf))
                .value_name("CONTIG_FILE")
                .required(true)
                .help("Input file with list of contig names"),
        )
        .arg(
            Arg::new("reference_file")
                .value_parser(value_parser!(PathBuf))
                .value_name("REFERENCE_FILE")
                .required(true)
                .help("Input FASTA file with reference sequence"),
        )
}

/// Handle command line options.  Set up Config structure
pub fn handle_cli() -> anyhow::Result<Config> {
    // Get matches from command line
    let m = cli_model().get_matches();

    // Setup logging
    init_log(&m);

    debug!("Processing command line options");

    // Read in contig list
    let ctg_hash = contig_hash_from_file(m.get_one::<PathBuf>("contig_file").unwrap())?;

    // Set up threads
    let nt = m
        .get_one::<NonZeroUsize>("threads")
        .map(|x| usize::from(*x))
        .unwrap_or_else(num_cpus::get);

    // Set up gc information from reference
    let block_size = u32::from(*m.get_one::<NonZeroU32>("block_size").unwrap());
    let gc_data = GcData::from_reference(
        m.get_one::<PathBuf>("reference_file").unwrap(),
        block_size,
        nt,
        &ctg_hash,
    )?;

    Ok(Config::new())
}
