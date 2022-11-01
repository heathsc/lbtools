use std::{
    num::{NonZeroU32, NonZeroUsize},
    path::PathBuf,
};

use clap::{
    crate_authors, crate_description, crate_name, crate_version, value_parser, Arg, ArgAction,
    Command,
};

use utils::{init_log, LogLevel};

use crate::{
    config::Config, contig::contig_hash_from_file, gc::GcData, sample::sample_vec_from_file,
};

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
            Arg::new("block_size")
                .short('b')
                .long("block-size")
                .value_parser(value_parser!(NonZeroU32))
                .value_name("INT")
                .default_value("10000")
                .help("Set block size in base pairs"),
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
            Arg::new("hts_threads")
                .short('@')
                .long("hts-threads")
                .value_parser(value_parser!(NonZeroUsize))
                .value_name("INT")
                .help("Set number of threads for sam/bam/cram reading [default: available cores]"),
        )
        .arg(
            Arg::new("readers")
                .short('R')
                .long("readers")
                .value_parser(value_parser!(NonZeroUsize))
                .value_name("INT")
                .help("Set maximum number of file readers operating simultaneously for sam/bam/cram reading [default: (threads + 3) / 4]"),
        )
        .arg(
            Arg::new("mapq")
                .short('Q')
                .long("mapq")
                .value_parser(value_parser!(u8))
                .value_name("INT")
                .default_value("0")
                .help("Minimum MAPQ for reads"),
        )
        .arg(
            Arg::new("qual")
                .short('q')
                .long("qual")
                .value_parser(value_parser!(u8))
                .value_name("INT")
                .default_value("0")
                .help("Minimum base quality"),
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
                .help("Set prefix for output file names"),
        )
        .arg(
            Arg::new("dir")
                .short('d')
                .long("dir")
                .value_parser(value_parser!(PathBuf))
                .value_name("PATH")
                .help("Set output directory [default: current directory]"),
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
            Arg::new("keep_duplicates")
                .short('k')
                .long("keep-duplicates")
                .action(ArgAction::SetTrue)
                .help("Do not remove duplicates"),
        )
        .arg(
            Arg::new("ignore_dup_flag")
                .short('D')
                .long("ignore-duplicate-flag")
                .action(ArgAction::SetTrue)
                .help("Ignore duplicate flag in input file"),
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
        .get_one::<NonZeroUsize>("n_tasks")
        .map(|x| usize::from(*x))
        .unwrap_or_else(num_cpus::get);

    let hts_threads = m
        .get_one::<NonZeroUsize>("hts_threads")
        .map(|x| usize::from(*x))
        .unwrap_or_else(num_cpus::get);

    let n_readers = m
        .get_one::<NonZeroUsize>("readers")
        .map(|x| usize::from(*x))
        .unwrap_or_else(|| (nt + 3) >> 2);

    let block_size = u32::from(*m.get_one::<NonZeroU32>("block_size").unwrap());

    let reference = m
        .get_one::<PathBuf>("reference_file")
        .expect("Missing reference file")
        .clone();

    // Set up gc information from reference
    let gc_data = GcData::from_reference(&reference, block_size, nt, &ctg_hash)?;

    let prefix = m
        .get_one::<String>("prefix")
        .expect("Missing default prefix")
        .clone();

    let samples = sample_vec_from_file(
        m.get_one::<PathBuf>("sample_file")
            .expect("Missing sample list file"),
    )?;

    let mut cfg = Config::new(samples, ctg_hash, gc_data, reference, prefix);

    if let Some(x) = m.get_one::<usize>("min_template_len") {
        cfg.set_min_template_len(*x)?
    }
    if let Some(x) = m.get_one::<usize>("max_template_len") {
        cfg.set_max_template_len(*x)?
    }

    if let Some(p) = m.get_one::<PathBuf>("dir") {
        cfg.set_output_dir(p)
    }

    if let Some(x) = m.get_one::<u8>("mapq") {
        cfg.set_min_mapq(*x)
    }
    if let Some(x) = m.get_one::<u8>("qual") {
        cfg.set_min_qual(*x)
    }
    if m.get_flag("keep_duplicates") {
        cfg.set_keep_duplicates()
    }
    if m.get_flag("ignore_dup_flag") {
        cfg.set_ignore_dup_flag()
    }

    cfg.set_hts_threads(hts_threads);

    cfg.set_block_size(block_size);
    cfg.set_n_tasks(nt);
    cfg.set_n_readers(n_readers);

    Ok(cfg)
}
