use std::{
    fmt,
    str::FromStr,
    io::BufRead,
};

use clap::ArgMatches;

/// LogLevel
///
/// Represents minimum level of messages that will be logged
///
#[derive(Debug, Clone, Copy)]
pub struct LogLevel {
    pub level: usize,
}

impl FromStr for LogLevel {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "error" => Ok(LogLevel { level: 0 }),
            "warn" => Ok(LogLevel { level: 1 }),
            "info" => Ok(LogLevel { level: 2 }),
            "debug" => Ok(LogLevel { level: 3 }),
            "trace" => Ok(LogLevel { level: 4 }),
            "none" => Ok(LogLevel { level: 5 }),
            _ => Err("no match"),
        }
    }
}

impl LogLevel {
    pub fn is_none(&self) -> bool {
        self.level > 4
    }
    pub fn get_level(&self) -> usize {
        if self.level > 4 {
            0
        } else {
            self.level
        }
    }
}

impl fmt::Display for LogLevel {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let level_str = ["error", "warn", "info", "debug", "trace", "none"];
        if self.level < 6 {
            write!(f, "{}", level_str[self.level])
        } else {
            write!(f, "unknown")
        }
    }
}

/// Initialize logging from command line arguments
pub fn init_log(m: &ArgMatches) {
    let verbose = m
        .get_one::<LogLevel>("loglevel")
        .copied()
        .unwrap_or_else(|| LogLevel::from_str("info").expect("Could not set loglevel info"));
    let quiet = verbose.is_none() || m.get_flag("quiet");
    let ts = m
        .get_one::<stderrlog::Timestamp>("timestamp")
        .copied()
        .unwrap_or(stderrlog::Timestamp::Off);

    stderrlog::new()
        .quiet(quiet)
        .verbosity(verbose.get_level())
        .timestamp(ts)
        .init()
        .unwrap();
}

/// Read in next line and split on tabs after trimming white space
pub fn get_next_line<'a, R: BufRead>(
    rdr: &mut R,
    buf: &'a mut String,
) -> anyhow::Result<Option<Vec<&'a str>>> {
    buf.clear();
    if rdr.read_line(buf)? == 0 {
        Ok(None)
    } else {
        Ok(Some(buf.trim().split('\t').collect()))
    }
}