mod cli;
mod config;
mod io;
mod process;
mod sample;

#[macro_use]
extern crate log;
#[macro_use]
extern crate anyhow;

use anyhow::Context;

fn main() -> anyhow::Result<()> {
    let cfg = cli::handle_cli().with_context(|| "Error processing command line arguments")?;
    process::process_samples(&cfg)
}
