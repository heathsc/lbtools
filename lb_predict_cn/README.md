## normalize_wgs_cov
Generate binned coverage from WGS experiments (SAM/BAM/CRAM files) normalized by GC content to give estimates of copy number 

 - [Introduction](#intro)
 - [Installation](#install)
 - [General usage](#usage)
   - [Command line options](#cli)
 - [Overview of operation](#overview)
 - [Changes](#changes)

## <a name="intro"></a>Introduction

Normalize_wgs_cov is a utility that takes SAM/BAM/CRAM files from WGS or WGBS experiments and generates binned coverage (with use specified bin sizes) that are
normalized for gc content within each bin.  The estimates are normalized to give an estimate of copy number, 
so that a 'normal' chromosome will normally have an estimated coverage of 2.

## <a name="install"></a>Installation

To compile you will need an up-to-date copy of rust.  This can be
installed locally following the instructions [here](https://www.rust-lang.org/learn/get-started).  
Note that if you have rust already installed you should update it
using ``rustup update`` before trying to compile baldur.

You will also need to have [htslib](https://github.com/samtools/htslib) installed in a place 
where the rust compiler can find it.  Note that *baldur* has been tested with htslib version 1.15, but older versions may work.

Clone the repository and then from the normalize_wgs_cov directory
use cargo to compile the application:
```
    git clone https://github.com/heathsc/normalize_wgs_cov.git
    cd normalize_wgs_cov
    cargo build --release
```
If you have an error during linkage about not being able to find libhts then you will need to specify the installation location of libhts
during the build process, for example:

    RUSTFLAGS="-L /opt/share/htslib/1.15/lib/"  cargo build --release

After a successful build the executable will be found in target/release/.  It
should be copied somewhere where it can be found by the shell.

Once installed, basic help can be found by invoking normalize_wgs_cov with
the -h flag.

## <a name="usage"></a>General usage

Normalize_was_cov requires three input files

 - The first has a list of the samples to process and the path to the
SAM/BAM/INPUT files.  There is one line per samples and at least 2 tab separated columns.  The sample name is expected
to be in column one and the file path in column 2.

 - The second file has a list of the contigs to be processed.  There is one file per contig, and 1 or 2 tab separated columns.
The first column has the contig name and the second column if present has a boolean value true / false which indicates
whether the contig should be used for the GC normalization.  If the second column is missing, yes is assumed.  In normal use
the autosomes should have 'yes' in the second column while the sex chromosomes and mitochondrian should have 'no'.
 - The third file is the reference FASTA file.  This should be indexed using samtools faidx.

All input files can be compressed.  If the reference file is compressed it should be compressed using bgzip; the other files
can be compressed with many programs such as bgzip, gzip, xz, zstd or bgzip2.

In addition to the 3 input files, there are many options to
allow specification of thresholds for base and mapping qualities, limits on acceptable template lengths etc.
A minimal command line invocation would therefore be:
```
normalize_wgs_cov sample_list.txt contig_list.txt reference.fasta.gz
```
And a more typical command line showing multiple options would be:
```
normalize_wgs_cov -q 20 -Q 10 sample_list.txt contig_list.txt reference.fasta.gz
```

The second line would filter the input data to only consider reads with MAPQ scores >=10 and where base qualities are >= 20.

### <a name="cli"></a>Command line options

normalize_wgs_cov has many command line options for controlling the operation of the calling process.

| Short | Long                  | Description                                           | Default           |
|-------|-----------------------|-------------------------------------------------------|-------------------|
| b     | block-size            | Size of blocks (bins)                                 | 1000              |
| Q     | mapq                  | MAPQ threshold                                        | 0                 |
| q     | qual                  | Minimum base quality                                  | 0                 |
| M     | min-template-len      | Set minimum template length                           | 0                 |
| m     | max-template-len      | Set maximum template length                           | 0                 |
| k     | keep-duplicates       | Do not remove duplicate reads                         |                   |
| D     | ignore-duplicate-flag | Ignore duplicate flag in input files                  |                   |
|||||
| p     | prefix                | Prefix for output files                               | cov               |
| d     | dir                   | Output directory                                      | current directory |
| t     | threads               | Number of calculation threads                         | available cores   |
| @     | hts-threads           | Number of threads for SAM/BAM/CRAM files              | available cores   |
| R     | readers               | Number of file readers                                | (threads + 3) / 4 |
| l     | loglevel              | Set log level (none, error, warn, info, debug, trace) | info              |


## <a name="overview"></a>Overview of workflow

- Read in reference file and calculate GC content of genomic bins
- Read in raw coverage data per sample
- For each sample, calculate the median coverage per GC content of bin (splitting GC content level into 128 equal bins)
- Perform locally weighted regression (LOESS) to generate smoothed estimates of coverage as a function of GC content
- Normalize each sample so that an average chromosome has an expected coverage level of 2
- Output raw and normalized coverage levels

## <a name="changes"></a>Changes

- 0.2.0 First public release
- 0.1.0 First commit

