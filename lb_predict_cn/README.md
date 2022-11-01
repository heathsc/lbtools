## lb_predict_cn
Predicts copy number in genomic bins from WGS experiments (SAM/BAM/CRAM files).

 - [Introduction](#intro)
 - [Installation](#install)
 - [General usage](#usage)
   - [Input files](#input_files)
   - [Output files](#output_files)
   - [Command line options](#cli)
 - [Overview of operation](#overview)
 - [Changes](#changes)

## <a name="intro"></a>Introduction

lb_predict_cn is a utility that takes SAM/BAM/CRAM files from WGS or WGBS experiments and
generate predictions of copy number (CN) that are normalized for GC content for genomic bins of user defined size.  
The lb_predict_cn utility can process multiple samples in parallel,
however each sample is normalized independently, only using information from the sample in question.

## <a name="install"></a>Installation

To compile you will need an up-to-date copy of rust.  This can be
installed locally following the instructions [here](https://www.rust-lang.org/learn/get-started).  
Note that if you have rust already installed you should update it
using ``rustup update`` before trying to compile baldur.

You will also need to have [htslib](https://github.com/samtools/htslib) installed in a place 
where the rust compiler can find it.  Note that *baldur* has been tested with htslib version 1.15, but older versions may work.

Clone the repository and then from the lb_predict_cn directory
use cargo to compile the application:
```
git clone https://github.com/heathsc/lb_predict_cn.git
cd lb_predict_cn
cargo build --release
```
If you have an error during linkage about not being able to find libhts then you will need to specify the installation location of libhts
during the build process, for example:

    RUSTFLAGS="-L /opt/share/htslib/1.15/lib/"  cargo build --release

After a successful build the executable will be found in target/release/.  It
should be copied somewhere where it can be found by the shell.

Once installed, basic help can be found by invoking lb_predict_cn with
the -h flag.

## <a name="usage"></a>General usage

### <a name="input_files"></a>Input files

lb_predict_cn requires three input files

 - The first has a list of the samples to process and the path to the
SAM/BAM/INPUT files.  There is one line per samples and at least 2 tab separated columns.  The sample name is expected
to be in column one and the file path in column 2.  An example of the sample file is shown below:
```
sample1  alignment/sample1.cram
sample2  alignment/sample2.cram
sample3  alignment/sample3.cram
sample4  alignment/sample4.cram
```

 - The second file has a list of the contigs to be processed and which contigs are to be used for GC normalization.
This intended use of this file is to list the chromosomes (i.e., not including unmapped contigs or alternate assemblies) 
and indicate the autosomes, as we do not want to use contigs that do not have expected copy number 2 to generate the GC normalization model.
The file has one line per contig, and 1 or 2 tab separated columns.
The first column has the contig name and the second column if present has a boolean value true / false which indicates
whether the contig should be used for the GC normalization.  If the second column is missing, yes is assumed.  In normal use
the autosomes should have 'yes' in the second column while the sex chromosomes and mitochondria should have 'no'.  An example of the contig file is shown below:
```
chr1    true
chr2    true
chr3    true
chr4    true
chr5    true
chr6    true
chr7    true
chr8    true
chr9    true
chr10   true
chr11   true
chr12   true
chr13   true
chr14   true
chr15   true
chr16   true
chr17   true
chr18   true
chr19   true
chr20   true
chr21   true
chr22   true
chrX    no
chrY    no
chrM    no
```
 - The third file is the reference FASTA file, which is used to generate the GC content per genomic bin. This should be indexed using samtools faidx.

All input files can be compressed.  If the reference file is compressed it should be compressed using bgzip; the other files
can be compressed with many programs such as bgzip, gzip, xz, zstd or bgzip2.

### <a name="output_files"></a>Output files

An output directory is made (if not already existing) for each sample, and within this directory an output file
is created per contig listed in the contig file.   By default, the sample output directories are created in the current 
directory, and each contig specific file will be names cov_*contig name*.txt i.e., cov_chr2.txt.
The behvaiour can be changed via the [Command line options](#cli), in particular look at the **dir** and **prefix** options.

The individual output files have a simple structure being tab delimited text files with 4 columns.
The 4 columns are:
 - contig name
 - mid-point of genomic bin
 - copy num estimate
 - average raw coverage within bin

A fragment of an example output file is shown below:
```
chr5    15715000        2.2442  31.4715
chr5    15725000        1.7038  18.4308
chr5    15735000        2.4320  34.1051
chr5    15745000        2.1085  32.5195
chr5    15755000        2.1093  25.1855
chr5    15765000        1.8203  25.5263
chr5    15775000        1.9811  25.7750
chr5    15785000        2.6654  42.2662
chr5    15795000        2.3168  34.6890
chr5    15805000        2.1840  32.6999
chr5    15815000        1.2509  15.6137
chr5    15825000        2.3367  32.7689
chr5    15835000        2.0017  28.0708
chr5    15845000        1.7501  23.6693
chr5    15855000        2.0719  25.8605
chr5    15865000        2.0551  30.7699
chr5    15875000        1.9828  27.8058
chr5    15885000        1.7974  26.0738
chr5    15895000        1.8685  28.8184
chr5    15905000        1.2978  16.1988
```
### <a name="cli"></a>Command line options

In addition to the 3 input files, there are many options to
allow specification of thresholds for base and mapping qualities, limits on acceptable template lengths etc.
A minimal command line invocation would therefore be:
```
lb_predict_cn sample_list.txt contig_list.txt reference.fasta.gz
```
And a more typical command line showing multiple options would be:
```
lb_predict_cn -q 20 -Q 10 sample_list.txt contig_list.txt reference.fasta.gz
```

The second line would filter the input data to only consider reads with MAPQ scores >=10 and where base qualities are >= 20.

lb_predict_cn has many command line options for controlling the operation of the calling process.

| Short | Long                  | Description                                           | Default           |
|-------|-----------------------|-------------------------------------------------------|-------------------|
| b     | block-size            | Size of blocks (bins)                                 | 10000             |
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

lb_predict_cn is multithreaded and can use multiple CPU cores for processing in parallel multiple samples or 
multiple contigs for a single sample, and can read in data for some samples while normalizing and generating the
output for other samples.  The multithreading behaviour is controlled by the **threads**, **hts-threads** 
and **readers** options.  The **threads** option sets the overall number of threads used for data processing; these
threads can be used for reading, normalizing or generating output as required.  By default, the number of threads
is set to the number of available threads in the system.  The **hts-threads** option set te number of threads used by 
htslib for reading SAM/BAM/CRAM files; depending on the speed of the storage system, htslib can use multiple
threads for reading a single contig from a single file, which can increase throughput.  Lastly the **readers** option is
used to specify the number of simultaneous read tasks.  By default this is set based on the number of threads; the optimum
setting of this will depend on the number of available cores and the performance of the storage system.

## <a name="overview"></a>Overview of workflow

- Read in reference file and calculate GC content of genomic bins
- Read in raw coverage data per sample
- For each sample, calculate the median coverage per GC content of bin (splitting GC content level into 128 equal bins)
- Perform locally weighted regression (LOESS) to generate smoothed estimates of coverage as a function of GC content
- Normalize each sample so that an average chromosome has an expected coverage level of 2
- Output estimated copy number and raw coverage per bin

## <a name="changes"></a>Changes

- 0.3.0 Rename as lb_predict_cn.  Reorganize into lbtools collection.
- 0.2.0 First public release
- 0.1.0 First commit

