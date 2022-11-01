## lbtools
A set of tools for estimating copy number from whole genome liquid biopsy sequencing experiments.

 - [Introduction](#intro)
 - [Installation](#install)
 - [Changes](#changes)

## <a name="intro"></a>Introduction

lbtools is a set of utilities for processing copy number estimates from whole genome liquid biopsy experiments.
There are currently three tools in the set:
 - [lb_predict_cn](https://github.com/heathsc/lbtools/tree/main/lb_predict_cn) which is used to generate single sample predictions of copy numnber from SAM/BAM/CRAM alignment files
  - [lb_norm](https://github.com/heathsc/lbtools/tree/main/lb_norm) which can normalize the estimates produced from 
lb_predict_cn using a set of control samples
 - [lb_region_test](https://github.com/heathsc/lbtools/tree/main/lb_region_test) is used to perform statistical tests of copy numbers
levels between samples and controls at specific genomic regions.

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
git clone https://github.com/heathsc/lb_tools.git
cd lb_tools
cargo build --release
```
If you have an error during linkage about not being able to find libhts then you will need to specify the installation location of libhts
during the build process, for example:

RUSTFLAGS="-L /opt/share/htslib/1.15/lib/"  cargo build --release

After a successful build the executable will be found in target/release/.  It
should be copied somewhere where it can be found by the shell.

Once installed, basic help can be found by invoking lb_predict_cn with
the -h flag.

## <a name="changes"></a>Changes

- 0.3.0 First public release

