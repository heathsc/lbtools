## lb_region_test
Perform statistical tests for difference in copy number at predfined regions, and generates estimates of ctDNA 
fraction based on the observations in these regions.

 - [Introduction](#intro)
 - [Installation](#install)
 - [General usage](#usage)
   - [Input files](#input_files)
   - [Output files](#output_files)
   - [Command line options](#cli)
 - [Changes](#changes)

## <a name="intro"></a>Introduction

lb_region_test is a utility that takes the CN predictions generated by [lb_predict_cn](https://github.com/heathsc/lbtools/tree/main/lb_predict_cn)
from a set of samples and an optional set of controls, and tests for difference in copy number 
at a predefined set of genomic regions.  Tests are carried out individually on each test sample compared to the control group.  Test regions
can be defined as having a known copy number change (i.e., single copy deletion) or not.  In the case of a known
or assumed copy number change, an estimate of the fraction of ctDNA is made based on the magnitude of the copy number change.

## <a name="install"></a>Installation
For installation see the instructions for the [lbtools package](https://github.com/heathsc/lbtools)

After a successful build the executable will be found in target/release/.  It
should be copied somewhere where it can be found by the shell.

Once installed, basic help can be found by invoking lb_region_test with
the -h flag.

## <a name="usage"></a>General usage

### <a name="input_files"></a>Input files

lb_region_test requires an input file with the list of samples to process with an indication of whether
each sample is a test or control sample.  Each line has 2 tab separated column with the first giving the
sample name and the second should have either 'Test' or 'Control', for example:
```
sample1 Test
sample2 Test
sample3 Control
sample4 Control
```
Note that case is ignored for the second column, and any prefix of Test or Control can be used, so 'T', 
'tEsT', 'co', 'contRoL' are all valid.

In addition to the sample file, lb_region_test also required a file giving the set of regions to be tested.  The
region file must be specified using the **region-list** option.  The format of the region file is a tab-separated
text file with one line per region.  Each line should have 3 or 4 columns:
 - Column 1 has a text description of the region
 - Column 2 has the chromosome name
 - Column 3 has one or more chromosome ranges separated by commas
 - Column 4 (if present) shows the expected copy number change for this region (-1 = single copy deletion, 1 = single copy duplication etc.)

An example region file is shown below:
```
1p36.32/1p36.11 chr1    3652516-3736201,26696015-26782104       -1
MYCN            chr2    15940550-15947004
11q22.3-q24.2   chr11   108223067-125676255                     -1
17q21.2-q21.33  chr17   42287547-51162168                        2

```

Files with the copy number predictions from 
[lb_predict_cn](https://github.com/heathsc/lbtools/tree/main/lb_predict_cn)
will be searched for according to the **input-dir** and **input-prefix** [command line options](#cli).
All input files can be compressed with many programs such as bgzip, gzip, xz, zstd or bgzip2.

### <a name="output_files"></a>Output file

A report file with the results of the test will be generated and output to <stdout>.  This behaviour 
can be changed with the **output** option that allows an output file to be specified.

The output file is a tab separated file with 9 columns:
 - Sample name
 - Description of the region being tested (from the region file)
 - Number of control samples with data on the test region
 - Standard deviation of control sample normalized coverage in the test region
 - Estimated copy number of sample in the test region
 - Estimated ctDNA fraction with 95% confidence limits for regions with a given expected copy number change
 - p value from t test (with n-1 df) of difference in copy number between the sample and the control group
 - Corrected p-value controlling for false discovery rate (Benjamini & Hchberg 1995)

### <a name="cli"></a>Command line options

In addition to the sample input file and region file, there are options to set the input directories
and prefixes, and the output file.

A minimal command line invocation would be:
```
lb_region_test -r region_list.txt sample_list.txt
```
And a more typical command line showing multiple options would be:
```
lb_predict_cn -P cov -D normalized -r region_list.txt -o report.tsv sample_list.txt
```

| Short | Long         | Description                                           | Default           |
|-------|--------------|-------------------------------------------------------|-------------------|
| P     | input-prefix | Prefix for input files                                | cov               |
| D     | input-dir    | Input directory                                       | current directory |
| o     | output       | Output file                                           | <stdout>          |
| r     | region-list  | File with list of regions to test                     |                   |
| l     | loglevel     | Set log level (none, error, warn, info, debug, trace) | info              |

## <a name="changes"></a>Changes

- 0.3.0 First public release
