## Aloft ##
Aloft is a [Bash](https://www.gnu.org/software/bash/) script script for running the excellent [Arriba](https://github.com/suhrig/arriba) RNA-Seq fusion detector in parallel, it takes charge of:

1. Downloading and compiling an [Arriba release](https://github.com/suhrig/arriba/releases) (if required).
2. Downloading references and annotation (via Arriba's [`download_references.sh`](https://github.com/suhrig/arriba/blob/master/download_references.sh))
3. Producing STAR indexes if need be.
4. Creating the known fusion events file from the COSMIC [Complete Fusion Export] (https://cancer.sanger.ac.uk/cosmic/download), and fixing none [HGNC](https://www.genenames.org/) compliant gene symbols.
5. Running the STAR aligner in parallel on a set of samples via [GNU parallel](https://www.gnu.org/software/parallel/).
6. Running Arriba in parallel - as above.
7. Running [SAMtools](http://www.htslib.org/) in parallel - enables viewing of STAR derived BAM with IGV.
8. Running Arriba's outstanding plotting script [`draw_fusions.R`](https://github.com/suhrig/arriba/blob/master/draw_fusions.R) in parallel over all samples as with previous stages.

This script was inspired by the [demo script](https://arriba.readthedocs.io/en/v1.0.1/workflow/#demo-script) `run_arriba.sh` supplied with Arriba. Aloft implements the recommended [Arriba workflow](https://arriba.readthedocs.io/en/v1.0.1/workflow/). The only difference being STAR alignment is output to disk rather than piped into Arriba, so that it can be subsequently sorted indexed and saved for manual inspection of fusions in say [IGV](http://software.broadinstitute.org/software/igv/).

### Requirements ###

* Linux or *nix like OS, with working make, [Wget](https://www.gnu.org/software/wget/), [Bash](https://www.gnu.org/software/bash/), [GNU sed](https://www.gnu.org/software/sed/), [GNU gawk](https://www.gnu.org/software/gawk/), [GNU grep](https://www.gnu.org/software/grep/), [gzip](https://www.gnu.org/software/gzip/) and [Perl](https://www.perl.org/). - Tested on [Ubuntu](https://www.ubuntu.com/) 18.04 LTS
* [STAR aligner](https://github.com/alexdobin/STAR)
* [SAMtools](http://www.htslib.org/)
* [GNU Parallel](https://www.gnu.org/software/parallel/)
* The COSMIC [Complete Fusion Export](https://cancer.sanger.ac.uk/cosmic/download) file.
* Some RNA-Seq [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files to analyse.

### Central config file ###
`aloftConfig.sh` contains various settings which will be used for execution along with comments, this Bash script is sourced by the main `aloft.sh` so it will inherit variables defined here.  Please review this before launching an analysis run.

### Sample sheet file ###
Samples are defined in a tab delimited flat file taking the form of:

```
sample_1	sample_1_R1.fastq.gz	sample_1_R2.fastq.gz
sample_2	sample_2_R1.fastq.gz	sample_2_R2.fastq.gz
sample_3	sample_3_R1.fastq.gz	sample_3_R2.fastq.gz
```

Here the first column defines the sample ID.  If more than one pair of FASTQ files exists for each sample simply add these on as extra tab delimited columns.  The path to these files should not be present in the sample sheet just their names.  As the path can be given below.

Having made such a file you, and reviewed the settings in `aloftConfig.sh` you can run the pipeline like so:

```
aloft.sh <tab delimited sample sheet> <input FASTQ path> <output dir for run> 
```

Out of the box aloft is configured to use 16 cores and will consume about 64GB of RAM during execution, the core count of various stages and number of concurrent jobs RAM allocated etc can be adjusted in `aloftConfig.sh`.