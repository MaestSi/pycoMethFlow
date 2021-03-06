# pycoMethFlow
**pycoMethFlow** is a [Nextflow](https://www.nextflow.io) pipeline for running [minimap2](https://github.com/lh3/minimap2) + [nanopolish](https://github.com/jts/nanopolish) call-methylation + [pycoMeth](https://github.com/snajder-r/pycoMeth) on multiple samples and across multiple infrastructures in a streamlined, portable and reproducible manner.

## Getting started

**Prerequisites**

* [Nextflow](https://nf-co.re/usage/installation)
* [Docker](https://docs.docker.com/engine/install/) or [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html)                                                                                  
                                                                                   
**Installation**

```
git clone https://github.com/MaestSi/pycoMethFlow.git
cd pycoMethFlow
chmod 755 *
```

## Usage

The pycoMethFlow pipeline requires you to open pycomethflow.conf configuration file and set the desired options. Then, you can run the pipeline using either docker or singularity environments just specifying a value for the -profile variable.

```
nextflow -c pycomethflow.conf run pycomethflow.nf \
--samples = "samples.txt" \
--results_dir = "results_dir" \
--reference = "file.fasta" \
--gff = "file.gff" \
-profile docker

    Mandatory argument:
    -profile                        Configuration profile to use. Available: docker, singularity
    Other mandatory arguments which may be specified in the pycomethflow.conf file
      --samples                     Path to tab separated sample sheet containing sample_name /path/to/file.fastq /path/to/fast5_dir /path/to/sequencing_summary.txt
      --results_dir                 Directory where results are stored
      --reference                   Reference file in fasta format
      --gff                         Annotation file in gff3 format
```

## Citation

For further information, please refer to the following manuscripts:

Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018 Sep 15;34(18):3094-3100. doi: 10.1093/bioinformatics/bty191. PMID: 29750242; PMCID: PMC6137996.

Simpson, J., Workman, R., Zuzarte, P. et al. Detecting DNA cytosine methylation using nanopore sequencing. Nat Methods 14, 407???410 (2017). https://doi.org/10.1038/nmeth.4184

PycoMeth: A toolbox for differential methylation testing from Nanopore methylation calls. Rene Snajder, Oliver Stegle, Marc Jan Bonder. bioRxiv 2022.02.16.480699; doi: https://doi.org/10.1101/2022.02.16.480699

[minimap2](https://github.com/lh3/minimap2)

[nanopolish](https://github.com/jts/nanopolish)

[pycoMeth](https://github.com/a-slide/pycoMeth)

[pycoMeth documentation](https://a-slide.github.io/pycoMeth/)
