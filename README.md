# pycoMethFlow
**pycoMethFlow** is a pipeline built using [Nextflow](https://www.nextflow.io) for running minimap2 + nanopolish --call-methylation + pycoMeth pipeline on multiple samples and across multiple infrastructure in a streamlined, portable and reproducible manner.

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
-profile docker|singularity \
--samples = "/path/to/samples.txt" \
--results_dir = "/path/to/results_dir" \
--reference = "/path/to/file.fasta" \
--gtf = "/path/to/file.gff"

```

## Citation

For further information, please refer to the following manuscripts:

Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018 Sep 15;34(18):3094-3100. doi: 10.1093/bioinformatics/bty191. PMID: 29750242; PMCID: PMC6137996.

Loman, N., Quick, J. & Simpson, J. A complete bacterial genome assembled de novo using only nanopore sequencing data. Nat Methods 12, 733–735 (2015). https://doi.org/10.1038/nmeth.3444

PycoMeth: A toolbox for differential methylation testing from Nanopore methylation calls. Rene Snajder, Oliver Stegle, Marc Jan Bonder. bioRxiv 2022.02.16.480699; doi: https://doi.org/10.1101/2022.02.16.480699

[minimap2](https://github.com/lh3/minimap2)

[nanopolish](https://github.com/jts/nanopolish)

[pycoMeth](https://github.com/snajder-r/pycoMeth)
