#!/usr/bin/env nextflow
/*
========================================================================================
                         MaestSi/pycoMethFlow
========================================================================================
 MaestSi/pycoMethFlow analysis pipeline.
 #### Homepage / Documentation
 https://github.com/MaestSi/pycoMethFlow
----------------------------------------------------------------------------------------
*/
def helpMessage() {
        log.info"""
    Usage:
    nextflow -c pycomethflow.conf run pycomethflow.nf --samples = "samples.txt" --results_dir = "results_dir" --reference = "file.fasta" --gtf = "file.gff" 
-profile docker
    Mandatory argument:
    -profile                        Configuration profile to use. Available: docker, singularity
    Other mandatory arguments which may be specified in the pycomethflow.conf file
      --samples                     Path to tab separated sample sheet containing sample_name /path/to/file.fastq /path/to/fast5_dir /path/to/sequencing_summary.txt
      --results_dir                 Directory where results are stored
      --reference                   Reference file in fasta format
      --gtf                         Annotation file in gtf format
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Input of sample names, conditions, and FAST5s path.
Channel
    .fromPath(params.samples)
    .splitCsv(header: true, sep:'\t')
    .map{ row-> tuple(row.sample, file(row.fastq)) }
    .set{samples_alignment}

// Alignment
process alignment {
    input:
        tuple val(sample), file(fastq) from samples_alignment

    output:
        val(sample) into alignment_nanopolish

    script:
    if(params.alignment)
    """
       	mkdir -p ${params.results_dir}/${sample}/alignment/
        
	minimap2 -ax map-ont -t ${task.cpus} ${params.reference} ${fastq} | samtools view -hSb | samtools sort -@ ${task.cpus} -o ${params.results_dir}/${sample}/alignment/minimap.bam
        samtools index -@ ${task.cpus} ${params.results_dir}/${sample}/alignment/minimap.bam
           
        ln -s ${params.results_dir}/${sample}/alignment/minimap.bam ./minimap.bam
        ln -s ${params.results_dir}/${sample}/alignment/minimap.bam.bai ./minimap.bam.bai
    """
    else
    """
        ln -s ${params.results_dir}/${sample}/alignment/minimap.bam ./minimap.bam
       	ln -s ${params.results_dir}/${sample}/alignment/minimap.bam.bai ./minimap.bam.bai
    """
}

// Nanopolish
process nanopolish {
    input:
    val(sample) from alignment_nanopolish

    output:
	tuple val(sample), file('nanopolish_cpg_methylation.tsv') into nanopolish_pycomethCpGAggregate

    script:
    if(params.nanopolish)
    """
       	mkdir -p ${params.results_dir}/${sample}/nanopolish/
        
	fastq=\$(grep \"^\"${sample}\"\t\" ${params.samples} | cut -f2)
        fast5=\$(grep \"^\"${sample}\"\t\" ${params.samples} | cut -f3)
        sequencing_summary=\$(grep \"^\"${sample}\"\t\" ${params.samples} | cut -f4)
        
	nanopolish index -d \${fast5} \${fastq} -s \${sequencing_summary}
        nanopolish call-methylation --reads \${fastq} --bam ${params.results_dir}/${sample}/alignment/minimap.bam --genome ${params.reference} --methylation cpg --threads ${task.cpus} > ${params.results_dir}/${sample}/nanopolish/nanopolish_cpg_methylation.tsv

        ln -s ${params.results_dir}/${sample}/nanopolish/nanopolish_cpg_methylation.tsv ./nanopolish_cpg_methylation.tsv
    """
    else
    """
        ln -s ${params.results_dir}/${sample}/nanopolish/nanopolish_cpg_methylation.tsv ./nanopolish_cpg_methylation.tsv
    """
}

// pycoMeth CGI Finder
process pycomethCGIFinder {
    input:

    output:
    file('CGI_Finder.bed') into pycomethCGIFinder_pycomethIntervalAggregate

    script:
    if(params.pycomethCGIFinder)
    """
       	mkdir -p ${params.results_dir}/pycometh/
        
	pycoMeth CGI_Finder -f ${params.reference} -b ${params.results_dir}/pycometh/CGI_Finder.bed -t ${params.results_dir}/pycometh/CGI_Finder.tsv -m ${params.CGI_Finder_m} -w ${params.CGI_Finder_w} -c ${params.CGI_Finder_c} -r ${params.CGI_Finder_r}

        if [[ ! -f "${params.results_dir}/pycometh/CGI_Finder.bed" ]]; then touch "${params.results_dir}/pycometh/CGI_Finder.bed"; fi
	
        ln -s ${params.results_dir}/pycometh/CGI_Finder.bed ./CGI_Finder.bed
        ln -s ${params.results_dir}/pycometh/CGI_Finder.tsv ./CGI_Finder.tsv
        
    """
    else
    """
        touch ./CGI_Finder.bed ./CGI_Finder.tsv
      
    """
}
process pycomethCpGAggregate {
    input:
    tuple val(sample), file('nanopolish_cpg_methylation') from nanopolish_pycomethCpGAggregate

    output:
    tuple val(sample), file ('CpG_Aggregate.bed') into pycomethCpGAggregate_pycomethIntervalAggregate
    file('CpG_Aggregate.tsv') into pycomethCpGAggregate_pycomethMethComp

    script:
    if(params.pycomethCpGAggregate)
    """
        mkdir -p ${params.results_dir}/${sample}/pycometh/

        pycoMeth CpG_Aggregate -i 'nanopolish_cpg_methylation' -f ${params.reference} -b ${params.results_dir}/${sample}/pycometh/CpG_Aggregate.bed -t ${params.results_dir}/${sample}/pycometh/CpG_Aggregate.tsv -s ${sample} -d ${params.CpG_Aggregate_d} -l ${params.CpG_Aggregate_l}
        
	if [[ ! -f "${params.results_dir}/${sample}/pycometh/CpG_Aggregate.bed" ]]; then touch "${params.results_dir}/${sample}/pycometh/CpG_Aggregate.bed"; fi
	
        ln -s ${params.results_dir}/${sample}/pycometh/CpG_Aggregate.bed ./CpG_Aggregate.bed
        ln -s ${params.results_dir}/${sample}/pycometh/CpG_Aggregate.tsv ./CpG_Aggregate.tsv
        
    """
    else
    """
        ln -s ${params.results_dir}/${sample}/pycometh/CpG_Aggregate.bed ./CpG_Aggregate.bed*
        ln -s ${params.results_dir}/${sample}/pycometh/CpG_Aggregate.tsv ./CpG_Aggregate.tsv*
  
    """
}

process pycomethIntervalAggregate {
    input:
    file('CGI_Finder.bed') from pycomethCGIFinder_pycomethIntervalAggregate
    tuple val(sample), file ('CpG_Aggregate.bed') from pycomethCpGAggregate_pycomethIntervalAggregate

    output:
    file('Interval_Aggregate.tsv') into pycomethIntervalAggregate_pycomethMethComp

    script:
    if(params.pycomethIntervalAggregate)
    """
        mkdir -p ${params.results_dir}/${sample}/pycometh/;
        if [[ ${params.pycomethCGIFinder} ]]; then
            pycoMeth Interval_Aggregate -i ${params.results_dir}/${sample}/pycometh/CpG_Aggregate.tsv -f ${params.reference} -b ${params.results_dir}/${sample}/pycometh/Interval_Aggregate.bed  -t ${params.results_dir}/${sample}/pycometh/Interval_Aggregate.tsv -n ${params.Interval_Aggregate_n} -m ${params.Interval_Aggregate_m} -s ${sample} -l ${params.Interval_Aggregate_l};
        else
            pycoMeth Interval_Aggregate -i ${params.results_dir}/${sample}/pycometh/CpG_Aggregate.tsv -f ${params.reference} -a ${params.results_dir}/${sample}/pycometh/CGI_Aggregate.bed  -b ${params.results_dir}/${sample}/pycometh/Interval_Aggregate.bed  -t ${params.results_dir}/${sample}/pycometh/Interval_Aggregate.tsv -n ${params.Interval_Aggregate_n}  -m ${params.Interval_Aggregate_m} -s ${sample} -l ${params.Interval_Aggregate_l} 
        fi
        
	if [[ ! -f "${params.results_dir}/${sample}/pycometh/Interval_Aggregate.tsv" ]]; then touch "${params.results_dir}/${sample}/pycometh/Interval_Aggregate.tsv"; fi
	
        ln -s ${params.results_dir}/${sample}/pycometh/Interval_Aggregate.bed ./Interval_Aggregate.bed;
        ln -s ${params.results_dir}/${sample}/pycometh/Interval_Aggregate.tsv ./Interval_Aggregate.tsv;
   
    """
    else
    """
        ln -s ${params.results_dir}/${sample}/pycometh/Interval_Aggregate.bed ./Interval_Aggregate.bed
        ln -s ${params.results_dir}/${sample}/pycometh/Interval_Aggregate.tsv ./Interval_Aggregate.tsv
    """
}

process pycomethMethComp {
    input:
    file('Interval_Aggregate.tsv*') from pycomethIntervalAggregate_pycomethMethComp.collect()
    file('CpG_Aggregate.tsv') from pycomethCpGAggregate_pycomethMethComp.collect()

    output:
    file ('Meth_Comp_Interval.tsv') into pycomethMethComp_pycomethCompReport

    script:
    if(params.pycomethMethComp)
    """
        mkdir -p ${params.results_dir}/pycometh/
        
	pycoMeth Meth_Comp -i Interval_Aggregate.tsv* -f ${params.reference} -b ${params.results_dir}/pycometh/Meth_Comp_Interval.bed -t ${params.results_dir}/pycometh/Meth_Comp_Interval.tsv -m ${params.Meth_Comp_m} -l ${params.Meth_Comp_l} --pvalue_adj_method ${params.Meth_Comp_pvalue_adj_method} --pvalue_threshold ${params.Meth_Comp_pvalue_threshold} ${params.Meth_Comp_only_tested_sites}
	pycoMeth Meth_Comp -i CpG_Aggregate.tsv* -f ${params.reference} -b ${params.results_dir}/pycometh/Meth_Comp_CpG.bed -t ${params.results_dir}/pycometh/Meth_Comp_CpG.tsv -m ${params.Meth_Comp_m} -l ${params.Meth_Comp_l} --pvalue_adj_method ${params.Meth_Comp_pvalue_adj_method} --pvalue_threshold ${params.Meth_Comp_pvalue_threshold} ${params.Meth_Comp_only_tested_sites}
	
	if [[ ! -f "${params.results_dir}/pycometh/Meth_Comp_Interval.tsv" ]]; then touch "${params.results_dir}/pycometh/Meth_Comp_Interval.tsv"; fi
        
	ln -s ${params.results_dir}/pycometh/Meth_Comp_Interval.bed ./Meth_Comp_Interval.bed
        ln -s ${params.results_dir}/pycometh/Meth_Comp_Interval.tsv ./Meth_Comp_Interval.tsv
        ln -s ${params.results_dir}/pycometh/Meth_Comp_CpG.bed ./Meth_Comp_CpG.bed
        ln -s ${params.results_dir}/pycometh/Meth_Comp_CpG.tsv ./Meth_Comp_CpG.tsv
    """
    else
    """
       ln -s ${params.results_dir}/pycometh/Meth_Comp_Interval.bed ./Meth_Comp_Interval.bed
       ln -s ${params.results_dir}/pycometh/Meth_Comp_Interval.tsv ./Meth_Comp_Interval.tsv
       ln -s ${params.results_dir}/pycometh/Meth_Comp_CpG.bed ./Meth_Comp_CpG.bed
       ln -s ${params.results_dir}/pycometh/Meth_Comp_CpG.tsv ./Meth_Comp_CpG.tsv
    """
}
process pycomethCompReport {
    input:
    file ('Meth_Comp_Interval.tsv') from pycomethMethComp_pycomethCompReport

    output:

    script:
    if(params.pycomethCompReport)
    """
        mkdir -p ${params.results_dir}/pycometh/
        
        pycoMeth Comp_Report -i 'Meth_Comp_Interval.tsv' -g ${params.gtf} -f ${params.reference} -o ${params.results_dir}/pycometh -n ${params.Comp_Report_n} -d ${params.Comp_Report_d} --pvalue_threshold ${params.Comp_Report_pvalue_threshold} --min_diff_llr ${params.Comp_Report_min_diff_llr} --n_len_bin ${params.Comp_Report_n_len_bin} ${params.Comp_Report_export_static_plots} ${params.Comp_Report_report_non_significant};
        ln -s ${params.results_dir}/pycometh/pycoMeth_summary_report.html ./pycoMeth_summary_report.html;
    """
    else
    """
        echo "Skipped"
    """
}
