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
    nextflow -c pycomethflow.conf run pycomethflow.nf --samples = "samples.txt" --results_dir = "results_dir" --reference = "file.fasta" --gff = "file.gff" 
-profile docker
    Mandatory argument:
    -profile                        Configuration profile to use. Available: docker, singularity
    Other mandatory arguments which may be specified in the pycomethflow.conf file
      --samples                     Path to tab separated sample sheet containing sample_name /path/to/file.fastq /path/to/fast5_dir /path/to/sequencing_summary.txt
      --results_dir                 Directory where results are stored
      --reference                   Reference file in fasta format
      --gff                         Annotation file in gff format
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

// Indexing
process indexing {
    input:

    output:
       file('chr_names.txt') into indexing_chrNames
    script:
    if (params.indexing)
    """
        mkdir -p ${params.results_dir}/pycometh/
        if [[ ! -f "${params.reference}.fai" ]]; then samtools faidx ${params.reference}; fi
        
        cat ${params.reference}.fai | cut -f1 > chr_names.txt

    """
    else
    """
        mkdir -p ${params.results_dir}/pycometh/
        cat ${params.reference}.fai | cut -f1 > chr_names.txt
    """
}

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
	val(sample) into nanopolish_pycomethMethSeg

    script:
    if(params.nanopolish)
    """
       	mkdir -p ${params.results_dir}/${sample}/nanopolish/
        
	    fastq=\$(grep \"^\"${sample}\"\t\" ${params.samples} | cut -f2)
        fast5=\$(grep \"^\"${sample}\"\t\" ${params.samples} | cut -f3)
        sequencing_summary=\$(grep \"^\"${sample}\"\t\" ${params.samples} | cut -f4)
        
	    nanopolish index -d \${fast5} \${fastq} -s \${sequencing_summary}
        nanopolish call-methylation --reads \${fastq} --bam ${params.results_dir}/${sample}/alignment/minimap.bam --genome ${params.reference} --methylation cpg --threads ${task.cpus} > ${params.results_dir}/${sample}/nanopolish/nanopolish_cpg_methylation.tsv

        meth5 create_m5 --input_paths ${params.results_dir}/${sample}/nanopolish/nanopolish_cpg_methylation.tsv --output_file ${params.results_dir}/${sample}/nanopolish/nanopolish_cpg_methylation.m5

        ln -s ${params.results_dir}/${sample}/nanopolish/nanopolish_cpg_methylation.m5 ./nanopolish_cpg_methylation.m5
    """
    else
    """
        ln -s ${params.results_dir}/${sample}/nanopolish/nanopolish_cpg_methylation.m5 ./nanopolish_cpg_methylation.m5

    """
}

indexing_chrNames
        .splitCsv()
        .set{chrNames_pycomethMethSeg}

// pycoMeth Meth_Seg
process pycomethMethSeg {
    input:
    val(chr) from chrNames_pycomethMethSeg
    val(sample) from nanopolish_pycomethMethSeg.collect()

    output:
    val(sample) into pycomethMethSeg_pycomethMethComp

    script:
    if(params.pycomethMethSeg)
    """
        mkdir -p ${params.results_dir}/pycometh/;
        chr_val=\$(echo "${chr}" | sed \'s/\\[//\' | sed \'s/\\]//\')
        samples_list=\$(echo "${sample}" | sed \'s/\\[//\' | sed \'s/\\]//\' | sed \'s/,/\\n/g\')

        m5=""
        for s in \$samples_list; do
            m5_curr=\$(find ${params.results_dir}/\$s/nanopolish | grep \"\\.m5\$\");
            m5=\$m5\" \"\$m5_curr;
        done
        
        chr_flag=1
        for file in \$(find ${params.results_dir}/*/nanopolish/ | grep nanopolish_cpg_methylation.tsv); do
            if ! grep -q \$chr_val \$file ; then chr_flag=""; fi
        done
        
        if [ ! -z "\$chr_flag" ]; then
            pycoMeth Meth_Seg -i \$m5 -c \$chr_val -b ${params.results_dir}/pycometh/Meth_Seg_\$chr_val  -t ${params.results_dir}/pycometh/Meth_Seg_\$chr_val.tsv -p ${task.cpus} --reader_workers ${task.cpus} -m ${params.Meth_Seg_m} -w ${params.Meth_Seg_w} ${params.Meth_Seg_print_diff_met};
        fi
       
    """
    else
    """
        mkdir -p ${params.results_dir}/pycometh/
        
        pycoMeth CGI_Finder -f ${params.reference} -b ${params.results_dir}/pycometh/CGI_Finder.bed -t ${params.results_dir}/pycometh/CGI_Finder.tsv -m ${params.CGI_Finder_m} -w ${params.CGI_Finder_w} -c ${params.CGI_Finder_c} -r ${params.CGI_Finder_r}

        ln -s ${params.results_dir}/pycometh/CGI_Finder.tsv ./CGI_Finder.tsv
        ln -s ${params.results_dir}/pycometh/CGI_Finder.bed ./CGI_Finder.bed

    """
}                                     

process pycomethMethComp {
    input:
    val(sample) from pycomethMethSeg_pycomethMethComp.collect()

    output:
    file ('Meth_Comp.tsv') into pycomethMethComp_pycomethCompReport
    val(sample) into pycomethMethComp_pycomethCompReport2

    script:
    if(params.pycomethMethComp)
    """
        mkdir -p ${params.results_dir}/pycometh/

        if [[ ${params.pycomethMethSeg} ]]; then
            cat ${params.results_dir}/pycometh/Meth_Seg_*.bedGraph >> ${params.results_dir}/pycometh/Meth_Seg.bedGraph
            cat ${params.results_dir}/pycometh/Meth_Seg_*.tsv >> ${params.results_dir}/pycometh/Meth_Seg.tsv
            ln -s ${params.results_dir}/pycometh/Meth_Seg.bedGraph ./Meth_Seg.bedGraph;
            ln -s ${params.results_dir}/pycometh/Meth_Seg.tsv ./Meth_Seg.tsv;
            #for s in \$samples_list; do
            #    cat ${params.results_dir}/\$s/pycometh/Meth_Seg_.*.bedGraph >> ${params.results_dir}/\$s/pycometh/Meth_Seg.bedGraph
            #    cat ${params.results_dir}/\$s/pycometh/Meth_Seg_.*.tsv >> ${params.results_dir}/\$s/pycometh/Meth_Seg.tsv
            #done
            #ln -s ${params.results_dir}/\$s/pycometh/Meth_Seg.bedGraph ./Meth_Seg.bedGraph;
            #ln -s ${params.results_dir}/\$s/pycometh/Meth_Seg.tsv ./Meth_Seg.tsv;
        else 
            ln -s ${params.results_dir}/pycometh/CGI_Finder.tsv ./CGI_Finder.tsv
            ln -s ${params.results_dir}/pycometh/CGI_Finder.bed ./CGI_Finder.bed
        fi

        samples_list=\$(echo "${sample}" | sed \'s/\\[//\' | sed \'s/\\]//\' | sed \'s/,/\\n/g\')
        m5=""
        for s in \$samples_list; do
            m5_curr=\$(find ${params.results_dir}/\$s/nanopolish | grep \"\\.m5\$\");
            m5=\$m5\" \"\$m5_curr;
        done
        
        if [[ ${params.pycomethMethSeg} ]]; then
            pycoMeth Meth_Comp -i \$m5 -a ${params.results_dir}/pycometh/Meth_Seg.bedGraph -s \$samples_list -f ${params.reference} -b ${params.results_dir}/pycometh/Meth_Comp.bed -t ${params.results_dir}/pycometh/Meth_Comp.tsv -m ${params.Meth_Comp_m} -l ${params.Meth_Comp_l} --pvalue_adj_method ${params.Meth_Comp_pvalue_adj_method} --pvalue_threshold ${params.Meth_Comp_pvalue_threshold} ${params.Meth_Comp_only_tested_sites}
        else
            pycoMeth Meth_Comp -i \$m5 -a {params.results_dir}/pycometh/CGI_Aggregate.bed -s \$samples_list -f ${params.reference} -b ${params.results_dir}/pycometh/Meth_Comp.bed -t ${params.results_dir}/pycometh/Meth_Comp.tsv -m ${params.Meth_Comp_m} -l ${params.Meth_Comp_l} --pvalue_adj_method ${params.Meth_Comp_pvalue_adj_method} --pvalue_threshold ${params.Meth_Comp_pvalue_threshold} ${params.Meth_Comp_only_tested_sites}
        fi
        
        ln -s ${params.results_dir}/pycometh/Meth_Comp.bed ./Meth_Comp.bed
        ln -s ${params.results_dir}/pycometh/Meth_Comp.tsv ./Meth_Comp.tsv
    """
    else
    """
       ln -s ${params.results_dir}/pycometh/Meth_Comp.bed ./Meth_Comp.bed
       ln -s ${params.results_dir}/pycometh/Meth_Comp.tsv ./Meth_Comp.tsv
    """
}

process pycomethCompReport {
    input:
    file ('Meth_Comp.tsv') from pycomethMethComp_pycomethCompReport
    val(sample) from pycomethMethComp_pycomethCompReport2.collect()
    
    output:

    script:
    if(params.pycomethCompReport)
    """
        mkdir -p ${params.results_dir}/pycometh/
        samples_list=\$(echo "${sample}" | sed \'s/\\[//\' | sed \'s/\\]//\' | sed \'s/,//g\')
        m5=""
        for s in \$samples_list; do
            m5_curr=\$(find ${params.results_dir}/\$s/nanopolish | grep \"\\.m5\$\");
            m5=\$m5\" \"\$m5_curr;
        done

        pycoMeth Comp_Report -i \$m5 -c 'Meth_Comp.tsv' -g ${params.gff} -f ${params.reference} -o ${params.results_dir}/pycometh -n ${params.Comp_Report_n} -d ${params.Comp_Report_d} --pvalue_threshold ${params.Comp_Report_pvalue_threshold} --min_diff_llr ${params.Comp_Report_min_diff_llr} --n_len_bin ${params.Comp_Report_n_len_bin} ${params.Comp_Report_export_static_plots} ${params.Comp_Report_report_non_significant};
        ln -s ${params.results_dir}/pycometh/pycoMeth_summary_report.html ./pycoMeth_summary_report.html;
    """
    else
    """
        echo "Skipped"
    """
}
