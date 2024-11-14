//
// Assembling Sanger sequencing data into contiguous consensus reads
//

process SANGERSEQ_ANALYSER {
    tag "$sample"
    label 'process_medium'

    input:
    tuple val(sample), file(trace_files)
	
    output:
    tuple val(sample), path("*consensus_sequence.fa"), emit: contig_fasta
    path "versions.yml"                              , emit: versions
	

    //Generate contiguous consensus sequence and QC metrics
    script:
    def batch_id   = params.batch_id ?: ""
    def trace_path = params.trace_path ?: ""
    """
    Rscript ${baseDir}/bin/sangerseq_batch_process.R \
         -p ${trace_path}/${batch_id} \
         -f ${trace_files[0]} \
         -r ${trace_files[1]} \
         -c ${params.trim_cutoff} \
         -l ${params.min_seq_len} \
         -d ${baseDir} \
         -o ${params.output_dir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version 2>&1))
    END_VERSIONS
    """
}