//
// Assembling Sanger sequencing data into contiguous consensus reads
//

process SANGERSEQ_ANALYSER {
    tag "$sample"
    label 'process_medium'
    //errorStrategy 'ignore'

    input:
    tuple val(sample), file(trace_files)
	
    output:
    path "versions.yml"             , emit: versions
	

    //Calculate tumor mutation burden per sample
    script:
    def batch_id   = params.batch_id ?: ""
    def trace_path = params.trace_path ?: ""
    """
    Rscript ${baseDir}/bin/sangerseq_batch_analyse.R \
         -p ${trace_path}/${batch_id} \
         -f ${trace_files[0]} \
         -r ${trace_files[1]} \
         -c ${params.trim_cutoff} \
         -l ${params.min_seq_len}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version 2>&1))
    END_VERSIONS
    """
}