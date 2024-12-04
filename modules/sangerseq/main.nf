//
// Assembling Sanger sequencing data into contiguous consensus reads
//

process SANGERSEQ_ANALYSER {
    tag "$sample"
    label 'process_medium'

    input:
    tuple val(sample), file(trace_files)
	
    output:
    tuple val(sample), path("*_qc_metrics.csv")             , emit: qc_report
    tuple val(sample), path("*consensus_sequence.fa")       , emit: contig_fasta
    tuple val(sample), path("*consensus_sequence.txt")      , emit: contig_seq
    tuple val(sample), path("*forward_sequence.txt")        , emit: fwd_seq
    tuple val(sample), path("*reverse_sequence.txt")        , emit: rev_seq
    tuple val(sample), path("*_F.chromatogram.pdf")         , emit: fwd_pherogram
    tuple val(sample), path("*_R.chromatogram.pdf")         , emit: rev_pherogram
    tuple val(sample), path("*_F.chromatogram_100bases.pdf"), emit: fwd_100bp_pherogram
    tuple val(sample), path("*_R.chromatogram_100bases.pdf"), emit: rev_100bp_pherogram
    path "versions.yml"                                     , emit: versions
	

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