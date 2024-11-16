// blastn
process BLASTN_HTML {
    tag "$sample"
    label 'process_medium'

	input:
	tuple val(sample), path(tsv_doc)
    tuple val(sample), path(html_doc)

	output:
	tuple val(sample), path("*_identity.html") , emit: html
    path "versions.yml"     				   , emit: versions

	script:
	"""
	${baseDir}/bin/reform_blast_html.py --tsv ${tsv_doc} --html ${html_doc}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PYTHON: \$(echo \$(python3 -V 2>&1))
    END_VERSIONS
	"""
}