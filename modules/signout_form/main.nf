// blastn
process SIGNOUT_FORM {
    tag "$sample"
    label 'process_medium'

	input:
	tuple val(sample), path(html)

	output:
	tuple val(sample), path("*.docx"), emit: docx
    path "versions.yml"     		 , emit: versions

	script:
	"""
	${baseDir}/bin/gen_signout_form.py -s ${sample}
	
	cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PYTHON: \$(echo \$(python3 -V 2>&1))
    END_VERSIONS
	"""
}