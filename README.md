# Workflow for assembling sanger seuqencing data into contigs

## Introduction
Accurate identification of infectious bacteria is a major challenge for clinical practice. The use of sequencing of the 16S ribosomal RNA (rRNA) gene has emerged as the preferred method for taxonomic classification and identification of bacteria at the species level. sangeranalyseR and sangerseqR are open-source R packages for processing Sanger sequencing data sangeranalyseR and sangerseqR are open-source R packages for processing Sanger sequencing data. sangeranalyseR provides a wide range of options for trimming reads, detecting secondary peaks, viewing chromatograms, aligning contigs, and outputs aligned and unaligned reads and contigs in FASTA format. Input data can be in either ABIF or FASTA format.

## Pipeline
1. Extracting ab1 files of the interest sequence (16S)
2. Separating sense and antisense sequences in two distinct data collections
3. Converting ab1 files to FASTQ to permit its use by other software tools
4. Trimming low quality ends of the sequences
5. Compute the reverse-complement for the antisense sequence only
6. Align sense and antisense sequences
7. Obtain a consensus sequence (which results the correspondance between nucleotides of the sense and the antisense sequences) for Blast search
8. Generate quality control metrics

## Software required:
1. [Nextflow](https://www.nextflow.io/docs/latest/)
Nextflow is a workflow engine for creating scalable, portable, and reproducible workflows.
2. [Python](https://www.python.org/)
Python is a object oriented and general-purpose programming language.
3. [R](https://www.r-project.org/)
R is a free software environment for statistical computing and graphics.
4. [SangerAnalyseR](https://sangeranalyser.readthedocs.io/en/latest/content/quickstart.html)
sangeranalseR is an R package that provides fast, flexible, and reproducible workflows for assembling your sanger seuqencing data into contigs.
5. [sangerseqR](https://github.com/jonathonthill/sangerseqR)
sangerseqR is an R packagefor analyzing Sanger Sequencing data files in R, including reading .ab1 files, making basecalls and plotting chromatograms.

## How to set up and run a workflow  
You will need to install a nextflow and docker prior to running the workflow.
- Input: the 16S rRNA gene sequencing data from both strands in a ABIF format and the batch ID of a sequencing run.
- Output: assembled consensus sequence, text-based FASTA files that store raw/trimmed sequencing reads, two-dimensional plots with the ordinate axis giving signal strength, and QC metrics.

To run this pipeline, enter 
```
nextflow run main.nf 
    --batch_id "2024_07_10_2024-07-10-15-00-51"  
    --trace_path "~/workspace/projects/sanger_seq"
    --trace_regex_suffix "2024-07-10-15-00-51_JO.ab1"
    --output_dir "~/workspace/projects/sanger_seq/2024_07_10_2024-07-10-15-00-51_sanger_seq_output"
``` 
from the project root directory.
## Directory Structure
```
sanger-seq
├── README.md
├── bin
│   └── sangerseq_batch_analyse.R
├── docker
│   ├── Dockerfile
│   └── requirements.txt
├── main.nf
├── modules
│   └── sangerseq
│       ├── main.nf
│       └── meta.yml
├── nextflow.config
└── workflows
    └── sangerseq_batch
        ├── main.nf
        ├── meta.yml
        └── nextflow.config
```