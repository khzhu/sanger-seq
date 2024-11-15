# Workflow for assembling sanger seuqencing data into contigs

## Introduction
Accurate identification of infectious bacteria is a major challenge for clinical practice. The 16S ribosomal RNA (rRNA) gene seqeuencing is usually chosen as the preferred method for taxonomic classification and identification of bacteria at the species level. sangeranalyseR and sangerseqR are among many open-source tools for processing Sanger sequencing data.
sangeranalyseR provides a wide range of options for trimming reads, detecting secondary peaks, viewing chromatograms, aligning contigs, and outputs aligned and unaligned reads and contigs in FASTA format.

## Pipeline
In order to reconstitute the entire sequence of the 16S rRNA gene from both strands, the .ab1 files need to be assembled, i.e.,
an overlap needs to be established and a contiguous consensus sequence (“contig”) needs to be generated. We are here to present a nextflow pipeline that assembles Sanger sequencing data into contiguous consensus reads. It performs the following operations on Sanger sequencing data:

1. Extract ab1 files of the interest sequence (16S)
2. Separate sense and antisense sequences in two distinct data collections
3. Convert ab1 files to FASTQ to permit its use by other software tools
4. Remove the primers from the sequence
5. Trim low quality ends of the sequences using trimmomatic’s sliding window method
6. Compute the reverse-complement for the antisense sequence only
7. Align sense and antisense sequences
8. Assemble a contiguous sequence for microbial nucleotide BLAST search
9. Generate quality control metrics
10. Perform a BLAST to identify bacterial species based on 16S rRNA sequences

## Trimming methods

M2 is like trimmomatic’s sliding window method. The method will cut N bases off, if average coverage of quality
for N bases is lower than the quality score (Q). M2CutoffQualityScore and M2SlidingWindowSize are
two parameters that control M2 trimming. Their default values are “20” and “10” respectively.
It cuts the read when the average quality of each 10-nt window falls below a quality score of 20.

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
You will need to install a nextflow and to pull a docker image from Azure CR prior to running the workflow.
- Input: the 16S rRNA gene sequencing data from both strands in a ABIF format, input trace file suffix, batch ID of a sequencing run, and the output directory.
- Output: assembled consensus sequence, text-based FASTA files that store raw/trimmed sequencing reads, two-dimensional plots with the ordinate axis giving signal strength, and QC metrics.

To run this pipeline, enter
```
nextflow run main.nf
        --batch_id 16s_test1
        --trace_path /Users/kz347/workspace/data/Sanger-Seq-Results
        --trace_regex_suffix 2024-10-23-14-21-48_sb.ab1
        --output_dir /Users/kz347/workspace/data/Sanger-Seq-Results/16s_test1_sanger_seq_output
        --blastn_db "rRNA_typestrains/16S_ribosomal_RNA"
        -bg
```
from the project root directory.
## Directory Structure
```
sanger-seq/
├── README.md
├── nextflow.config
├── bin
│   ├── alignment.R
│   ├── chromatogram.R
│   └── sangerseq_batch_process.R
├── docker-image
│   ├── Dockerfile
│   └── requirements.txt
├── main.nf
├── modules
│   ├── blastn
│   │   ├── main.nf
│   │   └── meta.yml
│   └── sangerseq
│       ├── main.nf
│       └── meta.yml
└── workflows
    ├── blastn_search
    │   ├── main.nf
    │   ├── meta.yml
    │   └── nextflow.config
    └── sangerseq_batch
        ├── main.nf
        ├── meta.yml
        └── nextflow.config
```