#!/usr/bin/env Rscript

#' title: "Wrapper functions for SangerAnalyse in R"
#' author: "kelseyz"
#' date: "October 08, 2024"

library(sangeranalyseR)
library(sangerseqR)
library(Biostrings)
library(DECIPHER)
library(stringr)
library(optparse)
library(R.utils)

options(warn=-1) #turn off warnings

option_list = list(
  make_option(c("-p", "--parent_dir"), type="character", 
              help="Path to trace files"),
  make_option(c("-o", "--output_dir"), type="character", 
              help="Output path"),
  make_option(c("-d", "--basedir"), type="character", 
              help="Path to codebase"),
  make_option(c("-f", "--forward"), type="character", 
              help="Forward trace filename"),
  make_option(c("-r", "--reverse"), type="character",
              help="Reverse trace filename"),
  make_option(c("-c", "--cutoff"), type="integer", default=20, 
              help="Cutoff quality score for trimming"),
  make_option(c("-l", "--seq_len"), type="integer", default=50, 
              help="Minimum sequence length for consensus reads"),
  make_option(c("-a", "--primer_fwd"), type="character", default="AGAGTTTGATCMTGGCTCAG", 
              help="forwad primer"),
  make_option(c("-b", "--primer_rev"), type="character", default="ACCAGGGTATCTAATCC", 
              help="reverse primer")
); 

# read options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
source (file.path(opt$basedir,"bin/chromatogram.R"))
source (file.path(opt$basedir,"bin/alignment.R"))
primer_fwd <- opt$primer_fwd
primer_rev <- opt$primer_rev

get_qc_report <- function(parent_dir, output_dir, sample_trace, 
                          sample_name, cutoff=20, window_size=10) {
  input_trace <- file.path(parent_dir, sample_trace)
  sangerseq_read <- SangerRead(readFeature = ifelse(grepl("_F",input_trace), 
                                                    "Forward Read", "Reverse Read"),
                               readFileName          = input_trace,
                               geneticCode           = GENETIC_CODE,
                               TrimmingMethod        = "M1",
                               M1TrimmingCutoff      = 0.01,
                               M2CutoffQualityScore  = NULL,
                               M2SlidingWindowSize   = NULL,
                               baseNumPerRow         = 100,
                               heightPerRow          = 200,
                               signalRatioCutoff     = 0.25,
                               showTrimmed           = TRUE)
  new_sanger_read <- updateQualityParam(sangerseq_read,
                                        TrimmingMethod       = "M2",
                                        M1TrimmingCutoff     = NULL,
                                        M2CutoffQualityScore = cutoff,
                                        M2SlidingWindowSize  = min(window_size,
                                            nchar(sangerseq_read@primarySeqRaw)))
  sangerseq_call <- MakeBaseCalls(new_sanger_read, signalRatioCutoff = 0.25)
  report <- "@"(sangerseq_call, QualityReport)
  seq_length <- as.integer(report@rawSeqLength)
  mean_qs <- as.integer(report@rawMeanQualityScore)
  qs_20plus <- as.integer(as.data.frame(table(report@qualityPhredScores>20))[2,2])
  signal_strength <- round(mean(apply(sangerseq_call@peakAmpMatrixRaw[,c(1,2)], 1, max)),2)
  
  # Plot Chromatogram
  pdf_file_name <- ifelse(grepl("_F",sample_trace), 
                          paste(sample_name,"_F.chromatogram.pdf",sep=""),
                          paste(sample_name,"_R.chromatogram.pdf",sep=""))
  pherogram(sangerseq_call, width = 100, height = 2, 
            trim5 = ifelse(str_detect(as.character(sangerseq_call@primarySeqRaw),
                                      primer_fwd), nchar(primer_fwd), 
                           ifelse (str_detect(as.character(sangerseq_call@primarySeqRaw), 
                                              primer_rev), nchar(primer_rev), 0)),
            trim3 = report@rawSeqLength-report@trimmedFinishPos,
            cex.mtext = 0.5, cex.base = 0.5,
            showcalls = "both",
            filename = file.path(".", pdf_file_name),
            showtrim=TRUE, showhets=TRUE)
  if ( seq_length >= 50 ) {
    pherogram(sangerseq_call, width = 100, height = 2, trim5 = ifelse(seq_length > 50, 50, seq_len),
                 trim3 = ifelse(seq_length-150 < 0, 0, seq_length-150),
                 cex.mtext = 0.5, cex.base = 0.5,
                 showcalls = "primary",
                 filename = file.path(".", gsub("chromatogram",
                                    "chromatogram_100bases",pdf_file_name)))
  }
  # Output raw sequences to a FASTA file
  fasta_file_name <- ifelse(grepl("_F",sample_trace), 
                            paste(sample_name,"_F.raw.fa",sep=""),
                            paste(sample_name,"_R.raw.fa",sep=""))
  writeXStringSet(DNAStringSet(c(sangerseq_call@primarySeqRaw)), 
                  filepath = file.path(output_dir, sample_name, fasta_file_name),
                  format = "fasta")
  c(seq_length,mean_qs,qs_20plus,signal_strength)
}

call_sangerseq <- function(parent_dir, output_dir, sample_trace_fwd, sample_trace_rev, cutoff=20, seq_len=50) {
  sample_name <- unlist(strsplit(sample_trace_fwd, "_F"))[1]
  sample_dir <- file.path(output_dir,sample_name)
  
  if (!file.exists(sample_dir)){
    dir.create(sample_dir)
  }
  
  # Extracting ab1 files of the interest sequence (16s)
  # Separating forward and reverse sequences in two distinct data collections
  qc_report_fwd <- get_qc_report(parent_dir, output_dir, sample_trace_fwd, sample_name)
  qc_report_rev <- get_qc_report(parent_dir, output_dir, sample_trace_rev, sample_name)
  print(paste(qc_report_fwd[2], qc_report_rev[2],sep=":"))
  contig_str <- get_consensus_seq(file.path(parent_dir, sample_trace_fwd), 
                                  file.path(parent_dir, sample_trace_rev),
                                  qc_report_fwd[2], qc_report_rev[2],
                                  primer_fwd, primer_rev)
  
  # Output assembled contiguous consensus sequence to a TXT file
  fasta_suffix <- c("consensus", "forward","reverse")
  match_rate = list()
  for (i in 1:length(contig_str)) {
    match_rate[[i]] <- round(str_count(contig_str[i],"T|C|G|A")/nchar(contig_str[i]),3)
    file_out <-file(paste("./",paste(sample_name,fasta_suffix[i],"sequence.txt",sep="_"),sep="/"))
    writeLines(contig_str[i], file_out)
    close(file_out)
  }
  
  # Save the consensus sequence as FASTA format for NCBI blast search
  fasta_file_name <- paste(unlist(strsplit(sample_name, "-"))[1],"consensus_sequence.fa",sep="_")
  writeXStringSet(DNAStringSet(c(contig_str)), 
                  filepath = file.path("./",fasta_file_name),
                  format = "fasta")
  
  # In order to reconstitute the entire sequence of the 16S rRNA gene from both strands 
  # the .ab1 files need to be assembled, i.e., an overlap needs to be established 
  # and a contiguous consensus sequence (“contig”) needs to be generated.
  # output quality control metrics to a CSV
  file_logger <- file(paste("./", paste(sample_name,"qc_metrics.csv",sep="_"),sep="/"),"a")
  writeLines(paste("Trace File Name","Read Length","Mean QS","QS20+","Contig Length",
                   "Signal Strength","Mismatch%",sep=","),file_logger)
  writeLines(paste(sample_trace_fwd,qc_report_fwd[1],
                     qc_report_fwd[2],qc_report_fwd[3],nchar(contig_str[1]),
                     qc_report_fwd[4],(1-match_rate[[2]]),sep=","),file_logger)
  writeLines(paste(sample_trace_rev,qc_report_rev[1],
                     qc_report_rev[2],qc_report_rev[3],nchar(contig_str[1]),
                     qc_report_rev[4],(1-match_rate[[3]]),sep=","),file_logger)
  close(file_logger)
}

call_sangerseq(opt$parent_dir, opt$output_dir, opt$forward, opt$reverse, opt$cutoff, opt$seq_len)