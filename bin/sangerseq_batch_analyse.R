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

options(warn=-1) #turn off warnings

option_list = list(
  make_option(c("-p", "--parent_dir"), type="character", 
              help="Path to trace files"),
  make_option(c("-f", "--forward"), type="character", 
              help="Forward trace filename"),
  make_option(c("-r", "--reverse"), type="character",
              help="Reverse trace filename"),
  make_option(c("-c", "--cutoff"), type="integer", default=0, 
              help="Cutoff quality score for trimming"),
  make_option(c("-l", "--seq_len"), type="integer", default=50, 
              help="Minimum sequence length for consensus reads")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

get_qc_report <- function(parent_dir, sample_trace, sample_name) {
  input_trace <- file.path(parent_dir, sample_trace)
  sangerseq_read <- SangerRead(readFeature = ifelse(grepl("_F",input_trace), 
                                              "Forward Read", "Reverse Read"),
                               readFileName          = input_trace,
                               geneticCode           = GENETIC_CODE,
                               TrimmingMethod        = "M1",
                               M1TrimmingCutoff      = 0.0001,
                               M2CutoffQualityScore  = NULL,
                               M2SlidingWindowSize   = NULL,
                               baseNumPerRow         = 100,
                               heightPerRow          = 200,
                               signalRatioCutoff     = 0.33,
                               showTrimmed           = TRUE)
  sangerseq_call <- MakeBaseCalls(sangerseq_read, signalRatioCutoff = 0.33)
  report <- "@"(sangerseq_call, QualityReport)
  seq_length <- as.integer(report@rawSeqLength)
  mean_qs <- as.integer(report@rawMeanQualityScore)
  qs_tab <- as.data.frame(table(report@qualityPhredScores>20))
  # Output raw sequences to a FASTA file
  writeFasta(sangerseq_call,
           outputDir         = file.path("/mnt",sample_name), #tempdir(),
           compress          = FALSE,
           compression_level = NA)
  c(seq_length,mean_qs,as.integer(qs_tab[2,2]))
}

process_sample_trace <- function(parent_dir, sample_trace, sample_name, read_length, seq_len) {
  # Extracting ab1 files of the interest sequence (e.g 16s)
  # Separating forward and reverse sequences in two distinct data collections
  input_trace <- file.path(parent_dir, sample_trace)
  sanger_read <- readsangerseq(input_trace)
  sanger_call <- makeBaseCalls(sanger_read, ratio = 0.33)
  pdf_file_name <- ifelse(grepl("_F",sample_trace), 
                      paste(sample_name,"_F.chromatogram.pdf",sep=""),
                      paste(sample_name,"_R.chromatogram.pdf",sep=""))
  output_dir = paste(parent_dir,"_sanger_seq_output",sep="")
  chromatogram(sanger_call, width = 100, height = 2, trim5 = 0, trim3 = 0,
               cex.mtext = 0.5, cex.base = 0.5,
               showcalls = "primary",
               filename = file.path(output_dir, sample_name, pdf_file_name))
  if ( read_length >= seq_len ) {
    chromatogram(sanger_call, width = 100, height = 2, trim5 = ifelse(seq_len > 50, 50, seq_len), 
               trim3 = ifelse(read_length-150 < 0, 0, read_length-150),
               cex.mtext = 0.5, cex.base = 0.5,
               showcalls = "primary",
               filename = file.path(output_dir, sample_name, 
                                    gsub("chromatogram","chromatogram_100bases",pdf_file_name)))
  }
  return (sanger_call)
}

call_sangerseq <- function(parent_dir, sample_trace_fwd, sample_trace_rev, cutoff, seq_len) {
  sample_name <- unlist(strsplit(sample_trace_fwd, "_F"))[1]
  sample_dir <- file.path('/mnt',sample_name)
  
  if (!file.exists(sample_dir)){
    dir.create(sample_dir)
  }
  
  # output quality control metrics to a CSV
  file_logger <- file(file.path(sample_dir,
                  paste(sample_name,"qc_metrics.csv",sep="_")),"a")
  writeLines(paste("Trace File Name","Read Length","Mean QS","QS20+","Contig Length",
               "Mismatch%","Match%",sep=","),file_logger)

  # Extracting ab1 files of the interest sequence (16s)
  # Separating forward and reverse sequences in two distinct data collections
  qc_report_fwd <- get_qc_report(parent_dir, sample_trace_fwd, sample_name)
  sanger_call_fwd <- process_sample_trace(parent_dir, sample_trace_fwd, 
                                sample_name, qc_report_fwd[1], seq_len)
  qc_report_rev <- get_qc_report(parent_dir, sample_trace_rev, sample_name)
  sanger_call_rev <- process_sample_trace(parent_dir, sample_trace_rev, 
                                sample_name, qc_report_rev[1], seq_len)

  
  # Compute the reverse-complement for the antisense sequence only
  reverse_seq <- DNAString(primarySeq(sanger_call_rev, string = TRUE))
  
  # Align forward and reverse sequences
  sanger_seq <- DNAStringSet(c(primarySeq(sanger_call_fwd, string = TRUE), 
                               as.character(reverseComplement(reverse_seq))))
  aligned_rna_seq <- AlignSeqs(sanger_seq)
  sense_seq <- as.character(aligned_rna_seq[1])
  antisense_seq <- as.character(aligned_rna_seq[2])
  
  # In order to reconstitute the entire sequence of the 16S rRNA gene from both strands 
  # the .ab1 files need to be assembled, i.e., an overlap needs to be established 
  # and a contiguous consensus sequence (“contig”) needs to be generated.
  if (qc_report_fwd[1] > seq_len & qc_report_rev[1] > seq_len) {
    my_sanger_contig <- SangerContig(
      inputSource = "ABIF",
      ABIF_Directory = parent_dir,
      contigName = sample_name,
      REGEX_SuffixForward = paste("_F",unlist(strsplit(sample_trace_fwd, "_F"))[2],sep=""),
      REGEX_SuffixReverse = paste("_R",unlist(strsplit(sample_trace_rev, "_R"))[2],sep=""),
      refAminoAcidSeq = "",
      TrimmingMethod = "M2",
      M1TrimmingCutoff = NULL,
      M2CutoffQualityScore = cutoff,
      M2SlidingWindowSize = 10,
      baseNumPerRow = 100,
      heightPerRow = 200,
      signalRatioCutoff = 0.33,
      showTrimmed = TRUE,
      processorsNum = 4)
    
    contig_seq <- as.character(my_sanger_contig@contigSeq)
    cat(contig_seq,file=file.path(sample_dir,
                          paste(sample_name,"consensus_sequence.txt",sep="_")))
    no_match_num <- str_count(sense_seq, pattern = "-")
    mis_match_rate <- round(no_match_num/nchar(contig_seq)*100,2)

    writeLines(paste(sample_trace_fwd,qc_report_fwd[1],
                qc_report_fwd[2],qc_report_fwd[3],nchar(contig_seq),
                mis_match_rate,(100-mis_match_rate),sep=","),file_logger)
    writeLines(paste(sample_trace_rev,qc_report_rev[1],
                qc_report_rev[2],qc_report_rev[3],nchar(contig_seq),
                mis_match_rate,(100-mis_match_rate),sep=","),file_logger)
  
    # Output assembled contiguous consensus sequence to a FASTA file
    writeFasta(my_sanger_contig,
              outputDir         = file.path(sample_dir),
              compress          = FALSE,
              compression_level = NA,
              selection         = "all")
  } else { #to cover negative controls
    writeLines(paste(sample_trace_fwd,qc_report_fwd[1],
                qc_report_fwd[2],qc_report_fwd[3],'NA',
                'NA','NA',sep=","),file_logger)
    writeLines(paste(sample_trace_rev,qc_report_rev[1],
                qc_report_rev[2],qc_report_rev[3],'NA',
                'NA','NA',sep=","),file_logger)
  }
  close(file_logger)
}

call_sangerseq(opt$parent_dir, opt$forward, opt$reverse, opt$cutoff, opt$seq_len)