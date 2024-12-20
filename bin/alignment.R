# mixed base code map
MIXED_BASE_MAP <- c(
  A="A",
  C="C",
  G="G",
  T="T",
  AC="M",
  AG="R",
  AT="W",
  CG="S",
  CT="Y",
  GT="K",
  ACG="V",
  ACT="H",
  AGT="D",
  CGT="B",
  ACGT="N"
)

process_sanger_trace <- function( input_trace, 
                                  sliding_window_size = 10,
                                  cutoff_qs = 20,
                                  base_num_per_row = 100,
                                  height_per_row=200, 
                                  signal_ratio_cutoff = 0.25) {
  sangerseq_read <- SangerRead(readFeature = ifelse(grepl("_F_",input_trace),
                                                    "Forward Read", "Reverse Read"),
                               readFileName          = input_trace,
                               geneticCode           = GENETIC_CODE,
                               TrimmingMethod        = "M1",
                               M1TrimmingCutoff      = 0.1,
                               M2CutoffQualityScore  = NULL,
                               M2SlidingWindowSize   = NULL,
                               baseNumPerRow         = base_num_per_row,
                               heightPerRow          = height_per_row,
                               signalRatioCutoff     = signal_ratio_cutoff,
                               showTrimmed           = TRUE)

  new_sanger_read <- updateQualityParam(sangerseq_read,
                                        TrimmingMethod       = "M2",
                                        M1TrimmingCutoff     = NULL,
                                        M2CutoffQualityScore = cutoff_qs,
                                        M2SlidingWindowSize  = min(sliding_window_size,
                                                  nchar(sangerseq_read@primarySeq)))

  sangerseq_call <- MakeBaseCalls(new_sanger_read, 
                                  signalRatioCutoff = signal_ratio_cutoff)
}

get_consensus_seq <- function(input_trace_fwd, 
                              input_trace_rev, 
                              qs_fwd, qs_rev,
                              primer_fwd="GAGTTTGATCATGGCTCAG",
                              primer_rev="ACCAGGGTATCTAATCC") {
  # align forward and reverse sanger reads
  sangerseq_call_fwd <- process_sanger_trace(input_trace_fwd)
  sangerseq_call_rev <- process_sanger_trace(input_trace_rev)
  primary_seq_raw_fwd <- as.character(sangerseq_call_fwd@primarySeq)
  primary_seq_raw_rev <- as.character(sangerseq_call_rev@primarySeq)
  primary_seq_trim_fwd <- primary_seq_raw_fwd
  primary_seq_trim_rev <- primary_seq_raw_rev
  if (str_detect(primary_seq_raw_fwd, primer_fwd)) {
    pos_fwd <- unlist(str_locate(pattern =primer_fwd,
                          as.character(primary_seq_raw_fwd)))[2]
    primary_seq_trim_fwd <- substr(primary_seq_raw_fwd,pos_fwd,
                                   sangerseq_call_fwd@QualityReport@trimmedFinishPos)
  } else {
    primary_seq_trim_fwd <- substr(primary_seq_raw_fwd,
                                   sangerseq_call_fwd@QualityReport@trimmedStartPos,
                                   sangerseq_call_fwd@QualityReport@trimmedFinishPos)
  }
  
  if (str_detect(primary_seq_raw_rev, primer_rev)) {
    pos_rev <- unlist(str_locate(pattern =primer_rev,
                                 as.character(primary_seq_raw_rev)))[2]
    primary_seq_trim_rev <- substr(primary_seq_raw_rev,pos_rev,
                                   sangerseq_call_rev@QualityReport@trimmedFinishPos)
  } else {
    primary_seq_trim_rev <- substr(primary_seq_raw_rev,
                                   sangerseq_call_rev@QualityReport@trimmedStartPos,
                                   sangerseq_call_rev@QualityReport@trimmedFinishPos)
  }

  if (qs_fwd < 20 & qs_rev >= 20 ) {
    contig_str <- as.character(reverseComplement(DNAString(primary_seq_trim_rev)))
  } else if (qs_fwd >= 20 & qs_rev < 20) {
    contig_str <- as.character(primary_seq_trim_fwd)
  } else {
    sanger_seq <- DNAStringSet(c(primary_seq_trim_fwd,
                    as.character(reverseComplement(DNAString(primary_seq_trim_rev)))))
    aligned_sanger_seq <- AlignSeqs(sanger_seq)
    seq_fwd <- as.character(aligned_sanger_seq[1])
    seq_rev <- as.character(aligned_sanger_seq[2])
    start_pos <- str_count(substr(seq_fwd, 0, 40),"-")
    end_pos <- str_count(substr(seq_rev, nchar(seq_rev)-40, nchar(seq_rev)),"-")
    contig_str <- substr(seq_fwd, start_pos+1,
                         nchar(seq_fwd)-end_pos)
    for (i in 1:str_count(contig_str, "-")) {
      offset <- unlist(str_locate(pattern = "-",contig_str))[1]
      rev_char <- substr(seq_rev, start = offset+start_pos, stop = offset+start_pos)
      substr(contig_str, start = offset, stop = offset) <- rev_char
    }
  }
  c(contig_str, primary_seq_trim_fwd, primary_seq_trim_rev)
}