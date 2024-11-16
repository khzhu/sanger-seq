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
  sangerseq_read <- SangerRead(readFeature = ifelse(grepl("_F",input_trace),
                                                    "Forward Read", "Reverse Read"),
                               readFileName          = input_trace,
                               geneticCode           = GENETIC_CODE,
                               TrimmingMethod        = "M1",
                               M1TrimmingCutoff      = 0.01,
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
                                                  nchar(sangerseq_read@primarySeqRaw)))

  sangerseq_call <- MakeBaseCalls(new_sanger_read, 
                                  signalRatioCutoff = signal_ratio_cutoff)
}

get_consensus_seq <- function(input_trace_fwd, 
                              input_trace_rev, 
                              qs_fwd, qs_rev,
                              primer_fwd="AGAGTTTGATCMTGGCTCAG",
                              primer_rev="ACCAGGGTATCTAATCC") {
  # align forward and reverse sanger reads
  sangerseq_call_fwd <- process_sanger_trace(input_trace_fwd)
  sangerseq_call_rev <- process_sanger_trace(input_trace_rev)
  primary_seq_raw_fwd <- as.character(sangerseq_call_fwd@primarySeqRaw)
  primary_seq_raw_rev <- as.character(sangerseq_call_rev@primarySeqRaw)
  primary_seq_trim_fwd <- primary_seq_raw_fwd
  primary_seq_trim_rev <- primary_seq_raw_rev
  if (str_detect(primary_seq_raw_fwd, primer_fwd)) {
    pos_fwd <- unlist(gregexpr(primer_fwd, primary_seq_raw_fwd))[1]
    primary_seq_trim_fwd <- substr(primary_seq_raw_fwd,pos_fwd+nchar(primer_fwd),
                                   nchar(primary_seq_raw_fwd))
  }
  
  if (str_detect(primary_seq_raw_rev, primer_rev)) {
    pos_rev <- unlist(gregexpr(primer_rev, primary_seq_raw_rev))[1]
    primary_seq_trim_rev <- substr(primary_seq_raw_rev,pos_rev+nchar(primer_rev),
                                   nchar(primary_seq_raw_rev))
  }

  if (qs_fwd < 20 & qs_rev > 20) {
    rev_comp_str <- as.character(reverseComplement(DNAString(primary_seq_trim_rev)))
    contig_str <- substr(rev_comp_str, nchar(primer_fwd)*2+1, nchar(rev_comp_str))
  } else {
    sanger_seq <- DNAStringSet(c(primary_seq_trim_fwd,
                    as.character(reverseComplement(DNAString(primary_seq_trim_rev)))))
    aligned_sanger_seq <- AlignSeqs(sanger_seq)
    seq_fwd <- as.character(aligned_sanger_seq[1])
    seq_rev <- as.character(aligned_sanger_seq[2])
    contig_str <- substr(seq_fwd, str_count(seq_fwd,"-")+1,nchar(seq_fwd)-str_count(seq_rev,"-"))
  }
}