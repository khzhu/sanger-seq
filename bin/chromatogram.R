pherogram <- function(obj, trim5=0, trim3=0, 
                         showcalls=c("primary", "secondary", "both", "none"), 
                         width=100, height=2, cex.mtext=0.8, cex.base=0.8, ylim=3, 
                         filename=NULL, showtrim=FALSE, showhets=TRUE) {
  originalpar <- par(no.readonly=TRUE)
  showcalls <- showcalls[1]
  traces <- obj@traceMatrix
  basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
  basecalls2 <- unlist(strsplit(toString(obj@secondarySeq), ""))
  aveposition <- rowMeans(obj@peakPosMatrix, na.rm=TRUE)
  basecalls1 <- basecalls1[1:length(aveposition)] 
  basecalls2 <- basecalls2[1:length(aveposition)] 
  if(showtrim == FALSE) {
    if(trim5+trim3 > length(basecalls1)) basecalls1 <- ""
    else basecalls1 <- basecalls1[(1 + trim5):(length(basecalls1) - trim3)]
    if(trim5+trim3 > length(basecalls2)) basecalls2 <- ""
    else basecalls2 <- basecalls2[(1 + trim5):(length(basecalls2) - trim3)]
    aveposition <- aveposition[(1 + trim5):(length(aveposition) - trim3)] 
  }
  indexes <- 1:length(basecalls1)
  trimmed <- indexes <= trim5 | indexes > (length(basecalls1) - trim3) # all 
  #false if not trimmed
  if (!is.null(trim3)) {
    traces <- traces[1:(min(max(aveposition, na.rm=TRUE) + 10, 
                            nrow(traces))), ]
  }
  if (!is.null(trim5)) {
    offset <- max(c(1, aveposition[1] - 10))
    traces <- traces[offset:nrow(traces),]
    aveposition <- aveposition - (offset-1)
  }
  maxsignal <- apply(traces, 1, max)
  ylims <- c(0, quantile(maxsignal, .75)+ylim*IQR(maxsignal))           
  p <- c(0, aveposition, nrow(traces))
  midp <- diff(p)/2
  starts <- aveposition - midp[1:(length(midp)-1)]
  starthets <- starts
  starthets[basecalls1 == basecalls2] <- NA
  ends <- aveposition + midp[2:(length(midp))]
  endhets <- ends
  endhets[basecalls1 == basecalls2] <- NA
  starttrims <- starts
  starttrims[!trimmed] <- NA
  endtrims <- ends
  endtrims[!trimmed] <- NA
  
  colortranslate <- c(A="green", C="blue", G="black", T="red")
  colorvector1 <- unname(colortranslate[basecalls1])
  colorvector1[is.na(colorvector1)] <- "purple"
  colorvector2 <- unname(colortranslate[basecalls2])
  colorvector2[is.na(colorvector2)] <- "purple"
  
  valuesperbase <- nrow(traces)/length(basecalls1)
  tracewidth <- width*valuesperbase
  breaks <- seq(1,nrow(traces), by=tracewidth)
  numplots <- length(breaks)
  if(!is.null(filename)) pdf(filename, width=8.5, height=height*numplots)
  # relative heights
  # margin values in the following order: bottom, left, top, right as value
  par(mar=c(0.5,1.8,2,0), mfrow=c(numplots, 1))
  basecallwarning1 = 0
  basecallwarning2 = 0
  j = 1
  for(i in breaks) {
    range <- aveposition >= i & aveposition < (i+tracewidth)
    starthet <- starthets[range] - tracewidth*(j-1)
    starthet[starthet < 0] <- 0
    endhet <- endhets[range] - tracewidth*(j-1)
    endhet[endhet > tracewidth] <- tracewidth
    lab1 <- basecalls1[range]
    lab2 <- basecalls2[range]
    pos <- aveposition[range] - tracewidth*(j-1)
    colors1 <- colorvector1[range]
    colors2 <- colorvector2[range]
    starttrim <- starttrims[range] - tracewidth*(j-1)
    endtrim <- endtrims[range] - tracewidth*(j-1)
    plotrange <- i:min(i+tracewidth, nrow(traces))
    bar_start = which(range)[1]
    bar_end = bar_start+length(lab1)-1
    quality_scores <- obj@QualityReport@qualityPhredScores[bar_start:bar_end]
    plot(traces[plotrange,1], type='n', ylim=ylims, ylab="", xaxt="n", 
         bty="n", xlab="", yaxt="n", , xlim=c(1,tracewidth))
    axis(1, at = seq(0, tracewidth, by = valuesperbase*10), las=2, labels = FALSE, cex=0.3)
    rug(x = seq(0, tracewidth, by = valuesperbase), ticksize = -0.02, side = 1)
    if (showhets==TRUE) {
      rect(starthet, 0, endhet, ylims[2], col='#D5E3F7', border='#D5E3F7')
    }
    if (showtrim==TRUE) {
      rect(starttrim, 0, endtrim, ylims[2], col='gray70', border='transparent', 
           density=15)
    }
    lines(traces[plotrange,1], col="green")
    lines(traces[plotrange,2], col="blue")
    lines(traces[plotrange,3], col="black")
    lines(traces[plotrange,4], col="red")
    mtext(paste("position",paste(as.character(bar_start),
          as.character(bar_end),sep="~"),sep=":"),side=2, line=0, cex=cex.mtext)
    for(k in 1:length(lab1)) {
      if (showcalls=="primary" | showcalls=="both") {
        if (is.na(basecalls1[1]) & basecallwarning1==0) {
          warning("Primary basecalls missing")
          basecallwarning1 = 1
        } 
        else if (length(lab1) > 0) {
          axis(side=3, at=pos[k], labels=quality_scores[k], 
               col.axis=colors1[k], cex.axis=0.48,
               family="mono", cex=cex.base, line=0, tick=FALSE,font = 2)
          axis(side=3, at=pos[k], labels=lab1[k], col.axis=colors1[k], cex.axis=0.9,
               family="mono", cex=cex.base, line=-0.5, tick=FALSE)
          if (quality_scores[k] < 10) {
            rect(pos[k], 0, pos[k], ylims[2], 
                 col=NULL, border='red', lwd = 1, lty=3)
          }
        }
      }
      if (showcalls=="secondary" | showcalls=="both") {
        if (is.na(basecalls2[1]) & basecallwarning2 == 0) {
          warning("Secondary basecalls missing")
          basecallwarning2 = 1
        } 
        else if (length(lab2) > 0) { 
          axis(side=3, at=pos[k], labels=lab2[k], col.axis=colors2[k], cex.axis=0.9,
               family="mono", cex=cex.base, line=-1.0, tick=FALSE)
        }
      }
    }
    j = j + 1
  }
  if(!is.null(filename)) {
    dev.off()
    message(paste("Chromatogram saved to", filename, 
                  "in the current working directory"))
  }
  else par(originalpar)
}