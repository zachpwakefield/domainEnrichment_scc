## function to obtain transcript, exon match for each exon input

getTranscript <- function(gtf = gtf, redExon = redExon, ex_type = ex_type, minOverlap = .5, swaps = F, cores = 8) {
  ## Start up Message
  print(paste("searching for ", ex_type, "...", sep = ""))
  rowOuts <- list()
  checks <- c()
  rc_out <- mclapply(1:length(redExon$geneR), mc.cores = cores, function(i) {
    hyb_stat <- "no" #HFE" "HLE" "no"
    
    
    ## Progress Messagge
    if ((i %% max(1, round(length(redExon[,1])/10, 0))) == 0) {
      print(paste(i, " exons or ", i/length(redExon[,1]), " completed", sep=""))
    }
    
    ## Reduce gtf for faster computation
    gtf_min <- gtf[gtf$geneID == redExon$geneR[i] & gtf$type == "exon" & 
                     gtf$classification %in% if (ex_type == "AFE") {c("first", "single_exon")} 
                   else if (ex_type == "ALE") {c("last", "single_exon")} 
                   else if (ex_type == "SE") {c("internal")},]
    
    
    ## Calculate Jaccard
    inEx <- seq(redExon$start[i], redExon$stop[i])
    gtfEx <- lapply(1:length(gtf_min$geneID), function(x) seq(gtf_min$start[x], gtf_min$stop[x]))
    un <- unlist(lapply(gtfEx, function(x) length(union(inEx, unlist(x)))))
    ins <- unlist(lapply(gtfEx, function(x) length(intersect(inEx, unlist(x)))))
    gtf_min$length_jacc <- ins/lengths(gtfEx)
    gtf_min$jaccard <- ins/un
    gtf_min <- gtf_min %>% arrange(desc(jaccard))
    
    ## Proceed through various checks to make sure matches are both identified and valid
    
    ## If no matches
    if (dim(gtf_min)[1] == 0) {checks <- 1
    rowOuts <- 0}
    
    ## If all bad matches
    else if (max(gtf_min$jaccard) < minOverlap) {checks <- 2
    rowOuts <- 0}
    
    ## If only one match
    else if (dim(gtf_min)[1] == 1) {checks <- 2
    rowOuts <- gtf_min$rownum[1]}
    
    ## If one best match
    else if (gtf_min$jaccard[1] > gtf_min$jaccard[2]) {
      checks <- 3
      rowOuts <- gtf_min$rownum[1]
    }
    
    ## If multiple best, use length of matched exon
    else if (gtf_min$jaccard[1] == gtf_min$jaccard[2]) {
      checks <- 4
      c_gtf <- gtf_min[gtf_min$jaccard == max(gtf_min$jaccard),] %>% arrange(desc(length_jacc))
      rowOuts <- c_gtf$rownum[1]
    }
    
    ## Other
    else {checks <- 5
    rowOuts <- 0}
    
    c(rowOuts, checks)
  })
  rowOuts <- unlist(lapply(rc_out, "[[", 1))
  checks <- unlist(lapply(rc_out, "[[", 2))
  ## Remove erroneous matches
  out_matched <- gtf[unlist(rowOuts)[unlist(rowOuts) != 0],]
  
  ## Data arrangement
  out_matched$input_id <- paste(redExon$geneR, ";", redExon$chr, ":", redExon$start, "-", redExon$stop, sep = "")[unlist(rowOuts) != 0]
  out_matched <- out_matched %>% relocate(input_id)
  tot_matched <- out_matched
  
  ## If looking for swaps, extract paired transcripts
  if (swaps == T) {
    out_matched <- out_matched %>% arrange(geneID) %>% dplyr::filter(geneID %in% names(table(out_matched$geneID))[table(out_matched$geneID) == 2])
  }
  return(list(out_matched = out_matched,
              typeBreakdown = checks,
              tot_matched = tot_matched
  ))
}

bedify <- function(matched = matched, saveBED = F, outname = outname, cores = 8) {
  tID <- matched[[1]]$transcriptID
  eiID <- matched[[1]]$input_id
  toBed <- list()
  toBed <- mclapply(1:length(tID), mc.cores = cores, function(i) {
    bed <- gtf[gtf$transcriptID == tID[i] & gtf$type == "exon",] %>% dplyr::select(chr, start, stop, transcriptID, geneID, strand)
    bed$eiID <- rep(eiID[i], length(gtf$geneID[gtf$transcriptID == tID[i] & gtf$type == "exon"]))
    bed$score <- 0
    bed$squish <- paste(bed$transcriptID, "#", bed$eiID, sep = "")
    
    colnames(bed) <- c('chrom','chromStart', 'chromEnd', 'transcriptID', "geneID", "strand", "eiID", "score", "name")
    bed <- bed %>% dplyr::select(chrom, chromStart, chromEnd, name, score, strand)
    if (unique(bed$strand) == "+") {
      bed$chromStart <- as.integer(bed$chromStart) - 1
    } else {
      bed$chromStart <- as.integer(bed$chromEnd) + 1
    }
    bed
  })
  toBed <- do.call(rbind, toBed)
  if (saveBED == T) {
    write.table(toBed, paste("./", outname, ".bed", sep = ""), sep = '\t', col.names = F, quote = F, row.names = F)
  }
  return(toBed)
}


proteinExtract_pipe <- function(files_dir, background = T, mOverlap = .5, saveOutput = F, inCores = 8, nC = 0, nE = 0, exon_type = "AFE") {
  
  if (background == T) {
    files <- paste(files_dir, list.files(files_dir)[grep('[.]exon', list.files(files_dir))], sep = "")
    cat(files)
    first_exons <- unique(unlist(lapply(files, function(x) {
      in_file <- read.delim(x)
      in_file <- in_file[in_file$ID == 'first',]
      paste(in_file$gene, ';', in_file$exon, ';',  in_file$strand, sep = "")})))
    redExon <- data.frame(geneR = unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 1)), split = '[.]'), "[[", 1)),
                          chr = unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 2)), split = '-'), "[[", 1)), split = ":"), "[[", 1)),
                          start = unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 2)), split = '-'), "[[", 1)), split = ":"), "[[", 2)),
                          stop = unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 2)), split = '-'), "[[", 2))
    )
  } else {
    df <- read.delim(files_dir, sep = " ")
    df.l <- lfc(df, numCont = nC, numExp = nE, exon_type = ex_type, cores = inCores)
    lfcPlot <- make_lfcPlot(df.l)
    
    redExon <- data.frame(geneR = unlist(lapply(strsplit(df.l$gene, split = "[.]"), "[[", 1)),
                          chr = sapply(strsplit(df.l$exon, split = ":"), "[[", 1),
                          start = sapply(strsplit(sapply(strsplit(df.l$exon, split = ":"), "[[", 2), split = "[-]"), "[[", 1),
                          stop = sapply(strsplit(sapply(strsplit(df.l$exon, split = ":"), "[[", 2), split = "[-]"), "[[", 2)
    )
    
  }
  
  
  colnames(redExon) <- c("geneR", "chr", "start", "stop")
  redExon$start <- as.numeric(redExon$start)
  redExon$stop <- as.numeric(redExon$stop)
  print("exon loaded...")
  
  
  matched <- getTranscript(gtf = gtf, redExon = redExon, ex_type = ex_type, minOverlap = mOverlap, swaps = !(background), cores = inCores)
  print("exons matched, bed-ifying...")
  bed <- bedify(matched, saveBED=F, outname = outname, cores = inCores)
  
  
  trans <- unlist(lapply(strsplit(unique(bed$name), "#"), "[[", 1))
  possT <- unlist(lapply(strsplit(bed$name, "#"), "[[", 1))
  
  
  ## Find annotated proteins for transcripts if possible
  protCode <- unlist(mclapply(trans, mc.cores = 8, function(x) {
    rc <- c_trans[which(c_trans == x)+1]
    if (length(rc) > 0) {
      rc[1]
    } else {"none"}
  }))
  
  proBed <- data.frame(id = unique(bed$name), strand = unlist(lapply(unique(bed$name), function(x) unique(bed$strand[bed$name == x][1])[1])), prot = protCode) %>% separate(id, c("transcript", "id"), "#") %>% separate("id", c("gene", "chr"), ";") %>% separate('chr', c('chr', 'coords'), ':') %>% separate('coords', c('start', 'stop'), '-')
  
  proFast <- c()
  if (length(proBed[,1]) %% 2 == 0) {
    subVal <- 1
  } else {subVal <- 0}
  for (i in seq(1, (length(proBed[,1])-subVal), by = 2)) {
    proFast <- c(proFast, paste(">", proBed$transcript[i], "#", proBed$gene[i], ";", proBed$chr[i], ":", proBed$start[i], "-", proBed$stop[i], ";", proBed$strand[i], sep = ""), 
                 proBed$prot[i], paste(">", proBed$transcript[i+1], "#", proBed$gene[i+1], ";", proBed$chr[i+1], ":", proBed$start[i+1], "-", proBed$stop[i+1], ";", proBed$strand[i+1], sep = ""), 
                 proBed$prot[i+1])
  }
  
  if (background == T) {
    if (saveOutput == T) {
      write.table(proBed, paste("./", "proteinOut.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = '\t') 
      write_lines(proFast, paste("./", "outFast.fa", sep = ""))
    }
    return(list(matched = matched,
                bed = bed,
                proBed = proBed,
                proFast = proFast))
  } else {
    protAlign <- list()
    protC <- c()
    pMatch <- c()
    alignType <- c()
    cate <- c()
    for (i in seq(from=1,to=(length(protCode)-1), by=2)) {
      if (protCode[i] == "none" | protCode[i+1] == "none") {
        if (protCode[i] == "none" & protCode[i+1] != "none") {
          protC <- c(protC, "nonPC", "PC")
          protAlign[[i]] <- "none"
          protAlign[[i]] <- "onePC"
          alignType <- c(alignType, "onePC")
        } else if (protCode[i] != "none" & protCode[i+1] == "none") {
          protC <- c(protC, "PC", "nonPC")
          protAlign[[i]] <- "onePC"
          alignType <- c(alignType, "onePC")
        } else {
          protC <- c(protC, "nonPC", "nonPC")
          protAlign[[i]] <- "none"
          alignType <- c(alignType, "noPC")
        }
        pMatch <- c(pMatch, 0)
      } else if (protCode[i] == protCode[i+1]) {
        protC <- c(protC, "Same", "Same")
        protAlign[[i]] <- msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])))
        pMatch <- c(pMatch, 1.04)
        alignType <- c(alignType, "Match")
      } else {
        protC <- c(protC, "Different", "Different")
        protAlign[[i]] <- msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE)
        
        minPc <- min(nchar(protCode[i]), nchar(protCode[i+1]))
        pMatch <- c(pMatch, table(unlist(lapply(strsplit(msaConsensusSequence(protAlign[[i]]), split = ""), function(x) x == "?")))[1]/min(nchar(protCode[i]), nchar(protCode[i+1])))
        if (nchar(paste(strsplit(msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]][nchar(strsplit(msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]]) > (.1*minPc)], collapse = "")) > .2*minPc)  {
          alignType <- c(alignType, "PartialMatch")
          # msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE
          # , file = paste(out_dir, "prettyAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_pm_prettyAlignment.pdf", sep = ""), output = "pdf")
        } else {
          alignType <- c(alignType, "FrameShift")
          # msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE
          # , file = paste(out_dir, "prettyAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_fs_prettyAlignment.pdf", sep = ""), output = "pdf")
        }
      }
    }
    
    print(table(protC))
    proBed$match <- protC
    proBed$prop <- rep(pMatch, each = 2)
    
    # Filled Density Plot
    (gdf <- ggplot(data.frame(dens = as.numeric(pMatch), type = alignType), aes(x = dens, fill = type)) + 
        geom_histogram(aes(y=..count../sum(after_stat(count))), colour = 1,
                       bins = 20) + geom_density(aes(y=.0005*after_stat(count)), color = 'black', fill = "coral2", bw = .1, alpha = .3) + 
        scale_fill_manual(values=c('noPC' = "azure4", 'Match' = "#E69F00", 'onePC' = "#56B4E9", 'FrameShift' = "pink", 'PartialMatch' = "deeppink4")) +
        theme_classic() + xlab("Alignment Score") + ylab("Fraction"))
    
    
    proBed$matchType <- rep(alignType, each = 2)
    
    if (saveOutput == T) {
      write.table(proBed, paste("./", "proteinOut.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = '\t') 
      write_lines(proFast, paste("./", "outFast.fa", sep = ""))
      
      pdf(file = paste(out_dir, "alignScores.pdf", sep = ""))
      print(gdf)
      dev.off()
    }
    return(list(matched = matched,
                bed = bed,
                proBed = proBed,
                proFast = proFast,
                df.l = df.l,
                gdf = gdf,
                deExons = lfcPlot))
  }
  
}


lfc <- function(de_df, numCont, numExp, exon_type, cores = 8) {
  de_df <- de_df[de_df$type == exon_type,]
  samps <- colnames(de_df)
  # unlist(lapply(strsplit(colnames(de_df)[-c(1:10, (10+numCont+numExp+1):dim(de_df)[2])], split = "[.]"), "[[", 5))
  
  lfc <- list()
  col <- list()
  
  
  lfc <- mclapply(1:length(de_df$gene), mc.cores = cores, function(i) {
    outlier <- which(unlist(lapply(samps, function(x) grepl(x, de_df$outlier[i]))))+10
    
    cont <- c(11:(11-1+numCont))
    cont <- cont[!(cont %in% outlier)]
    exp <- c((numCont+10+1):(11-1+numCont+numExp))
    exp <- exp[!(exp %in% outlier)]
    
    
    
    if (length(exp) == 0) {
      lfc_t <- log2((.01)/as.numeric(rowMeans(de_df[i,cont, drop = FALSE])+.01))
    } 
    if (length(cont) == 0) {
      lfc_t <- log2((as.numeric(rowMeans(de_df[i,exp, drop = FALSE]))+.01)/.01)
    }
    if (length(exp) == 0 & length(cont) == 0) {
      lfc_t <- 0
    }
    if (length(exp) != 0 & length(cont) != 0) {
      lfc_t <- log((as.numeric(rowMeans(de_df[i,exp, drop = FALSE]))+.01)/as.numeric(rowMeans(de_df[i,cont, drop = FALSE])+.01))
    }
    
    
    lfc_t
  })
  
  de_df$lfc <- unlist(lfc)
  de_df <- de_df[de_df$p_value >= 0,]
  col <- list()
  for (i in 1:length(de_df$gene)) {
    col[[i]] <- "#A7A9AC"#natparks.pals("Acadia", 15)[7]
    if (de_df$lfc[i] <= -1.0 & de_df$p_value[i] < .01) {
      col[[i]] <- "#FE9234" #natparks.pals("Acadia", 15)[15]
    }
    if (de_df$lfc[i] >= 1.0 & de_df$p_value[i] < .01) {
      col[[i]] <- "#00A79D"# natparks.pals("Acadia", 15)[1]
    }
    
  }
  de_df$col <- unlist(col)
  
  p.scDE <- de_df
  p.scDE$numOutliers <- str_count(p.scDE$outlier, "[.]")
  return (de = p.scDE)
}


make_lfcPlot <- function(lfc_df, num_thresh = 10) {
  hg38.conv <- readRDS("/projectnb2/evolution/zwakefield/proteinChange/pipeline/hg38_geneRef_conv.RDS")
  lfc_df$hgnc <- unlist(lapply(lfc_df$gene, function(x) {ifelse(unlist(lapply(strsplit(x, split = "[.]"), "[[", 1)) %in% hg38.conv$ens, 
                                                                unique(hg38.conv$hgnc[hg38.conv$ens == unlist(lapply(strsplit(x, split = "[.]"), "[[", 1))]), x)
    
  }))
  hg38.conv$hgnc[hg38.conv$ens %in% unlist(lapply(strsplit(lfc_df$gene, split = "[.]"), "[[", 1))]
  
  
  lab_thresh <- lfc_df %>% arrange(desc(abs(lfc)), p_value)
  lfc_df$hgnc[-log(lfc_df$p_value) <= -log(lab_thresh$p_value[num_thresh]) | abs(lfc_df$lfc) < abs(lab_thresh$lfc[num_thresh])]  <- ""
  
  (deExons <- ggplot(lfc_df, aes(x = lfc, y = -log(p_value), color = col, label = hgnc)) + geom_point(aes(shape = type), size = 2, color = lfc_df$col) +
      theme_classic() + ylab("-Log Adj P Value") + xlab("Log2FC")
    # +coord_cartesian(xlim = c(-100, 20))
    + geom_text(hjust=.2, vjust=0, size = 4)
  )
  print("lfc done!")
  return(deExons)
}