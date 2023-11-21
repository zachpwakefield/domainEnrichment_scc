#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(methods)
  library(ggpubr)
  library(xlsx)
  library(ggsignif)
  library(tidyverse)
  library(ComplexHeatmap)
  library(purrr)
  library(foreach)
  library(doParallel)
  library(pracma)
  library(hypeR)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)
  library(msa)
  library(stringr)
  library(ORFik)               
  library(GenomicFeatures)
  library(data.table)
  library(Biostrings) 
  library(parallel)
  library(readr)
})
args = commandArgs(trailingOnly=TRUE) # args <- '/projectnb2/evolution/zwakefield/proteinImpacts/mpc_cpc/'

### ### ### Domain Enrichment Analysis ### ### ###
cuD <- args[1]
tf <- list.files(cuD)
ipscan <- tf[grep("fa.tsv", tf)]
outFast <- tf[(!(grepl(".tsv", tf)) & !(grepl("xml", tf)) & !(grepl("gff3", tf)) & !(grepl("json", tf)) & grepl("outFast.fa", tf))]
outBed <- tf[grep("outBed.csv", tf)]


interproscan_results <- lapply(c("bg", "fg"), function(o) {
  ip <- read_delim(paste(cuD, ipscan[grep(o, ipscan)], sep = ""), col_names = F, delim = '\t')
  s <- read_csv(paste(cuD, outBed[grep(o, outBed)], sep =""))
  s$id <- paste(paste(paste(paste(s$transcript, s$gene, sep = "#"), s$chr, sep = ";"), paste(s$start, s$stop, sep = "-"), sep = ":"), s$strand, sep = ";")
  fa <- read_lines(paste(cuD, outFast[grep(o, outFast)], sep = ""))
  nFa <- fa[grep(">", fa)]
  gN <- gsub(">", "", nFa)
  geneL <- unlist(lapply(strsplit(unlist(lapply(strsplit(fa[grep(">", fa)], split = ";"), "[[", 1)), split = "#"), "[[", 2))
  ips <- ip %>% dplyr::filter(X4 %in% c("Pfam", "ProSitePatterns", "CDD", "PANTHER"))
  table(ips$X13)
  protInf <- list()
  upS <- seq(from=1,to=(length(gN)), by=1)
  for (j in upS) {
    protInf[[j]] <- ips$X13[ips$X1 == gN[j]]
    protInf[[j]] <- paste(protInf[[j]], collapse = ';')
  }
  protInf.o <- unlist(protInf)
  protInf.o[protInf.o == ""] <- "none"
  finOut <- data.frame(exon = gN,
                       protein = unlist(lapply(gN, function(x) s$prot[s$id == x])),
                       protInfor = protInf.o)
})

bg_dom <- unlist(strsplit(interproscan_results[[1]]$protInfor, split = ';'))
fg_dom <- unlist(strsplit(interproscan_results[[2]]$protInfor, split = ';'))
search <- unique(fg_dom)

successes <- lapply(search, function(x) c(sum(fg_dom == x), sum(bg_dom == x)))
pop_size <- length(bg_dom)
sample_size <- length(search)

pvals <- suppressWarnings(signif(stats::phyper(sapply(successes, "[[", 1)-1,
                                               sapply(successes, "[[", 2),
                                               pop_size-sapply(successes, "[[", 2),
                                               sample_size, lower.tail=FALSE), 3))


pvals.adj <- signif(stats::p.adjust(pvals, method="fdr"), 3)

data <- data.frame(domain = search,
                   pval = pvals,
                   fdr = pvals.adj,
                   sample_size = length(fg_dom),
                   sample_successes = sapply(successes, "[[", 1),
                   sample_prop = sapply(successes, "[[", 1) / length(fg_dom),
                   
                   pop_size = length(bg_dom),
                   pop_successes = sapply(successes, "[[", 2),
                   pop_prop = sapply(successes, "[[", 2) / length(bg_dom),
                   protein_contain = unlist(lapply(search, function(x) paste(interproscan_results[[2]]$exon[grep(x, interproscan_results[[2]]$protInfor)], collapse = "%")))
                   
)


data$rel_prop <- (data$sample_prop + .01) / (data$pop_prop + .01)
data <- data[data$domain != "none" & data$domain != "" & data$domain != "-",] %>% arrange(desc(rel_prop))
# data[data$fdr <= .05 & data$sample_successes > 10,] %>% arrange(desc(rel_prop))

write.csv(data, paste0(cuD, 'domainEnrichment.csv'))


