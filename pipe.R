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
})

args = commandArgs(trailingOnly=TRUE)
# args <- c('/projectnb2/evolution/zwakefield/proteinImpacts/', 'pc', 'AFE', '3', '3', "/projectnb2/evolution/zwakefield/Finished_Projects/CardiacDifferentiation/HIT_Stat_Xingpei_Analysis/HIT_Stat_Analysis/mpc_cpc", '/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/4-1744_T35/hit/', '/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/5-1748_T35/hit/', '/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/6-1749_T35/hit/', '/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/7-1744_T21/hit/', '/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/8-1748_T21/hit/', '/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/9-1749_T21/hit/')
gtf_path <- '/projectnb2/evolution/zwakefield/proteinChange/gtf_withInfo.csv'
curr_dir <- args[1]
outname <- args[2] #outname <- 'pc'
ex_type <- args[3] # ex_type <- "AFE"
nC <- as.numeric(args[4])#nC <- 3
nE <- as.numeric(args[5])#nE <- 3
exon_type <- ex_type
out_directory <- args[6]
ms <-  args[7]
ts <- args[8:length(args)] # total list of exon files, ts <- c('/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/10L-1/hit/','/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/10L-2/hit/','/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/10L-3/hit/')
# list of deExon files, ms <- c('/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/10L-1/hit/')


source('/projectnb2/evolution/zwakefield/proteinImpacts/init.R')
c_trans <- read_lines('/projectnb2/evolution/zwakefield/proteinImpacts/f_protein_code_from_gencodev43.txt')
gtf <- read.csv(gtf_path, header = T)[,-c(1)]



print("Extracting transcripts from Background Set...")
bg <- proteinExtract_pipe(files_dir = ts, background = T, mOverlap = .5, saveOutput = F, inCores = 8)
print("Extracting transcripts from Foreground Set...")
fg <- proteinExtract_pipe(files_dir = ms, background = F, mOverlap = .5, saveOutput = F, inCores = 8, nC = nC, nE = nE, exon_type = ex_type)

write_csv(bg$proBed, paste0(out_directory, "bgoutBed.csv"))
write_lines(bg$proFast, paste0(out_directory, "bgoutFast.fa"))
write_csv(bg$matched$out_matched,  paste0(out_directory, "bgmatched.csv"))
write_csv(bg$bed,  paste0(out_directory, "bgbed.csv"))


write_csv(fg$proBed, paste0(out_directory, "fgoutBed.csv"))
write_lines(fg$proFast, paste0(out_directory, "fgoutFast.fa"))
write_csv(fg$matched$out_matched,  paste0(out_directory, "fgmatched.csv"))
write_csv(fg$bed,  paste0(out_directory, "fgbed.csv"))
write_csv(fg$df.l,  paste0(out_directory, "fglfc.csv"))

pdf(file = paste0(out_directory, "alignPlot.pdf"))
fg$gdf
dev.off()
pdf(file = paste0(out_directory, "volcano.pdf"))
print(fg$deExons)
dev.off()

### ### ### deExon Analysis ### ### ### 

## deExon analysis (volcano, GO enrichment, etc)

### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ###



### ### ### Swap analysis ### ### ###

## RBP analysis
## TF binding analysis

### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ###


### ### ### Domain Domain Interaction Analysis ### ### ###











