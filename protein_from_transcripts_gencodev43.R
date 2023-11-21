## Script to make protein from transcript file:

gcpc_trans <- read_lines(paste('/projectnb2/evolution/zwakefield/proteinChange/pipeline/gencode.v43.pc_translations.fa', sep = ""))
index <- grep(">", gcpc_trans)
endInd <- c(index[-1], (length(gcpc_trans)+1))
new_trans <- list()
for (i in 1:length(index)) {
  new_trans[[i]] <- c(gcpc_trans[index[i]])
  new_trans[[i]] <- c(new_trans[[i]], paste(gcpc_trans[(index[i]+1):(endInd[i]-1)], collapse = ""))
}
c_trans <- unlist(new_trans)
write_lines(c_trans, '/projectnb2/evolution/zwakefield/proteinImpacts/protein_code_from_gencodev43.txt')


c_trans[seq(1, length(c_trans), by = 2)] <- unlist(lapply(strsplit(c_trans[seq(1, length(c_trans), by = 2)], split = "[|]"), 
                                                          function(y) unlist(strsplit(y[2], split = "[.]"))[1]))


write_lines(c_trans, '/projectnb2/evolution/zwakefield/proteinImpacts/f_protein_code_from_gencodev43.txt')
# read_lines('/projectnb2/evolution/zwakefield/proteinImpacts/protein_code_from_v43.txt')
# read_lines('/projectnb2/evolution/zwakefield/proteinImpacts/f_protein_code_from_v43.txt')