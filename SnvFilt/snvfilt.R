#!/usr/bin/R

library(data.table)

#read 
dt <- fread("~/project/Ecoli/input_test_fulltable.txt", sep="\t") #Chrom   Pos     Ref     E3_alt  E3_altnum       E3_cover        E4_alt  E4_altnum E4_alt ...
error_rate <- fread("~/project/Ecoli/error_rate.txt", sep="\t") #two columns: RefVar ErrorRate. Ex: AC      0.0003546

#save error rate to a vector
error_rate_vect <- error_rate$V2
names(error_rate_vect) <- error_rate$V1

#Change "E.coli" to "Ecoli"
dt[,Chrom:="Ecoli"]

#sample name in table, use these name to get information from full table
sampname <- c(paste0("E",3:14), "E3Merge", "E14Merge", "E3C")

#calculate pvalue and write out for each sample name
for (samp in sampname){
  cols <- c("Chrom", "Pos", "Ref", paste0(samp, c("_alt", "_cover", "_altnum")))
  dt_samp <- dt[,cols, with=F]
  colnames(dt_samp) <- c("Chrom", "position", "Ref", "Var", "Cov", "Alt")
  dt_samp <- dt_samp[Var %in% c("A", "T", "G", "C")]

  #calculate allele frequency
  dt_samp[, Freq:=Alt/Cov]

  #calculate p value according to bionominal model
  dt_samp[, Pvalue:= pbinom(Alt, Cov, error_rate_vect[paste0(Ref, Var)],0)]

  #filter: allele frequency <= 0.1, p value < 0.01
  result <- dt_samp[Freq <= 0.1 & Pvalue < 0.01]
 
  #change decimal digits and write out
  result$Freq <- paste0(round(result$Freq*100, 7), "%")
  result$Pvalue <- round(result$Pvalue, 7)
  fwrite(result, file=paste0("~/project/Ecoli/result/snv/", samp, ".snv"), sep="\t")
}

#clear memory
rm(dt, dt_samp, result, error_rate, error_rate_vect);gc()



