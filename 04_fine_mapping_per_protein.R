#!/usr/bin/env Rscript

## script to run SusieR for Olink regions
## Maik Pietzner 17/01/2023
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)

setwd("path")

## --> packages required <-- ##

require(data.table)
require(susieR)
require(doMC)
require(Rfast)

## --> import parameters <-- ##

olink <- args[1]
id    <- args[2]
chr.s <- args[3]
pos.s <- as.numeric(args[4])
pos.e <- as.numeric(args[5])

cat("Run Fine-Mapping using SuSiE with", olink, id, chr.s, pos.s, pos.e, "\n")

#-----------------------------------------#
##--     load regional assoc stats     --##
#-----------------------------------------#

cat("------------------------------\n")
cat("Reading summary statistics in \n")

## read the relevant data
res.olink        <- paste0("zcat path to stats",
                           gsub("invn_Xervv-1_30094", "invn_Xervv1_30094", olink),"/", gsub("invn_Xervv-1_30094", "invn_Xervv1_30094", olink),"_EPICTATAA_FEMA.txt.gz",
                           " | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e, 
                           " '{if(($2 == chr && $3 >= low && $3 <= upp) || NR==1) print $0}' -")
## import
res.olink        <- data.table::fread(cmd = res.olink, sep = "\t", header = T)

## drop variants present in only a subset (restrict to complete case analysis)
res.olink        <- subset(res.olink, TotalSampleSize == max(res.olink$TotalSampleSize))

## apply MAF filter?
res.olink$MAF    <- ifelse(res.olink$Freq1 > .5, 1 - res.olink$Freq1, res.olink$Freq1)
res.olink        <- subset(res.olink, MAF >= .001)

## create new MarkerName just to be sure
res.olink$snp.id <- apply(res.olink, 1, function(x){
  paste0("chr", as.numeric(x[2]), ":", as.numeric(x[3]), "_", paste(toupper(sort(x[4:5])), collapse = "_"))
})

cat("Found", nrow(res.olink), "entries \n")
cat("------------------------------\n")

#-----------------------------------------#
##--       import genotype data        --##
#-----------------------------------------#

cat("------------------------------\n")
cat("Reading snp data in\n")

## write file to obtain SNP dosages
write.table(res.olink$MarkerName, paste("tmpdir/tmp", olink, chr.s, pos.s, pos.e, "lst", sep="."), row.names = F, col.names = F, quote = F)

## import function to so
source("scripts/import_snp_dosages.R")
## import
ld        <- get.ld.info(olink, chr.s, pos.s, pos.e)
## ease downstream coding
snp.dat   <- ld[[1]]
snp.info  <- ld[[2]]
## delete and clean
rm(ld); gc(); gc()

## delete file no longer needed
system(paste("rm tmpdir/tmp", olink, chr.s, pos.s, pos.e, "*", sep="."))

## add rsid to Olink file, include alleles to align effect directions!! allele1 and alleleB are effect alleles;
## 'id' is the key to map to SNP dosages
res.olink        <- merge(res.olink, snp.info[, c("snp.id", "rsid", "id", "alleleA", "alleleB")])
## reassign effects
res.olink[, Effect := ifelse(toupper(Allele1) == alleleB, Effect, -Effect)]
## drop what is no longer needed
res.olink        <- res.olink[, c("id", "snp.id", "rsid", "MarkerName", "chr", "pos", "alleleA", "alleleB", "Effect", "StdErr", "Pvalue",
                           "TotalSampleSize", "MAF")]
## rename
names(res.olink) <- c("id", "snp.id", "rsid", "MarkerName", "chr", "pos", "NEA", "EA", "Effect", "StdErr", "Pvalue",
                      "TotalSampleSize", "MAF")

## compute LD-matrix as input for fine-mapping later on
ld               <- Rfast::cora(as.matrix(snp.dat[, res.olink$id]))

print(ld[1:5,1:5])

## identify missing entries
ii               <- apply(ld, 2, function(x) sum(is.na(x)))

## drop those with all missing
ld               <- ld[names(ii[ii < nrow(res.olink)]), names(ii[ii < nrow(res.olink)])]

## drop from results as well
res.olink        <- subset(res.olink, id %in% row.names(ld))

## SNPs of interest
snps             <- snp.info$id[snp.info$snp.id %in% res.olink$snp.id]

cat("Found", length(snps), "SNPs \n")
cat("------------------------------\n")

#------------------------------------------#
##--  run fine-mapping with asso stats  --##
#------------------------------------------#

## align order
res.olink <- as.data.table(res.olink)
res.olink <- res.olink[ order(pos) ]
## order ld matrix accordingly
ld        <- ld[ res.olink$id, res.olink$id]

## do in parallel
registerDoMC(6)

## run fine mapping to obtain 95%-credible sets (throws an error due to possible rounding errors in LD matrix)
set.seed(42)
## run through different values for L and catch possible errors
res.fine <- mclapply(2:10, function(x){
  
  ## run SuSiE
  tmp <- tryCatch(
    {
      susie_rss(res.olink$Effect/res.olink$StdErr, as.matrix(ld), L = x, coverage = .95, min_abs_corr=.1, max_iter = 10000)
    }, error=function(e){
      return(list(pip=rep(NA, nrow(res.olink)),
                  sets=list(cs=NA),
                  converged=F))
    })
  ## add L
  tmp$L <- x
  ## return
  return(tmp)

}, mc.cores=6)
## reduce to entries with at least one credible set
jj       <- unlist(lapply(res.fine, function(x) length(na.omit(x$sets$cs))))
## delete
res.fine <- res.fine[which(jj > 0)]

## proceed only if any
if(length(res.fine) > 0){
  
  #------------------------------------------#
  ##--        prune sets if needed        --##
  #------------------------------------------#
  
  ## do LD assessment
  ld.top <- mclapply(res.fine, function(x){
    
    ## add PIPs to ease selection of SNPs
    tmp.olink    <- summary(x)$vars
    ## subset to those in credible sets
    tmp.olink    <- as.data.table(subset(tmp.olink, cs > 0))
    ## only top SNPs
    tmp.olink    <- tmp.olink[order(cs, -variable_prob)]
    ## get only the top
    tmp.olink[, ind := 1:.N, by="cs"]
    tmp.olink    <- tmp.olink[ ind == 1]
    ## get the names
    tmp.olink[, id := sapply(tmp.olink$variable, function(x) rownames(ld)[x])]
    ## generate LD
    top.ld       <- ld[tmp.olink$id, tmp.olink$id, drop=F]^2
    ## identify possible sets in LD (r2>0.25); set diagonal to zero to ease downstream analysis
    diag(top.ld) <- 0
    top.ld       <- reshape2::melt(top.ld)
    ## subset to possible problematic candidates
    top.ld       <- subset(top.ld, value >= .25)
    ## return
    return(top.ld)

  }, mc.cores=6)
  
  ## find the maximum set of unrelated variants
  jj       <- unlist(lapply(ld.top, nrow))  
  ## subset accordingly
  if(sum( jj > 0) > 0){
    res.fine <- res.fine[-which(jj > 0)]
  }
  
  ## add to the results, if any
  if(length(res.fine) > 0){
    
    ## take only the last one
    res.fine  <- res.fine[[length(res.fine)]]
    
    ## add the information on pip and cs to the summary statistics
    tmp        <- summary(res.fine)$vars
    tmp$id     <- sapply(tmp$variable, function(x) rownames(ld)[x])
    names(tmp) <- c("variable", "pip", "cs", "id")
    res.olink  <- merge(res.olink, tmp[, c("id", "pip", "cs")], by="id")
    
    ## write all results to file
    write.table(res.olink, paste("output/fine.mapped", olink, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names = F)
    
    ## write top credible sets to separate file (N.B. variants not in any credible set are indicated by -1)
    res.olink  <- as.data.table(subset(res.olink, cs > 0))
    ## create indicator of strongest
    res.olink  <- res.olink[ order(cs, -pip)]
    res.olink[, ind := 1:.N, by="cs"]

    ## write all results to file
    write.table(res.olink[ind == 1], paste("lead_variants/cis.pqtl", olink, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names = F)
    
    }else{
    
    cat("\n-----------------------\n")  
    cat("Found no credible sets for", olink, id, chr.s, pos.s, pos.e, "\n")
    
    ## create empty file to be able to check for successful run jobs later on
    write.table(data.frame(id=NA, cs=NA), paste("lead_variants/cis.pqtl", olink, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names = F)
    
  }

  }else{

    cat("\n-----------------------\n")  
    cat("Found no credible sets for", olink, id, chr.s, pos.s, pos.e, "\n")
    
    ## create empty file to be able to check for successful run jobs later on
    write.table(data.frame(id=NA, cs=NA), paste("lead_variants/cis.pqtl", olink, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names = F)
    
}



