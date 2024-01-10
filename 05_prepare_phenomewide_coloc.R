################################################
#### Run cis-coloc with fine-mapped stats   ####
#### Maik Pietzner               19/01/2023 ####
################################################

rm(list=ls())
setwd("../path")
options(stringsAsFactors = F)
load(".RData")

## --> packages needed <-- ##
require(data.table)

###############################
#### import Olink targets  ####
###############################

## import Olink targets with at least one credible set
olink.targets <- fread("Olink.targets.with.cis.credible.set.20230119.txt")
## n = 1,314

## store input file for the pipeline
write.table(olink.targets[, c("pheno", "mrc_olink.id", "chr", "region_start", "region_end")], "Olink.coloc.targets.txt", sep="\t", col.names = F, row.names = F, quote = F)

###############################
####     import results    ####
###############################

#------------------------#
##--    SuSiE coloc   --##
#------------------------#

## get all output files
ii <- dir("../output/")
## subset to those with naive results
ii <- grep("susie", ii, value=T)

## collect everything
res.susie <- lapply(ii, function(x){
  
  ## import results
  tmp <- fread(paste0("../output/", x))
  ## add information to be able to map back to target information
  x   <- strsplit(x, "\\.")[[1]]
  ## add to the results
  tmp[, pheno := x[5]]
  tmp[, chr := as.numeric(x[6])]
  tmp[, region_start := as.numeric(x[7])]
  tmp[, region_end := as.numeric(x[8])]
  ## indicator whether successful
  tmp[, failed := F]
  ## return results
  return(tmp)
  
})
## combine everything again
res.susie <- do.call(plyr::rbind.fill, res.susie)

## combine with target information to see what is missing
res.susie <- merge(olink.targets, res.susie, by=c("pheno", "chr", "region_start", "region_end"), all.x = T)
res.susie <- as.data.table(res.susie)

# ## write to file what has been missed
# tmp       <- unique(rbind(res.naive[ is.na(failed), c("pheno", "mrc_olink.id", "chr", "region_start", "region_end")], res.susie[ is.na(failed), c("pheno", "mrc_olink.id", "chr", "region_start", "region_end")]))
# write.table(tmp, "Olink.coloc.targets.rerun.txt", sep="\t", col.names = F, row.names = F, quote = F)

###############################
####        apply QC       ####
###############################

## drop results w/o preserved pQTL
res.susie <- as.data.table(res.susie)
## related to any credible set
res.susie[, R2.max := pmax(R2.1, R2.2, R2.3, R2.3, R2.4, R2.5, R2.6, R2.7, R2.8, R2.9, R2.10, na.rm = T)]
res.susie <- res.susie[ R2.max >= .8 & nsnps >= 500 & ld.top >= .8 ]
## apply some pruning to outcome names
res.susie[, phenotype.ieu := gsub("\\\\", "", phenotype.ieu)]
res.susie[, phenotype.ieu := gsub("\"", "", phenotype.ieu)]
## drop findings not meeting stringent criteria
res.susie <- res.susie[ PP.H4.abf >= .8 ]
res.susie <- unique(res.susie)
## ensure marginal findings
res.susie <- res.susie[ as.numeric(p.pqtl) < 1e-5 & p.ieu < 1e-5 ]

## store both, but use susie for follow-up
write.table(res.susie, "Results.EPIC.Olink.cis.pQTL.phenome.wide.susie.coloc.20230315.txt", sep="\t", row.names = F)
