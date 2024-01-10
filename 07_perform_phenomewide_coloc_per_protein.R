#!/usr/bin/env Rscript

## script to run coloc for Olink regions
## Maik Pietzner 19/01/2023
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)

setwd("path")

## --> packages required <-- ##

require(data.table)
require(susieR)
require(coloc)
require(doMC)
require(Rfast)
require(ieugwasr)

## --> import parameters <-- ##

olink <- args[1]
id    <- args[2]
chr.s <- args[3]
pos.s <- as.numeric(args[4])
pos.e <- as.numeric(args[5])


cat("Run coloc with", olink, id, chr.s, pos.s, pos.e, "\n")

#-----------------------------------------#
##--     load regional assoc stats     --##
#-----------------------------------------#

cat("------------------------------\n")
cat("Reading summary statistics in \n")

## import regional statistics with credible set assignment 
res.olink <- fread(paste("path to fine-mapping results",
                         olink, chr.s, pos.s, pos.e, "txt", sep="."))

cat("Found", nrow(res.olink), "entries \n")
cat("------------------------------\n")

#-----------------------------------------#
##--       import genotype data        --##
#-----------------------------------------#

cat("------------------------------\n")
cat("Reading snp data in\n")

## TODO: replace sample list with the ones actually used

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
res.olink[, Effect := ifelse(EA == alleleB, Effect, -Effect)]
## drop what is no longer needed
res.olink        <- res.olink[, c("id", "snp.id", "rsid", "MarkerName", "chr", "pos", "alleleA", "alleleB", "Effect", "StdErr", "Pvalue",
                                  "TotalSampleSize", "MAF", "pip", "cs")]
## rename
names(res.olink) <- c("id", "snp.id", "rsid", "MarkerName", "chr", "pos", "NEA", "EA", "Effect", "StdErr", "Pvalue",
                      "TotalSampleSize", "MAF", "pip", "cs")

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
##--     obtain top SNPs and proxies    --##
#------------------------------------------#

## align order
res.olink  <- as.data.table(res.olink)
## order, keep in mind that cs -1 means not included in the credible set
res.olink  <- res.olink[ order(-cs, -pip) ]
## create indicator to select
res.olink[, ind := 1:.N, by="cs"]
## get top SNPs
top.snp    <- res.olink[ ind == 1 & cs > 0]
## add LD columns: be careful 'cs' does not mean that numbers match, but is a legacy from fine-map
for(j in top.snp$cs){
  res.olink[, paste0("R2.", j)] <- ld[ res.olink$id, top.snp$id[which(top.snp$cs == j)]]^2
}

## get all proxy SNPs
proxy.snps <- paste0("R2.", top.snp$cs)
proxy.snps <- apply(res.olink[, ..proxy.snps], 1, function(x){
  ifelse(sum(x >= .8) > 0, T, F)
})
## may include SNPs not in the credible set
proxy.snps <- res.olink[proxy.snps]

## store the LD pattern across top SNPs (convert to data frame to ease downstream operations)
ld.top.snps <- as.data.frame(res.olink)
ld.top.snps <- lapply(top.snp$cs, function(x){
  ## get all SNPs and corresponding LD
  tmp        <- paste0("R2.", x)
  tmp        <- ld.top.snps[which(ld.top.snps[, tmp] >= .8), c("MarkerName", "id", "rsid", tmp)]
  ## edit names
  names(tmp) <- c("MarkerName.proxy", "id.proxy", "rsid.proxy", "R2")
  print(tmp)
  ## add top SNP
  tmp        <- merge(as.data.frame(top.snp[x, c("MarkerName", "id", "rsid", "cs")]), tmp, suffix=c(".lead", ".proxy"))
  ## do some renaming to ease downstream coding
  names(tmp) <- c("MarkerName.lead", "id.lead", "rsid.lead", "cs", "MarkerName.proxy", "id.proxy", "rsid.proxy", "R2.proxy")
  ## return
  return(tmp)
})
## combine everything
ld.top.snps <- do.call(rbind, ld.top.snps)

#-----------------------------------------#
##--  perform PheWAS for all proxies   --##
#-----------------------------------------#

## query everything (omit FinnGen and BBJ as sample sizes are missing)
tmp.phewas <- NULL
while(is.null(dim(tmp.phewas))){
  tmp.phewas <- phewas(variants=ld.top.snps$rsid.proxy, pval=1e-6, batch=c("ebi-a", "ieu-a", "ieu-b", "met-a", "ubm-a", "ukb-b", "finn-b"))
}

## define break criteria here as well
if(nrow(tmp.phewas)>0){
  
  ## add respective lead SNPs from the credible sets
  tmp.phewas <- as.data.table(tmp.phewas)
  tmp.phewas <- merge(ld.top.snps, tmp.phewas, by.x=c("rsid.proxy"), by.y="rsid", suffixes = c(".snp", ".ieu"))
  ## careful: rsid now tags the lead variant of the respective credible set and rsid.proxy is the variant queried!!
  
  ## convert to data table
  tmp.phewas <- as.data.table(tmp.phewas)
  
  ## take the strongest association for each snp - datatype
  tmp            <- tmp.phewas[order(rsid.lead, trait, -abs(beta/se)),]
  tmp[, ind:= 1:.N, by=c("rsid.lead", "trait")]
  tmp            <- tmp[ind == 1]
  tmp$ind        <- NULL
  
  ## add protein stats 
  phewas.results <- merge(tmp, res.olink, by.x=c("MarkerName.proxy", "rsid.proxy"), by.y=c("MarkerName", "rsid") ,  suffix=c(".trait", ".pQTL"))
  
  #-----------------------------------------#
  ##--  run naive coloc for all results  --##
  #-----------------------------------------#
  
  res.naive <- lapply(1:nrow(phewas.results), function(x){
    
    cat("\n-------------------------------\n")
    
    print(phewas.results$trait[x])
    
    ## obtain summary stats and trait info
    reg             <- paste0(chr.s, ":", pos.s,"-", pos.e)
    ## ensure that the dataset is queried
    res.trait       <- NULL
    while(is.null(dim(res.trait))){
      cat("\n obtain GWAS stats \n")
      res.trait       <- associations(reg, phewas.results$id.trait[x])
    }
    ## just to make sure everything will work fine
    res.trait       <- subset(res.trait, !is.na(beta))
    ## make unique, since some data processing error
    res.trait       <- unique(res.trait)
    ## replace X
    res.trait$chr   <- as.numeric(gsub("X", 23, res.trait$chr))
    res.trait$se    <- abs(res.trait$se)
    ## drop zero SE's
    res.trait       <- subset(res.trait, se > 0)
    
    # print(head(res.trait))
    
    ## get meta information (careful some might miss those)
    tr.info         <- tryCatch(
      {
        gwasinfo(phewas.results$id.trait[x])
      }, error=function(e){
        return(NA)
      })
    
    ## merge (this might ignore INDELS)
    res.all                   <- merge(res.olink, res.trait, by=c("rsid", "chr"), suffixes = c(".pQTL", ".trait"))
    ## remove non-biallelelic variants
    ii                        <- table(res.all$rsid)
    res.all                   <- subset(res.all, rsid %in% names(ii[ii==1]))
    res.all                   <- as.data.frame(res.all)
    
    ## account for possible INDELs
    res.all[, c("ea", "nea")] <- t(apply(res.all[, c("ea", "nea")], 1, function(x){
      
      if(nchar(x[1]) > 1 | nchar(x[2]) > 1){
        ## replace
        if(nchar(x[1]) > nchar(x[2])){
          return(c("I", "D"))
        }else{
          return(c("D", "I"))
        }
      }else{
        return(x)
      }
      
    }))
    
    ## align effect estimates
    res.all$beta.aligned <- ifelse(res.all$EA == res.all$ea, res.all$beta, -res.all$beta)
    ## delete what is no longer needed
    res.all$ea           <- res.all$nea <- NULL
    
    #-----------------------------------------#
    ##-- 	         sanity check            --##
    #-----------------------------------------#
    
    ## top signal for the protein in all data
    ts      <- res.olink$id[which.max(abs(res.olink$Effect/res.olink$StdErr))]
    ## top signal for protein in the overlap
    is      <- res.all$id.pQTL[which.max(abs(res.all$Effect/res.all$StdErr))]           
    ## get the top SNP for the outcome
    io      <- res.all$id.pQTL[which.max(abs(res.all$beta/res.all$se))]
    
    ## conserved signal for phecode
    ld.sens <- ld[ts, is]^2
    ## ld between lead signals
    ld.ovl  <- ld[is, io]^2
    
    #-----------------------------------------#
    ##-- 	            run coloc            --##
    #-----------------------------------------#
    
    ## order by position
    res.all     <- as.data.table(res.all)
    res.all     <- res.all[order(pos)]
    
    ## prepare input
    D1          <- list(beta=res.all$Effect, varbeta=res.all$StdErr^2, 
                        type="quant", 
                        sdY=1,
                        N=max(res.all$TotalSampleSize),
                        MAF=res.all$MAF,
                        snp=res.all$id.pQTL,
                        position=1:nrow(res.all))
    
    ## try out
    if(!("ncase" %in% names(tr.info))){
      print("use quantitative")
      
      ## define N
      if(!is.null(dim(tr.info))){
        n <- tr.info$sample_size
      }else{
        n <- max(res.all$n, na.rm=T)
        ## adopt for FinnGen
        if(n == ""){
          n <- NULL
        }
      }
      
      D2          <- list(beta=res.all$beta.aligned, varbeta=res.all$se^2, 
                          type="quant", 
                          N=n,
                          sdY=1,
                          MAF=res.all$MAF,
                          snp=res.all$id.pQTL,
                          position=1:nrow(res.all))
    }else{
      ## binary outcome
      D2          <- list(beta=res.all$beta.aligned, varbeta=res.all$se^2, 
                          type="cc",
                          s=tr.info$ncase/(tr.info$ncontrol+tr.info$ncase), 
                          N=tr.info$sample_size,
                          MAF=res.all$MAF,
                          snp=res.all$id.pQTL,
                          position=1:nrow(res.all))
    }
    
    ## do naive coloc as well
    naive.coloc                    <- coloc.signals(D1, D2, method="single", p12=5e-6)
    
    ## add checks to the data
    naive.coloc$summary$ld.sens    <- ld.sens
    naive.coloc$summary$ld.overlap <- ld.ovl
    
    ## add the trait id and label
    naive.coloc$summary$id.ieu            <- phewas.results$id.trait[x]
    naive.coloc$summary$phenotype.ieu     <- phewas.results$trait[x]
    naive.coloc$summary$id.olink          <- olink
    
    #-----------------------------------------#
    ##-- 	        draw selected            --##
    #-----------------------------------------#
    
    if(naive.coloc$summary$PP.H4.abf > .7 | naive.coloc$summary$ld.overlap > .8){
      source("scripts/plot_locus_compare.R")
      png(paste0("graphics/", olink, ".", phewas.results$id.trait[x], ".", chr.s, ".", pos.s, ".", pos.e, ".png"), width=16, height=8, units="cm", res=300)
      par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=.01, cex.main=.6, font.main=2)
      ## more complex layout for gene assignment
      layout(matrix(c(1,1,1,2,3,4),3,2), heights = c(.43,.37,.2))
      plot.locus.compare(res.all, naive.coloc$summary, ld, a.vars=phewas.results$rsid.proxy[x])
      dev.off()
    }
    
    
    ## write results to file
    return(naive.coloc$summary)
    
  })
  res.naive <- do.call(rbind, res.naive)
  
  #-----------------------------------------#
  ##--  run susie coloc for all results  --##
  #-----------------------------------------#
  
  ## order just to make sure
  res.olink     <- res.olink[order(pos)]
  
  ## run coloc
  res.susie    <- lapply(1:nrow(phewas.results), function(x){
    
    cat("\n-------------------------------\n")
    
    print(phewas.results$id.trait[x])
    
    ## obtain summary stats and trait info
    reg             <- paste0(chr.s, ":", pos.s,"-", pos.e)
    ## ensure that the dataset is queried
    res.trait       <- NULL
    while(is.null(dim(res.trait))){
      res.trait       <- associations(reg, phewas.results$id.trait[x])
    }
    ## just to make sure everything will work fine
    res.trait       <- subset(res.trait, !is.na(beta))
    ## make unique, since some data processing error
    res.trait       <- unique(res.trait)
    ## replace X
    res.trait$chr   <- as.numeric(gsub("X", 23, res.trait$chr))
    
    # print(head(res.trait))
    
    ## get meta information (careful some might miss those)
    tr.info         <- tryCatch(
      {
        gwasinfo(phewas.results$id.trait[x])
      }, error=function(e){
        return(NA)
      })
    
    ## merge (this might ignore INDELS)
    res.all                   <- merge(res.olink, res.trait, by=c("rsid", "chr"), suffixes = c(".pQTL", ".trait"))
    ## remove non-biallelelic variants
    ii                        <- table(res.all$rsid)
    res.all                   <- subset(res.all, rsid %in% names(ii[ii==1]))
    res.all                   <- as.data.frame(res.all)
    
    ## account for possible INDELs
    res.all[, c("ea", "nea")] <- t(apply(res.all[, c("ea", "nea")], 1, function(x){
      
      if(nchar(x[1]) > 1 | nchar(x[2]) > 1){
        ## replace
        if(nchar(x[1]) > nchar(x[2])){
          return(c("I", "D"))
        }else{
          return(c("D", "I"))
        }
      }else{
        return(x)
      }
      
    }))
    
    ## align effect estimates
    res.all$beta.aligned <- ifelse(res.all$EA == res.all$ea, res.all$beta, -res.all$beta)
    ## delete what is no longer needed
    res.all$ea           <- res.all$nea <- NULL  
    
    ## compute MAF and exclude extreme outliers
    if(sum(!is.na(res.all$eaf)) > nrow(res.all) * .5){
      res.all$maf.ieu      <- ifelse(res.all$eaf > .5, 1-res.all$eaf, res.all$eaf)
      ## subset
      res.all              <- subset(res.all, abs(MAF - maf.ieu) <= .1)
    }else{
      res.all$maf.ieu      <- NA
    }
    
    #-----------------------------------------#
    ##-- 	        prepare coloc            --##
    #-----------------------------------------#
    
    ## order by position
    res.all     <- as.data.table(res.all)
    res.all     <- res.all[order(pos)]
    
    ## prepare input
    D1          <- list(beta=res.all$Effect, varbeta=res.all$StdErr^2, 
                        type="quant", 
                        sdY=1, 
                        N=max(res.all$TotalSampleSize, na.rm=T),
                        MAF=res.all$MAF,
                        snp=res.all$id.pQTL,
                        position=1:nrow(res.all),
                        ## subset the LD matrix to what is needed
                        LD=as.matrix(ld[res.all$id.pQTL, res.all$id.pQTL]))
    
    ## try out
    if(!("ncase" %in% names(tr.info))){
      print("use quantitative")
      
      ## define N
      if(!is.null(dim(tr.info))){
        n <- tr.info$sample_size
      }else{
        n <- max(res.all$n, na.rm=T)
        ## adopt for FinnGen
        if(n == ""){
          n <- NULL
        }
      }
      
      D2          <- list(beta=res.all$beta.aligned, varbeta=res.all$se^2, 
                          type="quant", 
                          N=n,
                          sdY=1,
                          MAF=res.all$MAF,
                          snp=res.all$id.pQTL,
                          position=1:nrow(res.all),
                          ## subset the LD matrix to what is needed
                          LD=as.matrix(ld[res.all$id.pQTL, res.all$id.pQTL]))
    }else{
      ## binary outcome
      D2          <- list(beta=res.all$beta.aligned, varbeta=res.all$se^2, 
                          type="cc",
                          s=tr.info$ncase/(tr.info$ncontrol+tr.info$ncase), 
                          N=tr.info$sample_size,
                          MAF=res.all$MAF,
                          snp=res.all$id.pQTL,
                          position=1:nrow(res.all),
                          ## subset the LD matrix to what is needed
                          LD=as.matrix(ld[res.all$id.pQTL, res.all$id.pQTL]))
    }
    
    ## Olink protein
    set.seed(42)
    susie.olink <- tryCatch(
      {
        runsusie(D1, 
                 # r2.prune = .25, 
                 max_iter = 10000,
                 L=ifelse(nrow(top.snp) == 1, 2, nrow(top.snp)))
      }, error=function(e){
        return(NA)
      })
    
    ## trait
    susie.trait <- tryCatch(
      {
        runsusie(D2, 
                 # r2.prune = .25, 
                 max_iter = 10000,
                 L=5)
      }, error=function(e){
        return(NA)
      })
    
    ## check whether both have outcome
    if(length(susie.olink) > 1 & length(susie.trait) > 1){
      
      ## additional check, whether both traits have at least one credible set to be tested
      if(length(summary(susie.olink)$cs) > 0 & length(summary(susie.trait)$cs) > 0){
        
        #-----------------------------------------#
        ##-- 	            run coloc            --##
        #-----------------------------------------#
        
        ## run coloc with susie input
        res.coloc        <- coloc.susie(susie.olink, susie.trait, p12 = 5e-6)
        
        #-----------------------------------------#
        ##--   cross-check with fine-mapping   --##
        #-----------------------------------------#
        
        ## add LD with fine-mapped variants 
        res.coloc        <- res.coloc$summary
        
        print(res.coloc)
        
        ## compute ld between selected lead hits
        res.coloc$ld.top <- apply(res.coloc[, c("hit1", "hit2"), drop=F], 1, function(k) ld[k[1], k[2]]^2)
        
        #-----------------------------------------#
        ##--     add additional information    --##
        #-----------------------------------------#
        
        ## add effect estimates (restrict to selected pQTLs, which means no effect estimates for the regional lead for the trait are taken forward)
        res.coloc                   <- merge(res.coloc, res.all, by.x = "hit1", by.y = "id.pQTL")
        
        ## add phecode
        res.coloc$id.olink          <- olink
        ## add trait
        res.coloc$id.ieu            <- phewas.results$id.trait[x]
        res.coloc$phenotype.ieu     <- phewas.results$trait[x]
        
        ## add rsids and LD with other variants
        res.coloc                   <- merge(res.coloc, res.olink[, c("id", "rsid")], by.x="hit2", by.y = "id", suffixes = c(".olink", ".ieu"))
        
        #-----------------------------------------#
        ##-- 	        draw selected            --##
        #-----------------------------------------#
        
        if(sum(res.coloc$PP.H4.abf > .7) > 0 | sum(res.coloc$ld.top > .8) > 0){

          source("scripts/plot_locus_compare.R")
          png(paste("graphics/susie", olink, phewas.results$id.trait[x], chr.s, pos.s, pos.e, "png", sep="."), width=16, height=8, units="cm", res=200)
          par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=.01, cex.main=.6, font.main=2)
          ## more complex layout for gene assignment
          layout(matrix(c(1,1,1,2,3,4),3,2), heights = c(.43,.37,.2))
          plot.locus.compare(res.all, res.coloc, ld, a.vars=unique(c(res.coloc$rsid.olink, res.coloc$rsid.ieu)))
          dev.off()
        }
        
        ## write results to file
        return(res.coloc)
        
      }else{
        return(NULL)
      }
      
    }else{
      return(NULL)
    }
    
  })
  res.susie        <- do.call(rbind, res.susie)
  
  #-----------------------------------------#
  ##--  	combine and store the data     --##
  #-----------------------------------------#
  
  ## combine both (N.B. naive results may include duplications if multiple independent variants are associated with the phenotype as well)
  phewas.results                   <- merge(unique(res.naive), phewas.results, by.x=c("id.ieu", "phenotype.ieu"), by.y=c("id.trait", "trait"))
  
  ## account for possible INDELs
  phewas.results                   <- as.data.frame(phewas.results)
  phewas.results[, c("ea", "nea")] <- t(apply(phewas.results[, c("ea", "nea")], 1, function(x){
    
    if(nchar(x[1]) > 1 | nchar(x[2]) > 1){
      ## replace
      if(nchar(x[1]) > nchar(x[2])){
        return(c("I", "D"))
      }else{
        return(c("D", "I"))
      }
    }else{
      return(x)
    }
    
  }))
  
  ## create aligned results
  phewas.results        <- as.data.table(phewas.results)
  phewas.results[, beta.trait := ifelse(EA == ea, beta, -beta)]
  ## keep only what is needed
  phewas.results        <- phewas.results[, c("id.olink", "id.ieu", "phenotype.ieu", "ld.sens", "ld.overlap", "nsnps", paste0("PP.H", 0:4, ".abf"),
                                              "MarkerName.lead", "rsid.lead", "MarkerName.proxy", "rsid.proxy", "EA", "NEA", "MAF", "Effect", "StdErr", "Pvalue", "beta.trait",
                                              "se", "p")]
  ## edit names
  names(phewas.results) <- c("id.olink", "id.ieu", "phenotype.ieu", "ld.sens", "ld.overlap", "nsnps", paste0("PP.H", 0:4, ".abf"),
                             "MarkerName.lead", "rsid.lead", "MarkerName.proxy", "rsid.proxy", "EA", "NEA", "MAF", "beta.pqtl", "se.pqtl", "p.pqtl", "beta.ieu",
                             "se.ieu", "p.ieu")
  
  ## write to file
  write.table(phewas.results, paste("output/results.pqtl.phewas.naive", olink, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names=F)
  
  ## keep only what is really needed
  if(!is.null(dim(res.susie))){
    res.susie        <- res.susie[, c("id.olink", "id.ieu", "phenotype.ieu", "idx1", "idx2", "ld.top", "nsnps",
                                      paste0("PP.H", 0:4, ".abf"), "MarkerName", "rsid.olink", "chr", "pos", "cs", "pip", "EA", "NEA", "MAF", "Effect", "StdErr", "Pvalue", 
                                      "beta.aligned", "se", "p", grep("R2", names(res.susie), value=T)), with=F]
    
    ## rename
    names(res.susie) <- c("id.olink", "id.ieu", "phenotype.ieu", "cs.set.phecode", "cs.set.ieu", "ld.top", "nsnps",
                          paste0("PP.H", 0:4, ".abf"), "MarkerName", "rsid", "chr", "pos", "cs", "pip", "EA", "NEA", "MAF", "beta.pqtl", "se.pqtl", "p.pqtl",
                          "beta.ieu", "se.ieu", "p.ieu", grep("R2", names(res.susie), value=T))
    
    ## write to file
    write.table(res.susie, paste("output/results.pqtl.phewas.susie", olink, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names=F)
  }else{
    
    cat("found no evidence for", olink, "on chr", chr.s, "between", pos.s, "and", pos.e, "\n")
    ## write to file
    write.table(data.frame(pheno=olink, id.ieu=NA), paste("output/results.pqtl.phewas.susie", olink, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names=F)
    
  }
    
}else{
  
  cat("found no evidence for", olink, "on chr", chr.s, "between", pos.s, "and", pos.e, "\n")
  ## write to file
  write.table(data.frame(pheno=olink, id.ieu=NA), paste("output/results.pqtl.phewas.susie", olink, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names=F)
  write.table(data.frame(pheno=olink, id.ieu=NA), paste("output/results.pqtl.phewas.naive", olink, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names=F)
  
  
}
