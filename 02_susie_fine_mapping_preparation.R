################################################
#### Run SuSIE fine-mapping for cis-regions ####
#### Maik Pietzner               17/01/2023 ####
################################################

rm(list=ls())
setwd("../path")
options(stringsAsFactors = F)
load(".RData")

## --> packages needed <-- ##
require(data.table)
require(readxl)

#################################
#### import Olink annotation ####
#################################

## import genomic annotations 
olink.targets       <- fread("~/rfs/rfs-epid-rfs-mpB3sSsgAn4/Studies/People/Maik/pGWAS_Olink_EPIC/Genomic.positions.Olink.Explore.3k.20220404.txt")

## add GWAS id
olink.targets$pheno <- paste0("invn_X", tolower(olink.targets$mrc_olink.id))

## check chromosome definition
table(olink.targets$chromosome_name)

#-------------------------------#
##-- unfold protein targets  --##
##-- with multiple genes     --##
#-------------------------------#

## make unique protein - gene combinations
tmp <- lapply(1:nrow(olink.targets), function(x){
  
  ## test whether there are more than one gene
  ii <- grep("\\|", olink.targets$chromosome_name[x])
  
  
  if(length(ii) > 0){
    
    print(olink.targets[x,])
    
    ## restore coordinates
    ensembl_gene_id  <- strsplit(olink.targets$ensembl_gene_id[x], "\\|")[[1]]
    uniprotswissprot <- strsplit(olink.targets$uniprotswissprot[x], "\\|")[[1]]
    hgnc_symbol      <- strsplit(olink.targets$hgnc_symbol[x], "\\|")[[1]]
    chromosome_name  <- strsplit(olink.targets$chromosome_name[x], "\\|")[[1]]
    start_position   <- strsplit(olink.targets$start_position[x], "\\|")[[1]]
    end_position     <- strsplit(olink.targets$end_position[x], "\\|")[[1]]
    
    print(ensembl_gene_id)
    
    ## return
    return(data.frame(olink.targets[x, c("OlinkID", "UniProt", "Assay", "mrc_olink.id", "Panel", "strand",
                                         "call.rate.belowLOD.samples", "Set", "source", "sl.seqid.v4", "olink.panel.id", "pheno"), drop=F],
                      ensembl_gene_id=ensembl_gene_id, uniprotswissprot=uniprotswissprot, hgnc_symbol=hgnc_symbol,
                      chromosome_name=chromosome_name, start_position=start_position, end_position=end_position))
  }else{
    return(olink.targets[x,])
  }
  
})
## combine
tmp           <- do.call(rbind, tmp)
## n=2961 in total to test

## replace
olink.targets <- tmp
rm(tmp); gc()

## define regions
olink.targets$region_start <- as.numeric(pmax(0,as.numeric(olink.targets$start_position) - 5e5))
olink.targets$region_end   <- as.numeric(olink.targets$end_position) + 5e5

## make chromosome numeric
olink.targets$chr          <- as.numeric(gsub("X", 23, olink.targets$chromosome_name))

## drop two, that didn't map to coordinates
olink.targets              <- subset(olink.targets, !is.na(region_start))

##############################################
####   export input for fine-map script   ####
##############################################

## do some renaming
olink.targets$pheno <- gsub("\\.", "", olink.targets$pheno)

## drop X-chromosome (no GWAS has been run)
write.table(subset(olink.targets, chr != 23)[, c("pheno", "mrc_olink.id", "chr", "region_start", "region_end")], "Olink.fine.map.targets.txt", sep="\t", col.names = F, row.names = F, quote = F)

## write lables to file
write.table(olink.targets, "Labels.Olink.Explore.expanded.20230315.txt", sep="\t", row.names=F)

##############################################
####            import results            ####
##############################################

## get all lead variants
ii       <- dir("../lead_variants/")

## import
res.pqtl <- lapply(ii, function(x){
  
  ## import results
  tmp <- fread(paste0("../lead_variants/", x))
  ## split
  x   <- strsplit(x, "\\.")[[1]]
  
  ## add information
  tmp$pheno        <- x[3]
  tmp$chr          <- as.numeric(x[4])
  tmp$region_start <- as.numeric(x[5])
  tmp$region_end   <- as.numeric(x[6])
  tmp$run          <- "yes"
  
  ## return
  return(tmp)
  
})
## combine
res.pqtl <- do.call(plyr::rbind.fill, res.pqtl)
## combine with target information
res.pqtl <- merge(subset(olink.targets, chr != 23), res.pqtl, by=c("pheno", "chr", "region_start", "region_end"), all.x = T)

## find out what needs to be rerun
res.pqtl[ is.na(run)]
## drop X-chromosome (no GWAS has been run)
write.table(res.pqtl[ is.na(run), c("pheno", "mrc_olink.id", "chr", "region_start", "region_end")], 
            "Olink.fine.map.targets.rerun.txt", sep="\t", col.names = F, row.names = F, quote = F)
## two proteins not included ARHGAP25_21255 and SMAD1_21180, two HLA not fine-mapped

## drop what did not map
res.pqtl <- res.pqtl[ !is.na(cs) ]

## create identifier
olink.targets[, fine.map.id := paste(pheno, chr, region_start, region_end, sep=".")]
res.pqtl[, fine.map.id := paste(pheno, chr, region_start, region_end, sep=".")]

## write to file
write.table(res.pqtl, "Results.fine.mapping.cis.pQTL.Olink.EPIC.MA.20230120.txt", sep="\t", row.names=F)

## create input files for phenome-wide coloc
write.table(subset(olink.targets, fine.map.id %in% res.pqtl$fine.map.id), "../../02_phewas_coloc_Olink_EPIC/input/Olink.targets.with.cis.credible.set.20230119.txt",
            sep="\t", row.names = F)
