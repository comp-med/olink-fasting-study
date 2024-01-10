###################################################
#### Analysis of Olink fasting study           ####
#### Maik Pietzner                  11/11/2022 ####
###################################################

rm(list=ls())
setwd("../path")
options(stringsAsFactors = F)
load(".RData")

## packages needed
require(data.table)
require(lmerTest)
require(igraph)
require(doMC)
require(colorspace)
require(basicPlotteR)
require(navmix)
require(gprofiler2)

########################################
####      import QCed data set      ####
########################################

## import overall data
npx.data         <- fread("Fasting.proteomics.Olink.Explore.20240109.txt")
## label for proteins
prot.label       <- fread("Fasting.proteomics.Olink.Explore.label.20240109.txt")

##########################################
####       create normalization       ####
##########################################

## create common scale across proteins, use day -4 as a reference
fast.norm <- lapply(prot.label$OlinkID, function(x){
  ## get the relevant data
  tmp <- fast.dat[which(fast.dat$t.point == -2), x]
  return(data.frame(OlinkID=x, mean.value=mean(tmp, na.rm=T), sd.value=sd(tmp, na.rm=T)))
})
## combine again
fast.norm <- do.call(rbind, fast.norm)

## apply to the data
for(j in 1:nrow(fast.norm)){
  ## scale
  fast.dat[, fast.norm$OlinkID[j]] <- (fast.dat[, fast.norm$OlinkID[j]] - fast.norm$mean.value[j])/fast.norm$sd.value[j]
}
## check whether it worked
sd(fast.dat[which(fast.dat$t.point == -2), fast.norm$OlinkID[45]])
## seem to have worked

##########################################
####       simple linear model        ####
##########################################

## load function to do so
source("../scripts/mixed_effect_regression.R")

## convert t.point to factor 
fast.dat$t.factor       <- as.factor(fast.dat$t.point) 

## time point analysis
res.fast.linear         <- mixed.anova(fast.dat, "t.factor", prot.label$OlinkID, "+ (1|participant)")
## add protein label
res.fast.linear         <- merge(prot.label, res.fast.linear, by.x="OlinkID", by.y="outcome")
## how many FDR significant overall
res.fast.linear$fdr.aov <- p.adjust(res.fast.linear$pval.aov.t.factor, method = "BH")
nrow(subset(res.fast.linear, fdr.aov < .05))
## n=1044 proteins, almost half of all proteins...

## convert to data table
res.fast.linear         <- as.data.table(res.fast.linear)

#--------------------------------------#
##--        cluster results         --##
#--------------------------------------#

## try package https://github.com/aj-grant/navmix
# devtools::install_github("aj-grant/navmix")
require(navmix)

## get relevant data
cl.data           <- as.data.frame(subset(res.fast.linear, fdr.aov < .05, select = c("OlinkID", grep("beta|se|pval", names(res.fast.linear), value=T))))
## give rownames to ease clustering
rownames(cl.data) <- cl.data$OlinkID
## compute z-scores
for(j in c(0:7, 10)){
  cl.data[, paste0("zscore.", j)] <- cl.data[, paste0("beta.exposure", j)]/cl.data[, paste0("se.exposure", j)]
}

## run navmix on the output
cluster_out               <-  navmix(as.matrix(cl.data[, paste0("zscore.", c(0:7,10))]), K=20, plot=T, plot_radial = T)
## total of 19 clusters (18 plus one for noise)
plot(1:20, cluster_out$BIC) ## indicates marginal improvement beyond 9; repeat with max of 9
## set seed to get reproducible results
set.seed(42)
cluster.fasting           <- navmix(as.matrix(cl.data[, paste0("zscore.", c(0:7,10))]), K=9, plot=T, select_K = F, reorder_traits = F)
## looks very sensible! ten proteins in the noise cluster 

## extract the assignment
tmp.cl                    <- cluster.fasting$fit$z
## add protein information
tmp.cl                    <- data.frame(OlinkID=row.names(cl.data), cluster=cluster.fasting$fit$z)
## add to fasting results (keep non-significant proteins as well!)
res.fast.linear           <- merge(res.fast.linear, tmp.cl, all.x=T)
## order again for plotting
res.fast.linear           <- res.fast.linear[order(fdr.aov)]

## --> create label for each cluster <-- ##
cl.label.fast   <- data.frame(cluster=c(1:10),
                              label=c("Late response with continuous decrease",
                                      "Early strong positive response",
                                      "Weak response, overshoot",
                                      "Response mid study",
                                      "Early decline, normal before end",
                                      "Late response with moderate decrease",
                                      "Early decline with plateau",
                                      "Late decline, overshoot",
                                      "Strong increase, sustained",
                                      "Noise"),
                              cluster.paper=c("cluster 5",
                                              "cluster 1",
                                              "cluster 8",
                                              "cluster 2",
                                              "cluster 9",
                                              "cluster 6",
                                              "cluster 4",
                                              "cluster 7",
                                              "cluster 3",
                                              "cluster 10"))

#-----------------------------------#
##--   create plots by cluster   --##
#-----------------------------------#

# pdf("../graphics/Clustering.Fasting.navmix.K.9.pdf", width = 6.3, height = 6.3)
png("../graphics/Clustering.Fasting.navmix.K.9.png", width=16, height=16, units="cm", res=900)
par(mar=c(2,2,1.5,.5), mgp=c(1,.2,0), tck=.01, cex.axis=.6, cex.lab=.6, mfrow=c(3,3), lwd=.5)

## loop over number of clusters from k-Means (omit noise cluster)
for(i in c(2,4,9,7,1,6,8,3,5)){
  
  ## get names of metabolites to plot
  mts <- subset(res.fast.linear, cluster == i)$OlinkID
  
  ## define counter
  k   <- 0
  
  ## create empty plot
  plot(seq(-2,10, length.out=200), seq(-6,6, length.out=200), type="n",
       xlab="time [days]", xaxt="n", yaxt="n",
       ylab="Z-score"
  )
  ## indicate fasting peroid
  pm <- par("usr")
  rect(0, pm[3], 7, pm[4], border = NA, col="grey95")
  abline(h=0, lwd=.4, col="black")
  ## add axis
  axis(1, at=c(-2,0,1:7,10), lwd=.5); axis(2, lwd=.5)
  
  ## create colour vector
  cl.vec <- c(ggsci::pal_startrek()(6)[c(1:4,6)], rep("grey70", length(mts)))
  
  ## start with the ones not to highlight
  if(length(mts) > 5){
    ## increase counter
    k <- 5
    ## start loop
    for(j in mts[6:length(mts)]){
      
      ## increase counter
      k     <- k+1
      ## aggregate mean value
      foo.m <- do.call(data.frame, aggregate(as.formula(paste0(j, "~ t.point")), fast.dat, mean))
      ## plot predicted values
      points(foo.m[,1], foo.m[,2], type="l", col=adjust_transparency(cl.vec[k], .3), lwd=.3) 
      
    }
  }
  
  ## reset counter
  k <- 0
  
  ## subset to remaining mts
  mts <- mts[1:(min(5,length(mts)))]
  
  ## start loop
  for(j in mts){
    
    ## increase counter
    k     <- k+1
    ## aggregate mean value
    foo.m <- do.call(data.frame, aggregate(as.formula(paste0(j, "~ t.point")), fast.dat, mean))
    ## plot predicted values
    points(foo.m[,1], foo.m[,2], type="l", col=adjust_transparency(cl.vec[k], .5), lwd=1) 
    
    ## legend
    legend("topleft", lty=1, pch=NA, cex=.5, bty="n", 
           legend = sapply(mts, function(x) prot.label$Assay[which(prot.label$OlinkID == x)]),
           col=cl.vec[1:length(mts)])
    
  }
  
  ## header
  mtext(paste0(cl.label.fast$label[which(cl.label.fast$cluster == i)], " - ", cl.label.fast$cluster.paper[which(cl.label.fast$cluster == i)]), cex=.4)
  # legend(x="bottomright", legend=i, cex=.7)
}
dev.off()

## add label to fasting results
res.fast.linear <- merge(res.fast.linear, cl.label.fast, all.x=T)

write.table(res.fast.linear, "Results.Fasting.linear.simple.20230130.txt", sep="\t", row.names=F)

##########################################
####         pathway enrichment       ####
##########################################

## load packages to perform enrichment
require(gprofiler2)

## create background list of genes
prot.back       <- na.omit(unique(unlist(lapply(prot.label$hgnc_symbol, function(x){
  x <- strsplit(x, "\\|")[[1]]
  return(x)
}))))
## n = 2937 protein coding genes

## get all protein coding genes differential
prot.diff       <- unique(unlist(lapply(subset(res.fast.linear, fdr.aov < .05)$OlinkID, function(x){
  ## get all relevant genes
  x <- subset(prot.label, OlinkID == x)$hgnc_symbol
  ## resolve multiple mappings
  return(strsplit(x, "\\|")[[1]])
}))) ## n = 1020 unqiue genes

## try g:Profiler
res.enrich.fast <- gost(query = prot.diff, 
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE, 
                        user_threshold = 0.05, correction_method = "fdr", 
                        domain_scope = "annotated", custom_bg = prot.back, 
                        numeric_ns = "", sources = c("KEGG", "REAC"), as_short_link = FALSE)
## keep only what is needed
res.enrich.fast <- as.data.table(res.enrich.fast$result)
res.enrich.fast[, fc := (intersection_size/term_size)/(query_size/effective_domain_size)]

#------------------------------#
##--     do by cluster      --##
#------------------------------#

## apply by cluster 
res.enrich.fast.cluster <- lapply(1:10, function(x){
  
  ## get all proteins unique altered in that cluster
  prot.diff  <- unique(unlist(lapply(subset(res.fast.linear, fdr.aov < .05 & cluster == x)$OlinkID, function(k){
    ## get all relevant genes
    k <- subset(prot.label, OlinkID == k)$hgnc_symbol
    ## resolve multiple mappings
    return(strsplit(k, "\\|")[[1]])
  })))
  print(prot.diff)
  ## g:Profiler
  res.enrich <- gost(query = prot.diff, 
                     organism = "hsapiens", ordered_query = FALSE, 
                     multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                     measure_underrepresentation = FALSE, evcodes = TRUE, 
                     user_threshold = 0.05, correction_method = "fdr", 
                     domain_scope = "annotated", custom_bg = prot.back, 
                     numeric_ns = "", sources = c("KEGG", "REAC"), as_short_link = FALSE)
  
  ## return results
  if(!is.null(res.enrich$result) > 0){
    return(data.frame(cluster=x, res.enrich$result))
  }
})
## combine results
res.enrich.fast.cluster <- do.call(rbind, res.enrich.fast.cluster)
## add cluster label
res.enrich.fast.cluster <- merge(cl.label.fast, res.enrich.fast.cluster)
## convert to data table
res.enrich.fast.cluster <- as.data.table(res.enrich.fast.cluster)
## add fold change
res.enrich.fast.cluster[, fc := (intersection_size/term_size)/(query_size/effective_domain_size)]

#------------------------------#
##--  create common figure  --##
#------------------------------#

source("../scripts/plot_enrichment.R")

png("../graphics/Enrichment.Fasting.by.cluster.20230111.png", units = "cm", width = 16, height = 7, res = 900)
## graphical parameters
par(mar=c(1.5,.25,1.5,.25), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=.01, mfrow=c(2,4), bty="n")

## overall
plot.enrich(res.enrich.fast, compress = T)
## add header
mtext("All proteins", cex=.4)

## cluster
for(j in 1:6){
  ## plot cluster
  plot.enrich(res.enrich.fast.cluster,  s.set = j, compress = T)
  ## add header
  mtext(cl.label.fast$label[which(cl.label.fast$cluster == j)], cex=.4)
}
dev.off()

#---------------------------#
##--    simple summary   --##
#---------------------------#

## do simple counting graph
res.count <- lapply(c(1:7,10), function(j){
  
  ## get all positive findings
  pf <- nrow(res.fast.linear[ fdr.aov < .05 & eval(as.name(paste0("pval.exposure",j))) < .05/1044 & eval(as.name(paste0("beta.exposure",j))) > 0])
  ## get all inverse findings
  nf <- nrow(res.fast.linear[ fdr.aov < .05 & eval(as.name(paste0("pval.exposure",j))) < .05/1044 & eval(as.name(paste0("beta.exposure",j))) < 0])
  ## return results
  return(data.frame(day=j, positive.find=pf, negative.pf=nf))
  
})
## combine and create figure
res.count <- do.call(rbind, res.count)

## graphical device
pdf("../graphics/Overview.fasting.proteins.20230119.pdf", width = 1.575, height = 1.575)
# png("../graphics/Overview.fasting.proteins.20230119.png", width = 8, height = 8, units = "cm", res=300)
## graphical parameters
par(mar=c(1.5,1.5,.5,.2), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.6, tck=-.01, lwd=.5)
## empty plot
plot(c(.5, nrow(res.count)+.5), c(-360, 50), xlab="Day", ylab="Number of proteins changing",
     xaxt="n", yaxt="n", type="n")
## add axis
axis(1, at=1:nrow(res.count), labels=c(1:7,10), lwd=.5)
axis(2, lwd=.5, at=c(-300,-200,-100,0,50), labels = c(300,200,100,0,50))
## add protein targets, positive change c("#00A4CC", "#F95700")
rect(1:nrow(res.count)-.4, 0, 1:nrow(res.count)+.4, res.count$positive.find, lwd=.3, col="#00A4CC")
## add protein targets, inverse change c("#00A4CC", "#F95700")
rect(1:nrow(res.count)-.4, 0, 1:nrow(res.count)+.4, -res.count$negative.pf, lwd=.3, col="#F95700")
## add baseline
abline(h=0, lwd=.5)
## add legend
legend("bottomleft", pch=22, pt.cex=.8, pt.lwd=.2, pt.bg=c("#00A4CC", "#F95700"),
       legend = c("increase", "decrease"), bty="n", cex=.5)
## close device
dev.off()

## open graphical device
png("../graphics/Summary.Fasting.volcano.20230124.png", width = 11.3, height = 11.3, res = 300, units = "cm")
par(mar=c(1.5,1.5,.5,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5)

## create colour gradient
cl.gradient <- c(RColorBrewer::brewer.pal(8, "OrRd"), "", "", "#F95700")
## plot
plot(-log10(fdr.aov) ~ b.high, res.fast.linear, pch=21, lwd=.1, cex=(res.fast.linear$t.high+1)/10,
     bg=ifelse(res.fast.linear$fdr.aov < .05, cl.gradient[res.fast.linear$t.high+1], adjust_transparency(cl.gradient[res.fast.linear$t.high+1], .5)), 
     xlab="Largest change from baseline (s.d. units)", xaxt="n", yaxt="n",
     ylab=expression(log[10](q-value)))
## abline
abline(h=-log10(.05), lwd=.5, lty=3); abline(v=0, lwd=.5, lty=2);
axis(1, lwd=.5); axis(2, lwd=.5)
legend("bottomright", cex=.5, bty="n", pch=21, pt.cex=c(1:8,11)/10, pt.lwd=.2, 
       pt.bg=cl.gradient[c(1:8,11)], legend = paste0("Day", c(0:7,10)))
## annotate some points
res.fast.linear[, tmp.diff := sqrt(b.high^2 + (-log10(fdr.aov))^2)]
## order 
res.fast.linear <- res.fast.linear[ order(-tmp.diff)]
## add those
tmp             <- res.fast.linear[1:15,]
tmp             <- tmp[ Assay != "LEP"]
addTextLabels(tmp$b.high, -log10(tmp$fdr.aov), tmp$Assay, cex.label = .4, col.background=rgb(0,0,0,.3), 
              col.label="white", lwd = .2, keepLabelsInside = F, cex.pt = .7)

## label some candidates
tmp             <- res.fast.linear[ Assay %in% c("GHRL", "BDNF", "TSHB", "FGF21", "LEP", "ADIPOQ")]
addTextLabels(tmp$b.high, -log10(tmp$fdr.aov), tmp$Assay, cex.label = .4, col.background=adjust_transparency("#00A4CC", .6), 
              col.label="white", lwd = .2, keepLabelsInside = F, cex.pt = .5)

dev.off()

##########################################
####           pathway figures        ####
##########################################

## add rownames
tmp           <- as.data.frame(res.fast.linear)
## drop duplicates
ii            <- table(tmp$ensembl_gene_id)
tmp           <- subset(tmp, ensembl_gene_id %in% names(ii[ii == 1]))
rownames(tmp) <- tmp$ensembl_gene_id

## load package
require("pathview") ## "#00A4CC", "#F95700"
tmp <- pathview(gene.data = tmp[, paste0("beta.exposure", 1:7), drop=F], pathway.id = "04610",
                species = "hsa", out.suffix = "coagulation", gene.idtype = "ENSEMBL",
                limit = c(-1,1), low = "#00A4CC", mid = "grey90", high = "#F95700", 
                node.sum = "median", cex=.5)

## second pathway
tmp <- pathview(gene.data = tmp[, paste0("beta.exposure", 1:7), drop=F], pathway.id = "04512",
                species = "hsa", out.suffix = "ECM", gene.idtype = "ENSEMBL",
                limit = c(-1,1), low = "#00A4CC", mid = "grey90", high = "#F95700", 
                node.sum = "median", cex=.5)

## third pathway
tmp <- pathview(gene.data = tmp[, paste0("beta.exposure", 1:7), drop=F], pathway.id = "04512",
                species = "hsa", out.suffix = "ECM", gene.idtype = "ENSEMBL",
                limit = c(-1,1), low = "#00A4CC", mid = "grey90", high = "#F95700", 
                node.sum = "median", cex=.5)

## fourth pathway
tmp <- pathview(gene.data = tmp[, paste0("beta.exposure", 1:7), drop=F], pathway.id = "04979",
                species = "hsa", out.suffix = "cholesterol", gene.idtype = "ENSEMBL",
                limit = c(-1,1), low = "#00A4CC", mid = "grey90", high = "#F95700", 
                node.sum = "median", cex=.5)


##########################################
####        targeted figures          ####
##########################################

#-------------------------------#
##-- proteins for CAD figure --##
#-------------------------------#

## plot time gradient for proteins
## "OID30753" "OID30731" "OID20300" "OID30212" "OID30547" "OID31416" "OID30747" "OID20188" "OID20235" "OID20621"

## only one LMOD1 (OID30212)

pdf("../graphics/CAD.fasting.proteins.trajectory.20230731.pdf", width=3.15, height=2)
par(mar=c(1.5,1.5,.5,.5), mgp=c(.6,0,0), tck=.01, cex.axis=.5, cex.lab=.5, lwd=.5)

## get names of metabolites to plot
mts <- c("OID30753", "OID30731", "OID20300", "OID30212", "OID30747", "OID20188", "OID20235", "OID20621")

## define counter
k   <- 0

## create empty plot
plot(seq(-1,12, length.out=200), seq(-5,2, length.out=200), type="n",
     xlab="time [days]", xaxt="n", yaxt="n",
     ylab="Z-score"
)
## indicate fasting peroid
pm <- par("usr")
rect(0, pm[3], 7, pm[4], border = NA, col="grey95")
abline(h=0, lwd=.4, col="black")
## add axis
axis(1, at=c(0:7,10), lwd=.5); axis(2, lwd=.5)

## create colour vector
cl.vec <- RColorBrewer::brewer.pal(8, "Set2")

## start loop
for(j in mts){
  ## increase counter
  k     <- k+1
  ## background line
  points(c(0:7,10)+k/30-.2, res.fast.linear[which(res.fast.linear$OlinkID == j), paste0("beta.exposure", c(0:7,10))], 
         pch=21, lwd=.6, col=adjust_transparency(cl.vec[k],.5), type="l", cex=.5)
  ## confidence intervals
  arrows(c(0:7,10)+k/30-.2, unlist(res.fast.linear[which(res.fast.linear$OlinkID == j), paste0("beta.exposure", c(0:7,10))] - 1.96 * res.fast.linear[which(res.fast.linear$OlinkID == j), paste0("se.exposure", c(0:7,10))]),
         c(0:7,10)+k/30-.2, unlist(res.fast.linear[which(res.fast.linear$OlinkID == j), paste0("beta.exposure", c(0:7,10))] + 1.96 * res.fast.linear[which(res.fast.linear$OlinkID == j), paste0("se.exposure", c(0:7,10))]),
         lwd=.4, length = 0, col=adjust_transparency(cl.vec[k],.5))
  ## point estimates from regression models
  points(c(0:7,10)+k/30-.2, res.fast.linear[which(res.fast.linear$OlinkID == j), paste0("beta.exposure", c(0:7,10))], 
         pch=21, lwd=.3, bg=cl.vec[k], cex=.5)
}

## add legend
# legend("bottomleft", lty=1, pch=NA, cex=.5, bty="n",
#        legend = sapply(mts, function(x) prot.label$Assay[which(prot.label$OlinkID == x)]),
#        col=cl.vec[1:length(mts)])

## order by estimates at last timepoint
tmp <- data.frame(var=mts, lab=sapply(mts, function(x) prot.label$Assay[which(prot.label$OlinkID == x)]),
                  beta=sapply(mts, function(x) res.fast.linear$beta.exposure10[which(res.fast.linear$OlinkID == x)]),
                  col=cl.vec)
## order
tmp <- tmp[order(-tmp$beta), ]
## add text
text(10.5, seq(max(tmp$beta), -1.8, length.out=nrow(tmp)), pos=4,
     labels=tmp$lab, cex=.5, col=tmp$col, offset = 0)
## add arrows
arrows(10.1, tmp$beta, 10.4, seq(max(tmp$beta), -1.8, length.out=nrow(tmp)), lwd=.3, col=tmp$col,
       length=0)


dev.off()
