###################################################
#### pGWAS results for fasting proteins        ####
#### Maik Pietzner                  15/02/2023 ####
###################################################

rm(list=ls())
setwd("../path")
options(stringsAsFactors = F)
load(".RData")

## packages needed
require(data.table)
require(igraph)
require(readxl)
require(MendelianRandomization)
require(doMC)

#########################################
####      import results needed      ####
#########################################

## import results from fasting study
res.fast  <- fread("ST2 in Excel")

## create indicator at which day the protein passes significance
res.fast[, day.sig := apply(res.fast[, paste0("pval.exposure", c(1:7,10)), with=F], 1, function(x){
  ii <- which(x < .05/1044)
  if(length(ii) > 0){
    return(min(ii))
  }else{
   return(99) 
  }
})]

## import coloc results
res.susie <- fread("ST4 in Excel")
## create unique identifier
res.susie[, id := paste(OlinkID, chr, region_start, region_end, cs, sep="_")]
res.susie <- unique(res.susie)

#########################################
####          create network         ####
#########################################

## create column for protein increasing effect
res.susie[, beta.ieu.incr := sign(beta.pqtl) * beta.ieu]

## create object
e.tmp       <- res.susie[, c("id", "id.ieu", "PP.H4.abf", "beta.ieu.incr", "OlinkID")] 
## add information when significant
e.tmp       <- merge(e.tmp, res.fast[, c("OlinkID", "day.sig")])
## order
e.tmp       <- e.tmp[, c("id", "id.ieu", "OlinkID", "PP.H4.abf", "beta.ieu.incr", "day.sig")]
## assign first occurence of phenotypes
e.tmp       <- as.data.table(e.tmp)
e.tmp       <- e.tmp[ order(id.ieu, day.sig)]
e.tmp[, ind := 1:.N, by="id.ieu"]
## add to nodes later
tmp         <- e.tmp[ind == 1, c("id.ieu", "day.sig")]

## now vertices
tmp1        <- data.frame(name=unique(e.tmp$id.ieu), type="phenotype") ## n = 643
## add additional information
tmp1        <- merge(tmp1, unique(res.susie[, c("id.ieu", "phenotype.ieu")]), by.x="name", by.y="id.ieu")
## add first occurence
tmp1        <- merge(tmp1, tmp, by.x="name", by.y = "id.ieu")
names(tmp1) <- c("name", "type", "label", "day.sig")
## add protein suffix to allow for overlap between risk factors and protein targets
tmp2        <- data.frame(name=unique(e.tmp$id), type="protein") ## n = 338
## add additional information
tmp2        <- merge(tmp2, unique(res.susie[, c("id", "OlinkID")]), by.x="name", by.y="id")
tmp2        <- merge(tmp2, res.fast[, c("OlinkID", "Assay", "day.sig")], by="OlinkID")
names(tmp2) <- c("OlinkID", "name", "type", "label", "day.sig")
## combine into one
v.tmp       <- plyr::rbind.fill(tmp1, tmp2)
## assign new colours to proteins (may add cluster assignment ....)
v.tmp$cl    <- ifelse(v.tmp$type == "protein", "white", "grey80")
## edit colour
col.grad    <- colorRampPalette(c("#00A4CC", "white", "#F95700"))(120)
v.tmp$cl    <- sapply(v.tmp$OlinkID, function(x){
  ## protein or not
  if(is.na(x)){
    return("white")
  }else{
    return(col.grad[60 + res.fast$b.high[which(res.fast$OlinkID == x)]*10])
  }
})
## create another colour gradient based on each day
res.fast <- as.data.frame(res.fast)
for(j in c(0:7,10)){
  v.tmp[, paste0("cl.", j)] <- sapply(v.tmp$OlinkID, function(x){
    ## protein or not
    if(is.na(x)){
      return("white")
    }else{
      return(col.grad[60 + res.fast[which(res.fast$OlinkID == x), paste0("beta.exposure", j)]*10])
    }
  })
}

## create the graph
fast.net    <- graph_from_data_frame(e.tmp, directed = F, v.tmp)

## define layout
l           <- layout_with_fr(fast.net, grid = "nogrid")
l           <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

## create PDF
pdf("../graphics/Fast.protein.risk.factors.pdf", width = 6.3, height = 6.3)
par(mar=rep(1,4), tck=-.01, cex.axis=.4, cex.lab=.4, mgp=c(.6,0,0), bty="n", yaxs="i", xaxs="i", lwd=.1)

## results from naive application
plot(fast.net,
     vertex.size=ifelse(V(fast.net)$type == "protein", 1.5, 1.5),
     vertex.color=V(fast.net)$cl,
     vertex.frame.color=ifelse(V(fast.net)$category == "Protein", "grey20", V(fast.net)$cl),
     edge.width=.1,
     edge.lty=ifelse(sign(E(fast.net)$beta.ieu.incr) == 1, 1, 2),
     ## omit edges not included in the network
     edge.color="grey50",
     vertex.label.cex = .1,
     vertex.label.color="black",
     vertex.shape=ifelse(V(fast.net)$type == "protein", "square", "circle"),
     # vertex.label.family="Helvetica",
     vertex.label = V(fast.net)$label,
     layout=l*.95,
     rescale=F)

dev.off()

#---------------------------------------#
##--      time resolved network      --##
#---------------------------------------#

pdf("../graphics/Fast.protein.risk.factors.time.resolved.version.4.pdf", width = 6.3, height = 6.3)
par(mar=rep(.5,4), mfrow=c(3,3), lwd=.05, bty="n", bg="white")

## one empty place holder
plot(1,1, type="n", xaxt="n", yaxt="n", xlab="", ylab="")

for(j in c(2:8, 100)){
  
  ## dummy plot
  plot(fast.net, layout=l*.95, rescale=F, vertex.label=NA)
  
  ## background
  pm <- par("usr")
  rect(-1.1,-1.1,1.1,1.1, border = "black", col="grey80", lwd=.5)
  
  ## new plot
  par(new=T)
  
  ## plot entire network
  plot(fast.net,
       vertex.size=ifelse(V(fast.net)$day.sig >= j, 0, ifelse(V(fast.net)$type == "protein", 3, 3)),
       vertex.color=ifelse(V(fast.net)$day.sig >= j, NA, get.vertex.attribute(fast.net)[[paste0("cl", ifelse(j == 100, "", paste0(".", j-1)))]]),
       vertex.frame.color=ifelse(V(fast.net)$day.sig >= j, NA, rgb(0,0,0,.2)),
       vertex.frame.width=.1,
       edge.width=ifelse(E(fast.net)$day.sig >= j, 0, .05),
       # edge.lty=ifelse(sign(E(fast.net)$beta.ieu.incr) == 1, 1, 2),
       edge.lty=ifelse(E(fast.net)$day.sig >= j, 0, 1),
       ## omit edges not included in the network
       edge.color=rgb(0,0,0,.4),
       vertex.label.cex = .1,
       vertex.label.color="grey90",
       vertex.shape=ifelse(V(fast.net)$type == "protein", "square", "circle"),
       # vertex.label.family="Helvetica",
       # vertex.label = ifelse(V(fast.net)$day.sig >= j, "", V(fast.net)$label),
       vertex.label=NA,
       layout=l*.95,
       rescale=F)
  
  ## surrounding rectangle
  rect(-1.1,-1.1,1.1,1.1, lwd=.5)
  
  ## legend for protein and phenotype
  if(j == 2){
    legend("bottomleft", lty=c(0,0,1), bty="n", cex=.5, pch=c(22,21, NA), pt.cex=.8,
           pt.lwd=.3, pt.bg=c("white", "white"), legend = c("Protein", "Phenotype", "PP>80%"), 
           text.col="black")
  }
  
  ## header for the figure
  if(j %in% 2:8){
    mtext(paste("Proteins changed before not until day", j-1), cex=.4, line=-.6, col="black")
  }else{
    mtext("All proteins changed", cex=.4, line=-.6, col="black")
  }
  
}
dev.off()

#---------------------------------------#
##--   export for cytoscape session  --##
#---------------------------------------#

## write relevant data to files to create cytoscape session
write.table(e.tmp, "Edges.Fasting.cis.coloc.phenome.20230221.txt", sep="\t", row.names = F, quote=F)
write.table(v.tmp, "Nodes.Fasting.cis.coloc.phenome.20230221.txt", sep="\t", row.names = F, quote=F)

#########################################
####  some reporting for the draft   ####
#########################################

#------------------------#
##-- add MR estimates --##
#------------------------#

## go through all
res.susie <- lapply(1:nrow(res.susie), function(x){
  
  ## get the relevant data
  mri <- mr_input(bx = res.susie$beta.pqtl[x], bxse = res.susie$se.pqtl[x],
                  by = res.susie$beta.ieu[x], byse = res.susie$se.ieu[x])
  ## run the IVW
  mri <- mr_ivw(mri)
  ## return relevant information
  return(data.table(res.susie[x,], beta.mr=mri@Estimate, se.mr=mri@StdError, pval.mr=mri@Pvalue))
})
## combine again
res.susie <- do.call(rbind, res.susie)

#------------------------------------------------#
##-- indicate whether fasting reverses effect --##
#------------------------------------------------#

## add other annotations
res.susie <- merge(res.susie, res.fast[, c("OlinkID", "cluster.paper", "label", "t.high", "b.high", "pval.aov.t.factor", "fdr.aov")], by="OlinkID")

## create indicator
res.susie[, ind.direction := ifelse(sign(beta.mr) == sign(b.high), "exaggerating", "compensating")]
table(res.susie$ind.direction)
# compensating exaggerating 
# 652          597 

## write to file
write.table(res.susie, "Results.cis.coloc.susie.fasting.candidates.pruned.20230222.txt", sep="\t", row.names=F)

## most frequently coloc phenotype
tail(sort(table(res.susie$phenotype.ieu)))
## Coronary artery disease
res.susie[ phenotype.ieu == "Coronary artery disease"]
res.susie[ phenotype.ieu == "Myocardial infarction"]

#------------------------------------------------#
##--      examples coronary artery disease    --##
#------------------------------------------------#

## pick out examples for coronary artery disease (ebi-a-GCST005194)
res.tmp <- res.susie[ id.ieu == "ebi-a-GCST005194"]
## order
res.tmp <- res.tmp[ order(Assay, OlinkID)]

pdf("../graphics/Fasting.protein.CAD.coloc.evidence.20230228.pdf", width = 3.15, height = 2)
## graphical parameters
par(mar=c(2,5,.5,.5), mgp=c(1,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, yaxs="i", xaxs="i", lwd=.5)
## empty plot
plot(c(-.3,.5), c(.5,nrow(res.tmp)+.5), type="n", xlab="Odds ratio (95%-CI) for CAD per 1 s.d. increase\nin genetically predicted protein levels",
     xaxt="n", yaxt="n", ylab="")
axis(1, at=log(c(.75,1,1.25)), labels = c(.75,1,1.25), lwd=.5)
## retangles to divide
pm <- par("usr")
rect(pm[1], 1:nrow(res.tmp)-.5, pm[2], 1:nrow(res.tmp)+.5, border=NA, col=c("white", "grey90"))
abline(v=0, lwd=.5)
## add confidence intervals
arrows(res.tmp$beta.mr - 1.96 * res.tmp$se.mr, 1:nrow(res.tmp), res.tmp$beta.mr + 1.96 * res.tmp$se.mr, 
       1:nrow(res.tmp), lwd=.3, length=0)
## add point estimates
points(res.tmp$beta.mr, 1:nrow(res.tmp), pch=22, lwd=.3, bg=ifelse(res.tmp$pval.mr < .05, "grey40", "grey80"), cex=.5)
## add names
text(pm[1], 1:nrow(res.tmp), labels = paste0(res.tmp$Assay, " - ", res.tmp$rsid), pos=2, xpd=NA, cex=.5)
## close device
dev.off()

## some numbers for the manuscript
jj <- which(res.susie$Assay == "HYOU1" & res.susie$phenotype.ieu == "Coronary artery disease")
exp(res.susie$beta.mr[jj]); exp(res.susie$beta.mr[jj] - 1.96*res.susie$se.mr[jj]); exp(res.susie$beta.mr[jj] + 1.96*res.susie$se.mr[jj])

## correct for case/control imbalance in the estimate (from https://www.ahajournals.org/doi/full/10.1161/CIRCRESAHA.117.312086?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org)
i.cases  <- 88122 + 34541
i.contrl <- 261984 + 166544
## convert beta estimate log(odds ratio) = β / (µ * (1 - µ)), µ = case fraction
beta.new <- res.susie$beta.ieu[jj]/(i.cases/i.contrl*(1-i.cases/i.contrl))
se.new   <- res.susie$se.ieu[jj]/(i.cases/i.contrl*(1-i.cases/i.contrl))
## mr estimate
mr.hyou1 <- mr_input(bx = res.susie$beta.pqtl[jj], bxse = res.susie$se.pqtl[jj],
                     by = beta.new, byse = se.new)
## get mr results
mr.hyou1 <- mr_ivw(mr.hyou1)
## get estiamtes for the paper
exp(mr.hyou1@Estimate); exp(mr.hyou1@Estimate - 1.96*mr.hyou1@StdError); exp(mr.hyou1@Estimate + 1.96*mr.hyou1@StdError)
mr.hyou1@Pvalue
