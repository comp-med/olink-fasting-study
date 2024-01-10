###########################################
###########################################
#### function to compute linear regression
#### with categories and additional inter-
#### action and adjustment for age and sex

mixed.anova <- function(dat, expo, outc, formel){
  
  ## 'dat'  -- data frame containing all needed values
  ## 'expo' -- vector of exposure of interest
  ## 'outc' -- vector of outcomes of interest
  ## 'form' -- string of the right hand side of the formula (confounder and random effects)
  
  ## load packages for mixed models
  require("lmerTest")
  
  ## run across all exposures and all outcomes of interest
  res <- lapply(outc, function(i){
    
    ## loop over all exposures of interest
    tmp <- lapply(expo, function(j){
      ## create formula
      ff <- paste0(i, "~", j, formel)
      print(ff)
      ## run the model
      ml <- lmer(as.formula(ff), data=dat)
      ff <- summary(ml)$coefficients
      ## return needed information
      ff <- data.frame(outcome=i, exposure=j, ind=row.names(ff), beta=ff[,1], se=ff[,2], pval=ff[,5])
      ## reshape the data
      ff <- reshape(ff, idvar = c("outcome", "exposure"), timevar = "ind", direction = "wide")
      ## exchange names to allow combination afterwards
      names(ff) <- gsub(j, "exposure", names(ff))
      ## add number of observations
      ff$n <- nrow(na.omit(dat[, c(expo, outc)]))
      ## add ANOVA p-values
      ml                                      <- anova(ml)
      ff[, paste0("pval.aov.", rownames(ml))] <- ml$`Pr(>F)`
      return(ff)
    })
    ## combine into one data frame
    tmp <- do.call(rbind, tmp)
    return(tmp)
  })
  ## combine into one data frame
  res <- do.call(rbind, res)
  print(head(res))
  ## get only information needed
  res <- res[, c("outcome", "n", grep("exposure|aov", names(res), value=T))]
  return(res)
}