#-----------------------------------------------------------------------------------------------------------------------------------
# Wrapper function to produce design-consistent estimates from including survey weights in stereotype model with rrvglm()
# Reference paper: "Applying survey-weighted proportional odds models for unbiased inference in outcome-dependent samples with ordinal outcomes"
# Author of the code: Aya Mitani
# Version: v1.0
#------------------------------------------------------------------------------------------------------------------------------------

### FUNCTION ###

svy_sm <- function(formula, design){
  
  vars <- intersect(all.vars(formula), colnames(design))
  surveydata <- model.frame(design)[, vars, drop=FALSE]
  
  pwts <- weights(design, "sampling")
  surveydata$pwts <- pwts
  
  p <- length(attr(terms(formula), "term.labels")) # number of predictors in formula
  
  attach(surveydata)
  fit <- rrvglm(formula, multinomial, weights = pwts)
  detach(surveydata)
  sfit <- summary(fit)
  coefs <- sfit@coef3[,1]
  invinf <- sfit@cov.unscaled
  scores <- weights(fit, deriv = TRUE, type = "working")$deriv
  cons <- do.call(cbind, constraints(fit))
  mf <- model.frame(formula, data = as.data.frame(surveydata))
  Y <- model.response(mf)
  K <- nlevels(factor(Y)) # number of categories for the response variable
  mmat <- model.matrix(fit)
  case_index <- as.numeric(gsub("^(.+):.*", "\\1", rownames(mmat)))
  mmatsum <- t(t(rowsum(mmat, case_index, reorder = FALSE )) / colSums(cons))
  scorethetabeta <- (((scores / pwts) %*% cons) * mmatsum)
  betahat <- coef(fit)[K:(K+(p-1))]
  xmat <- as.matrix(fit@x[,-1])
  xbeta <- xmat %*% betahat
  scorephi <- residuals(fit, type = "response")[,c(-1,-K)] * as.vector(xbeta)
  scores <- cbind(scorephi, scorethetabeta)
  inffuns <- (scores) %*% invinf
  varmat <- vcov(svytotal(inffuns, design))
  
  rval <- list(coef = coefs, fit = fit, var = varmat, design = design, call = sys.call())
  rval
  
}

### END OF FUNCTION ###

