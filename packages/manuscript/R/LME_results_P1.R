## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# function to run linear mixed models with repeated measurements including baseline score as dependent variables
LME_results_P1 <- function(dataset, protein_names, ind_cov, ID, fixedFormula, randomFormula, method) {
  # could fit constrained LDA—assume baseline in the treatment groups (without treatment effect in the model) or regular LDA—allowed different baseline scores in the treatment groups) (Liu et al. 2009),
  # it used myLME function which could take in different fixed effect formula and random effect formula
  # print(paste("nrow of dataset",nrow(dataset)))
  varnames <- names(dataset)
  index1 <- which(varnames %in% protein_names) # index of each protein
  ID_index <- grep(ID, varnames) # index of subject/patient ID
  # print(paste("patientID index",ID_index))
  ind_cov_index <- which(varnames %in% ind_cov) # index for ind covariates
  # print(paste("ind covariate index",ind_cov_index))
  results1 <- data.frame()
  # print(dim(results1))
  j <- 0
  for (i in index1) { # loop through all proteins
    print(paste("protein index", i))
    j <- j + 1
    data2_protein <- dataset[, c(ID_index, i, ind_cov_index)]
    print(dim(data2_protein))
    print(names(data2_protein))
    colnames(data2_protein)[2] <- "protein" # to have the same name in the regression equation
    data2_protein <- data2_protein |> dplyr::filter(!is.na(protein))

    lme.fit <- myLME(dataset = data2_protein, fixedF = fixedFormula, randomF = randomFormula, method = method) # random intercept model
    Ttable <- as.data.frame(unlist(summary(lme.fit)$tTable))
    yhat.mar <- fitted(lme.fit, level = 0)
    yhat.cond <- fitted(lme.fit, level = 1)
    resid.mar <- resid(lme.fit, type = "response", level = 0) # marginal residuals
    resid.cond <- resid(lme.fit, type = "response", level = 1) # conditional residuals

    Ttable$protein <- names(dataset)[i] # protein name
    Ttable$protein_ID <- index1[j]
    Ttable$coef.names <- rownames(Ttable)
    Ttable$shapirotest.pvalue <- shapiro.test(resid.mar)$p.value
    Ttable$resid.mar.mean <- mean(resid.mar)
    Ttable$resid.mar.corr <- cor(yhat.mar, resid.mar)
    Ttable$outliers <- sum(resid.mar < (mean(resid.mar) - 1.5 * (quantile(resid.mar, 0.75) - quantile(resid.mar, 0.25))) | resid.mar > (mean(resid.mar) + 1.5 * (quantile(resid.mar, 0.75) - quantile(resid.mar, 0.25))))
    Ttable$R2 <- cor(data2_protein[, 2], yhat.mar)

    Ttable$shapirotest.pvalue.cond <- shapiro.test(resid.cond)$p.value
    Ttable$resid.cond.mean <- mean(resid.cond)
    Ttable$resid.cond.corr <- cor(yhat.cond, resid.cond)
    Ttable$outliers.cond <- sum(resid.cond < (mean(resid.cond) - 1.5 * (quantile(resid.cond, 0.75) - quantile(resid.cond, 0.25))) | resid.cond > (mean(resid.cond) + 1.5 * (quantile(resid.cond, 0.75) - quantile(resid.cond, 0.25))))
    Ttable$R2.cond <- cor(data2_protein[, 2], yhat.cond)
    Ttable$AIC <- AIC(lme.fit)
    results1 <- dplyr::bind_rows(results1, Ttable)

    # print(dim(results1))
  }
  colnames(results1)[5] <- "p.value"
  return(results1)
}
