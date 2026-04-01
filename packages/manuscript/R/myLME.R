## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# function to input fixed factors, random factors and run the LMM
myLME <- function(fixedF, randomF, dataset, method) {
  lme.fit <- nlme::lme(fixed = as.formula(fixedF), random = as.formula(randomF), data = dataset, method = method, na.action = na.omit)
  return(lme.fit)
}
