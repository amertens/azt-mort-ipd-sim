
############################################
# fit.rma
##############################################

# Documentation: fit.rma
# Usage: fit.rma(data,age,ni,xi,measure,nlab, method = "REML")
# Description: Take an input dataframe with each row summarizing the N's and cases from a single
#             cohort and then calculated the parameter of interest specified by measure, pooled
#             across all cohorts
# age (optional): variable level to filter the data by, selecting a specific age category


# Args/Options:
# data: a data frame with variables study, country, agecat, and outcome-specific summary measures for ni and xi
# ni: name of the variable that is the total count of observations
# xi:  name of the variable
# measure: fed into rma() function; character string indicating the type of data supplied to the function. "PLO" by default for logit transformed proportion
# nlab: optional label of the count of observations i.e. "children" to be appended to the formatted counts/
# method: fed into rma() function; haracter string specifying whether a fixed- ("FE") or a random/mixed-effects model ("REML") should be fitted.

# random effects function, save results nicely

fit.rma <- function(data, ni, xi = NULL, yi = NULL, vi = NULL, measure = "PLO", nlab = "", method = "REML", age=NULL) {
  if(!is.null(age)){
    data <- data[data$agecat==age,]
    data$age <- age
  }
  
  mode_continuous <- !is.null(yi) | !is.null(vi)
  mode_binary <- !is.null(xi)
  if (mode_binary & mode_continuous) stop("can only do binary or continuous")
  # check if measure=="PLO" - default for all parameters bounded between 0 and 1 (prevalence, cumulative incidence)
  # because then output needs back transformation
  if (measure == "PLO" & mode_binary) {
    # If only one row of the data, no need for pooling (single study, often seens when looping over agecats),
    # so just use the escalc() function that underlies RMA to calculate the parameter for the single study
    if (nrow(data) == 1) {
      fit <- NULL
      try(fit <- escalc(data = data, ni = data[[ni]], xi = data[[xi]], method = method, measure = measure, append = T))
      data <- fit
      data$se <- sqrt(data$vi)
      out <- data %>%
        ungroup() %>%
        mutate(nstudies = 1, nmeas = data[[ni]]) %>%
        mutate(
          est = plogis(yi),
          lb = plogis(yi - 1.96 * se),
          ub = plogis(yi + 1.96 * se),
          nmeas.f = paste0("N=", format(sum(data[[ni]]), big.mark = ",", scientific = FALSE), " ", nlab),
          nstudy.f = paste0("N=", nstudies, " studies"),
          method.used=method,
          ptest.f = sprintf("%0.0f", est)
        ) %>%
        subset(., select =c(nstudies, nmeas, est, se, lb, ub, nmeas.f, nstudy.f, method.used, agecat, ptest.f)) %>%
        as.tibble()
      rownames(out) <- NULL
      # If input is more than 1 row (multiple studies), pool across studies with rma() function from metafor package
    } else {
      # Check if 0 cases of the outcome
      # Use FE model if all 0 counts because no heterogeneity and rma.glmm fails
      if(!is.null(xi)){
        if (sum(data[[xi]]) == 0) method <- "FE"
      }
      fit <- NULL
      method_fit <- method
      try(fit <- rma(
        data = data,
        ni = data[[ni]],
        method = method,
        xi = data[[xi]],
        measure = measure
      ))
      if(is.null(fit)){try(fit <- rma(
        data = data,
        ni = data[[ni]],
        method = "ML",
        xi = data[[xi]],
        measure = measure
      ))
        method_fit <- "ML"
      }
      if(is.null(fit)){try(fit <- rma(
        data = data,
        ni = data[[ni]],
        method = "DL",
        xi = data[[xi]],
        measure = measure
      ))
        method_fit <- "DL"
      }
      # Create formatted dataset to return
      out <- data %>%
        ungroup() %>%
        summarise(
          nstudies = length(unique(study)),
          nmeas = sum(data[[ni]])
        ) %>%
        mutate(
          est = plogis(fit$beta),
          se = fit$se,
          lb = plogis(fit$beta - 1.96 * fit$se),
          ub = plogis(fit$beta + 1.96 * fit$se),
          nmeas.f = paste0("N=", format(sum(data[[ni]]), big.mark = ",", scientific = FALSE), " ", nlab),
          nstudy.f = paste0("N=", nstudies, " studies"),
          Qstat=fit$QE,
          tau2=fit$tau2,
          I2=fit$I2,
          method.used=method_fit
        ) %>%
        as.tibble()
      rownames(out) <- NULL
    }
  } else {
    
    if(measure == "IR"){
      # If measure if IR for incidence rate, use `to="if0all"` argument to add .5 to all cells of the 2x2 table if one is 0 so rma() works
      to <- "if0all"
      fit <- NULL
      method_fit <- method
      try(fit <- rma(
        data = data,
        ti = data[[ni]],
        method = method,
        xi = data[[xi]],
        measure = measure,
        to=to
      ))
      if(is.null(fit)){try(fit <- rma(
        data = data,
        ti = data[[ni]],
        method = "ML",
        xi = data[[xi]],
        measure = measure
      ))
        method_fit <- "ML"
      }
      if(is.null(fit)){try(fit <- rma(
        data = data,
        ti = data[[ni]],
        method = "DL",
        xi = data[[xi]],
        measure = measure
      ))
        method_fit <- "DL"
      }
      
    }else{
      
      
      # If measure other than PLO or IR is chosen:
      to <- "only0"
      
      fit <- NULL
      method_fit <- method
      
      if (mode_binary) {
        try(fit <- rma(
          data = data,
          ni = data[[ni]],
          method = method,
          xi = data[[xi]],
          measure = measure,
          to=to
        ))
        if(is.null(fit)){try(fit <- rma(
          data = data,
          ni = data[[ni]],
          method = "ML",
          xi = data[[xi]],
          measure = measure
        ))
          method_fit <- "ML"
        }
        if(is.null(fit)){try(fit <- rma(
          data = data,
          ni = data[[ni]],
          method = "DL",
          xi = data[[xi]],
          measure = measure
        ))
          method_fit <- "DL"
        }
      }
      if (mode_continuous) {
        try(fit <- rma(
          data = data,
          mi = data[[yi]],
          sdi = sqrt(data[[vi]]),
          ni = data[[ni]],
          method = method,
          measure = "MN"
        ))
        if(is.null(fit)){try(fit <- rma(
          data = data,
          mi = data[[yi]],
          sdi = sqrt(data[[vi]]),
          ni = data[[ni]],
          method = "ML",
          measure = "MN"))
          method_fit <- "ML"
        }
        if(is.null(fit)){try(fit <- rma(
          data = data,
          mi = data[[yi]],
          sdi = sqrt(data[[vi]]),
          ni = data[[ni]],
          method = "DL",
          measure = "MN"))
          method_fit <- "DL"
        }
      }
    }
    out <- data %>%
      ungroup() %>%
      summarise(
        nstudies = length(unique(study)),
        nmeas = sum(data[[ni]])
      ) %>%
      mutate(
        est = fit$beta,
        se = fit$se,
        lb = fit$ci.lb,
        ub = fit$ci.ub,
        nmeas.f = paste0(
          "N=", format(sum(data[[ni]]), big.mark = ",", scientific = FALSE),
          " ", nlab
        ),
        nstudy.f = paste0("N=", nstudies, " studies"),
        method.used=method_fit,
        Qstat=fit$QE,
        tau2=fit$tau2,
        I2=fit$I2
      )
    
  }
  return(out)
}





# Meta-analysis of log IRR estimates from negative binomial models
nb_meta <- function(log_irr, se_log_irr, study_labels = NULL, method = "REML"){
  # log_irr: Vector of log IRR estimates from neg. binomial models
  # se_log_irr: Standard errors of log IRR estimates
  # study_labels: Optional vector of study identifiers
  
  # Fit random-effects model
  res=NULL
  try(res <- rma(yi = log_irr,sei = se_log_irr,method = method,test = "knha"))
  if(is.null(res)){ try(res <- rma(yi = log_irr, sei = se_log_irr, method = "ML",test = "knha"))}
  if(is.null(res)){ try(res <- rma(yi = log_irr, sei = se_log_irr, method = "DL",test = "knha"))}
  if(is.null(res)){ try(res <- rma(yi = log_irr, sei = se_log_irr, method = "HE",test = "knha"))}
  
  # Back-transform results to IRR scale
  irr <- exp(c(b = res$b,
               ci.lb = res$ci.lb,
               ci.ub = res$ci.ub))
  
  # Calculate prediction interval
  pi_lb <- res$b - qt(0.975, df = res$k-1) * 
    sqrt(res$tau2 + res$se^2)
  pi_ub <- res$b + qt(0.975, df = res$k-1) * 
    sqrt(res$tau2 + res$se^2)
  
  pi <- exp(c(pi.lb = pi_lb, pi.ub = pi_ub))
  
  return(data.frame(irr = irr[1], lb = irr[2], ub=irr[3], pi_lb = pi_lb, pval=res$pval, pi_ub=pi_ub))
}

# Example usage:
# results <- nb_meta(
#     log_irr = c(-0.35, -0.42, -0.28),  # from glm.nb() coefficients
#     se_log_irr = c(0.12, 0.15, 0.11),  # from summary(glm.nb())
#     study_labels = c("Study 1", "Study 2", "Study 3")
# )
# 
# # Examine results
# summary(results$model)
# 
# # Forest plot
# forest(results$model,
#        transf = exp,
#        refline = 1,
#        xlab = "Incidence Rate Ratio")