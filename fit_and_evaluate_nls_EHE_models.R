# R-script
# R version 4.0.2 (2020-06-22)
# File encoding: UTF-8
#
# Aim:  Fit and evaluation of two distinct nonlinear
#       least squares models (MM vs. INH)
#       on the substrate-dependent catalytic behavior
#       of three extracellular hydrolytic enzymes (EHE):
#       BG:   1,4-β-glucosidase
#       NAG:  1,4-β-N-acetyl-glucosaminidase
#       AP:   acid phosphomonoesterase
#   Further abbreviations:
#       MM, 2P:   common Michealis-Menten model (2-parametric)
#       INH, 3P:  Substrate-INHibition model - an extension of the
#                 Michealis-Menten model (3-parametric)
# 
# Authors:  Martin Zwanzig¹ and Alexander Tischer²
#           ¹Technische Universität Dresden
#           ²Friedrich-Schiller Universität Jena
#
# Date: 03.01.2022

# R-package requirements ----
library(drc)
library(readxl)
library(tidyr)
library(nlme)
library(nlstools)
library(dplyr)

# Self-defined functions ----

# ...used to evaluate model performance:

# root mean square error (RMSE)
# = general estimate of the absolute error (here in units of reaction rate R) (best 0, worse ∞)
RMSE <- function(fit, obs) {
  obs1 <- obs[obs > 0]; fit1 <- fit[obs > 0] # avoiding deflation by non-reaction data
  sqrt(sum((fit1 - obs1)^2) / length(obs1))
}
# mean absolute percentage error (MAPE %)
# = general estimate of the relative error (best 0, worse ∞)
MAPE <- function(fit, obs) {
  obs1 <- obs[obs > 0]; fit1 <- fit[obs > 0] # in order to avoid division by 0
  (100 / length(obs1)) * sum(abs(obs1 - fit1) / obs1)
}
# bias % (comparable between different models)
# = relative, systematic under- or overestimation (best 0, worse ∞)
biasPERCENT <- function(fit, obs) {
  obs1 <- obs[obs > 0]; fit1 <- fit[obs > 0] # in order to avoid division by 0
  (100 / length(obs1)) * sum((fit1 - obs1) / obs1)
}
# bias as the mean of model residuals (not comparable between different models)
# = absolute, systematic under- or overestimation (here in units of reaction rate R) (best 0, worse ∞)
biasMEANRES <- function(fit, obs) {mean(fit - obs)} 
# R^2, the ‘fraction of variance explained by the model’ (best 1, worse 0)
Rsquared <- function(fit, obs) {cor(fit, obs, method = "pearson")^2}
# percent relative standard error (PRSE) of model parameters
# = relative parameter identifiability / uncertainty (best 0, worse ∞)
PRSE <- function(coeff.sd, coeff.mean) {100 * coeff.sd / coeff.mean}

# ... to calculate model tipping points of the 3P model (inh)
Rmax <- function(Vmax, Km, Ki) { Vmax / (1+2*sqrt(Km/Ki)) }
SRmax <- function(Km, Ki){ sqrt(Km * Ki) }

# Prepare ----

options(max.print=10000)
# Load data
rdf <- read_xlsx("./data/pk1_hydrolases_sep.xlsx", sheet = "pk1_hydrolases_high") # ___________________CHECK !
rdf <- rdf[order(rdf$Sub_conc),]

for (ddd in 1:3){ # EHE's are analyzed one after the other
  enzyme_ddd <- c("BG", "NAG", "AP")[ddd] # one of the EHE's is choosen
  
  # Create subset based on the data of the focal EHE:
  df <- subset(rdf, select = c("Sample_ID", "Sub_conc", enzyme_ddd))
  # check data:
  str(df)
  names(df) # check names of enzymes and substrate concentration
  names(df)[1] <- "ID" # rename ID (here in line 1) as "ID"
  names(df)[2] <- "S" # rename substrate concentration (here in line 2) as "S"
  names(df)[3] <- "R" # rename activity / reaction rate (here in line 3) as "R"
  df$ID <- as.factor(df$ID) # make ID a factor
  
  if(enzyme_ddd == "AP"){ # only AP is given in micromol in the database
    df$R <- df$R * 1000 # conversion from micromol to nanomol
  }
  
  # Create container ----
  # for simulation results and evaluation
  # for MM model:
  params <- NULL # where parameter estimates will be stored
  mpstats <- NULL # model performance stats
  bootstrap_est <- NULL
  # for INH model:
  params_inh <- NULL # where parameter estimates will be stored
  mpstats_inh <- NULL # model performance stats
  bootstrap_inh_est <- NULL
  # anova
  aov_comp <- NULL
  
  # graphics showing fits and residuals of single samples are stored in one multipage pdf file:
  pdf(paste("./output/model_figures_high_Enzymes_",enzyme_ddd,"_nls.pdf", sep = ""),
      width = 6, height = 7.5)
  
  # analyze data of each EHE sample by sample (including replicates: 3 per level of S)
  ID_split <- split(df, df$ID) 
  for (i in 1:length(ID_split)){
    df_i <- data.frame(ID_split[[i]]) 
    
    print(paste(enzyme_ddd, df_i$ID[1], "- processing sample",i,"/", length(ID_split),
                "of data subset",ddd,"/ 3 - start:", Sys.time()))
    
    # Fit MM model (full S range) ----
    # ... to full data:
    drm_i <- nls(R ~ ((Vmax*S)/(Km+S)), # classical Michaelis-Menten model
                 data = df_i, start = c(Vmax = max(df_i$R), Km = max(df_i$S/3)),
                 control = nls.control(maxiter = 50, tol = 1e-05))
    x1 <- coef(drm_i)[1]
    x2 <- coef(drm_i)[2]
    
    # Fit MM model (restricted S range) ----
    # ... restricted to ranges of S showing no decrease in R
    agg <- aggregate(df_i, list(df_i$S), mean) # calculate mean R at every S-level (mean of 3 observations)
    max_pre_decrease_S <- agg$S[agg$R == max(agg$R)] # observed maximum S level not associated with a decline in average R levels
    drm_i_pdS <- nls(R ~ ((Vmax*S)/(Km+S)), # classical Michaelis-Menten model
                     data = df_i[df_i$S <= max_pre_decrease_S,], # fit to potentially restricted S range only
                     start = c(Vmax = max(df_i$R), Km = max(df_i$S/3)),
                     control = nls.control(maxiter = 50, tol = 1e-05))
    
    # Fit INH model (full S range) ----
    drm_i_inh <- NULL
    drm_i_inh <- nls(R ~ (Vmax * S / (Km + S + (S^2) / Ki)), # Michaelis-menten model extension considering substrate inhibition
                     data=df_i,
                     # algorithm = default: Gauss-Newton algorithm
                     start = list(Vmax = x1, Km = x2, Ki = max(df_i$S)),
                     trace = FALSE, # set trace = TRUE for more control
                     control = nls.control(maxiter = 100, tol = 1e-4, warnOnly = TRUE))
    
    if(!drm_i_inh$convInfo$isConv){ # if INH did not converge -> another algorithm is used
      drm_i_inh <- nls(R ~ (Vmax * S / (Km + S + (S^2) / Ki)), data=df_i,
                       start = list(Vmax = x1, Km = x2, Ki = max(df_i$S)),
                       upper = c(Inf,Inf,1e7), # upper bounds of parameter values ... CHECK upper bounds !
                       lower = c(1e-16,1e-16,1e-16), # lower bounds of parameter values
                       trace = FALSE, # set trace = TRUE for more control
                       control = nls.control(maxiter = 100, tol = 1e-4, warnOnly = TRUE),
                       algorithm = "port") # ‘nl2sol’ algorithm from the Port library
      if(!drm_i_inh$convInfo$isConv){ # if INH still not converges -> set initial Ki to 1e7 (upper bound)
        drm_i_inh <- nls(R ~ (Vmax * S / (Km + S + (S^2) / Ki)), data=df_i,
                         start=list(Vmax = x1, Km = x2, Ki = 1e7),
                         upper = c(Inf,Inf,1e7), # upper bounds of parameter values
                         lower = c(1e-16,1e-16,1e-16), # lower bounds of parameter values
                         trace = FALSE, # set trace = TRUE for more control
                         control = nls.control(maxiter = 100, tol = 1e-4),
                         algorithm = "port") # ‘nl2sol’ algorithm from the Port library
      }
    }
    
    # Check if the INH model could converge or not
    if(!drm_i_inh$convInfo$isConv){drm_i_inh <- NULL}
    print(paste("3P INH model convergence: ",drm_i_inh$convInfo$isConv))
    
    # ANOVA ----
    INHbetter_aov_pvalue <- NA
    if(!is.null(drm_i_inh)){ # if the INH model could converge
      aob <- anova(drm_i, drm_i_inh)
      aov_comp <- rbind(aov_comp, data.frame(
        Enzyme = enzyme_ddd, # Name of enzyme
        sample = df_i$ID[1], # respective sample
        maxS = max(df_i$S), # maximum S concentration tested
        max_pre_decrease_S = max_pre_decrease_S, # maximum S concentration before R drops (if at all)
        Sinhibition = max_pre_decrease_S < max(df_i$S), # TRUE if R drops at high S
        RSS_MM = aob$`Res.Sum Sq`[1], # RSS of MM model
        RSS_INH = aob$`Res.Sum Sq`[2], # RSS of INH model
        pvalue = aob$`Pr(>F)`[2], # p value indicating if INH model is better than MM
        INHbetter = aob$`Pr(>F)`[2] < 0.05) # p value indicating if INH model is better than MM
      )
      INHbetter_aov_pvalue <- aob$`Pr(>F)`[2]
    }
    
    # Statistics MM model ----
    
    # ... parameter estimates ----
    params <- rbind(params, data.frame(Enzyme = enzyme_ddd, # Name of enzyme
                                       sample = df_i$ID[1], # respective sample
                                       Vmax = summary(drm_i)$parameters["Vmax","Estimate"],
                                       Vmax_se =  summary(drm_i)$parameters["Vmax","Std. Error"],
                                       Km = summary(drm_i)$parameters["Km","Estimate"],
                                       Km_se = summary(drm_i)$parameters["Km","Std. Error"],
                                       max_pre_decrease_S = max_pre_decrease_S,
                                       Sinhibition = max_pre_decrease_S < max(df_i$S),
                                       INHbetter_aov_pvalue = INHbetter_aov_pvalue))
    
    # ... evaluation critieria ----
    PRSEs <- PRSE(coeff.sd = summary(drm_i)$coefficients[,2],
                  coeff.mean = summary(drm_i)$coefficients[,1])
    mpstats <- rbind(mpstats, data.frame(Enzyme = enzyme_ddd, # Name of enzyme
                                         sample = df_i$ID[1], # respective sample
                                         RMSE = RMSE(fit = predict(drm_i), obs = df_i$R),
                                         RMSEp1 = ifelse(!is.null(drm_i_inh), RMSE(fit = predict(drm_i)[df_i$S <= max_pre_decrease_S], obs = df_i$R[df_i$S <= max_pre_decrease_S]), NA),
                                         RMSEp2 = ifelse(!is.null(drm_i_inh), RMSE(fit = predict(drm_i)[df_i$S > max_pre_decrease_S], obs = df_i$R[df_i$S > max_pre_decrease_S]), NA),
                                         MAPE = MAPE(fit = predict(drm_i), obs = df_i$R),
                                         MAPEp1 = ifelse(!is.null(drm_i_inh), MAPE(fit = predict(drm_i)[df_i$S <= max_pre_decrease_S], obs = df_i$R[df_i$S <= max_pre_decrease_S]), NA),
                                         MAPEp2 = ifelse(!is.null(drm_i_inh), MAPE(fit = predict(drm_i)[df_i$S > max_pre_decrease_S], obs = df_i$R[df_i$S > max_pre_decrease_S]), NA),
                                         biasPERCENT = biasPERCENT(fit = predict(drm_i), obs = df_i$R),
                                         biasPERCENTp1 = ifelse(!is.null(drm_i_inh), biasPERCENT(fit = predict(drm_i)[df_i$S <= max_pre_decrease_S], obs = df_i$R[df_i$S <= max_pre_decrease_S]), NA),
                                         biasPERCENTp2 = ifelse(!is.null(drm_i_inh), biasPERCENT(fit = predict(drm_i)[df_i$S > max_pre_decrease_S], obs = df_i$R[df_i$S > max_pre_decrease_S]), NA),
                                         biasMEANRES = biasMEANRES(fit = predict(drm_i), obs = df_i$R),
                                         biasMEANRESp1 = ifelse(!is.null(drm_i_inh), biasMEANRES(fit = predict(drm_i)[df_i$S <= max_pre_decrease_S], obs = df_i$R[df_i$S <= max_pre_decrease_S]), NA),
                                         biasMEANRESp2 = ifelse(!is.null(drm_i_inh), biasMEANRES(fit = predict(drm_i)[df_i$S > max_pre_decrease_S], obs = df_i$R[df_i$S > max_pre_decrease_S]), NA),
                                         Rsquared = Rsquared(fit = predict(drm_i), obs = df_i$R),
                                         PRSEmax = max(PRSEs),
                                         PRSE_Vmax = PRSEs["Vmax"],
                                         PRSE_Km = PRSEs["Km"],
                                         RMSEp1pdS = RMSE(fit = predict(drm_i_pdS)[df_i$S <= max_pre_decrease_S], obs = df_i$R[df_i$S <= max_pre_decrease_S]),
                                         MAPEp1pdS = MAPE(fit = predict(drm_i_pdS)[df_i$S <= max_pre_decrease_S], obs = df_i$R[df_i$S <= max_pre_decrease_S]),
                                         biasPERCENTp1pdS = biasPERCENT(fit = predict(drm_i_pdS)[df_i$S <= max_pre_decrease_S], obs = df_i$R[df_i$S <= max_pre_decrease_S]),
                                         biasMEANRESp1pdS = biasMEANRES(fit = predict(drm_i_pdS)[df_i$S <= max_pre_decrease_S], obs = df_i$R[df_i$S <= max_pre_decrease_S]),
                                         Sinhibition = max_pre_decrease_S < max(df_i$S),
                                         INHbetter_aov_pvalue = INHbetter_aov_pvalue))
    
    # ... bootstrapping ----
    bootstrap_i <-  nlsBoot(drm_i, niter = 999) # bootstrapped parameter estimation of 2P model
    bootstrap_est <- rbind(bootstrap_est, data.frame(Enzyme = enzyme_ddd, # Name of enzyme
                                                     sample = df_i$ID[1], # respective sample
                                                     Vmax.boot = bootstrap_i$estiboot[1,1],
                                                     Vmax.boot.SE = bootstrap_i$estiboot[1,2],
                                                     Vmax.boot.median = bootstrap_i$bootCI[1,1],
                                                     Vmax.boot.CIlow = bootstrap_i$bootCI[1,2],# 2.5% quantile
                                                     Vmax.boot.CIhigh = bootstrap_i$bootCI[1,3],# 97.5% quantile
                                                     Km.boot = bootstrap_i$estiboot[2,1],
                                                     Km.boot.SE = bootstrap_i$estiboot[2,2],
                                                     Km.boot.median = bootstrap_i$bootCI[2,1],
                                                     Km.boot.CIlow = bootstrap_i$bootCI[2,2],# 2.5% quantile
                                                     Km.boot.CIhigh = bootstrap_i$bootCI[2,3]))# 97.5% quantile
    
    
    
    # Statistics INH model ----
    
    if(!is.null(drm_i_inh)){ # if the INH model could converge
      
      # ... parameter estimates ----
      
      # calculate Rmax (see 'Rmax'-function on top)
      Rmax_i <- Rmax(Vmax = coef(drm_i_inh)["Vmax"],
                     Km = coef(drm_i_inh)["Km"],
                     Ki = coef(drm_i_inh)["Ki"])
      # calculate SRmax (see 'SRmax'-function on top)
      SRmax_i <- SRmax(Km = coef(drm_i_inh)["Km"],
                       Ki = coef(drm_i_inh)["Ki"])
      Sseq <- seq(0, SRmax_i, length.out = 10000)
      SRmax_i_05 <- NA
      SRmax_i_05 <- Sseq[match(TRUE, predict(drm_i_inh, newdata = data.frame(S = Sseq)) > Rmax_i/2)] # the first match of increasing R
      Sseq_highS <- seq(SRmax_i, SRmax_i * 1000, length.out = 100000)
      SRmax_i_05_highS <- NA
      SRmax_i_05_highS <- Sseq_highS[match(TRUE, predict(drm_i_inh, newdata = data.frame(S = Sseq_highS)) <= Rmax_i/2)] # the first match of decreasing R
      if(is.na(SRmax_i_05_highS)){ # if no match so far, try even higher values
        Sseq_highS <- seq(SRmax_i * 1000, SRmax_i * 100000, length.out = 100000)
        SRmax_i_05_highS <- Sseq_highS[match(TRUE, predict(drm_i_inh, newdata = data.frame(S = Sseq_highS)) <= Rmax_i/2)] # the first match of decreasing R
      }
      params_inh <- rbind(params_inh, data.frame(Enzyme = enzyme_ddd, # Name of enzyme
                                                 sample = df_i$ID[1], # respective sample
                                                 Vmax = summary(drm_i_inh)$parameters["Vmax","Estimate"],
                                                 Vmax_se =  summary(drm_i_inh)$parameters["Vmax","Std. Error"],
                                                 Km = summary(drm_i_inh)$parameters["Km","Estimate"],
                                                 Km_se = summary(drm_i_inh)$parameters["Km","Std. Error"],
                                                 Ki = summary(drm_i_inh)$parameters["Ki","Estimate"],
                                                 Ki_se = summary(drm_i_inh)$parameters["Ki","Std. Error"],
                                                 Rmax = Rmax_i,
                                                 SRmax = SRmax_i,
                                                 SRmax05 = SRmax_i_05,
                                                 SRmax05highS = SRmax_i_05_highS,
                                                 max_pre_decrease_S = max_pre_decrease_S,
                                                 Sinhibition = max_pre_decrease_S < max(df_i$S),
                                                 INHbetter_aov_pvalue = INHbetter_aov_pvalue))
      
      
      # ... evaluation critieria ----
      
      abc <- NULL
      abc <- summary(drm_i_inh)
      if(!is.null(abc)){
        PRSEs_inh <- PRSE(coeff.sd = summary(drm_i_inh)$coefficients[,2], coeff.mean = summary(drm_i_inh)$coefficients[,1])
      }else{PRSEs_inh <- data.frame(Vmax = NA, Km = NA, Ki = NA)}
      mpstats_inh <- rbind(mpstats_inh, data.frame(Enzyme = enzyme_ddd, # Name of enzyme
                                                   sample = df_i$ID[1], # respective sample
                                                   RMSE = RMSE(fit = predict(drm_i_inh), obs = df_i$R),
                                                   RMSEp1 = RMSE(fit = predict(drm_i_inh)[df_i$S <= max_pre_decrease_S], obs = df_i$R[df_i$S <= max_pre_decrease_S]),
                                                   RMSEp2 = RMSE(fit = predict(drm_i_inh)[df_i$S > max_pre_decrease_S], obs = df_i$R[df_i$S > max_pre_decrease_S]),
                                                   MAPE = MAPE(fit = predict(drm_i_inh), obs = df_i$R),
                                                   MAPEp1 = MAPE(fit = predict(drm_i_inh)[df_i$S <= max_pre_decrease_S], obs = df_i$R[df_i$S <= max_pre_decrease_S]),
                                                   MAPEp2 = MAPE(fit = predict(drm_i_inh)[df_i$S > max_pre_decrease_S], obs = df_i$R[df_i$S > max_pre_decrease_S]),
                                                   biasPERCENT = biasPERCENT(fit = predict(drm_i_inh), obs = df_i$R),
                                                   biasPERCENTp1 = biasPERCENT(fit = predict(drm_i_inh)[df_i$S <= max_pre_decrease_S], obs = df_i$R[df_i$S <= max_pre_decrease_S]),
                                                   biasPERCENTp2 = biasPERCENT(fit = predict(drm_i_inh)[df_i$S > max_pre_decrease_S], obs = df_i$R[df_i$S > max_pre_decrease_S]),
                                                   biasMEANRES = biasMEANRES(fit = predict(drm_i_inh), obs = df_i$R),
                                                   biasMEANRESp1 = biasMEANRES(fit = predict(drm_i_inh)[df_i$S <= max_pre_decrease_S], obs = df_i$R[df_i$S <= max_pre_decrease_S]),
                                                   biasMEANRESp2 = biasMEANRES(fit = predict(drm_i_inh)[df_i$S > max_pre_decrease_S], obs = df_i$R[df_i$S > max_pre_decrease_S]),
                                                   Rsquared = Rsquared(fit = predict(drm_i_inh), obs = df_i$R),
                                                   PRSEmax = max(PRSEs_inh),
                                                   PRSE_Vmax = PRSEs_inh["Vmax"],
                                                   PRSE_Km = PRSEs_inh["Km"],
                                                   PRSE_Ki = PRSEs_inh["Ki"],
                                                   RMSEp1pdS = NA,
                                                   MAPEp1pdS = NA,
                                                   biasPERCENTp1pdS = NA,
                                                   biasMEANRESp1pdS = NA,
                                                   Sinhibition = max_pre_decrease_S < max(df_i$S),
                                                   INHbetter_aov_pvalue = INHbetter_aov_pvalue))
      
      # ... bootstrapping ----
      # by default, bootstrapping is not performed due to convergence problems of too many samples, in particular of AP
      #ifelse(enzyme_ddd == "AP", bootstrap_i <- FALSE, bootstrap_i <- TRUE)
      bootstrap_i = FALSE
      if (bootstrap_i){
        bootstrap_i_inh <-  nlsBoot(drm_i_inh, niter = 999) # bootstrapped parameter estimation of 3P model
        bootstrap_inh_est <- rbind(bootstrap_inh_est, data.frame(Enzyme = enzyme_ddd, # Name of enzyme
                                                                 sample = df_i$ID[1], # respective sample
                                                                 Vmax.boot = bootstrap_i_inh$estiboot[1,1],
                                                                 Vmax.boot.SE = bootstrap_i_inh$estiboot[1,2],
                                                                 Vmax.boot.median = bootstrap_i_inh$bootCI[1,1],
                                                                 Vmax.boot.CIlow = bootstrap_i_inh$bootCI[1,2],# 2.5% quantile
                                                                 Vmax.boot.CIhigh = bootstrap_i_inh$bootCI[1,3],# 97.5% quantile
                                                                 Km.boot = bootstrap_i_inh$estiboot[2,1],
                                                                 Km.boot.SE = bootstrap_i_inh$estiboot[2,2],
                                                                 Km.boot.median = bootstrap_i_inh$bootCI[2,1],
                                                                 Km.boot.CIlow = bootstrap_i_inh$bootCI[2,2], # 2.5% quantile
                                                                 Km.boot.CIhigh = bootstrap_i_inh$bootCI[2,3],# 97.5% quantile
                                                                 Ki.boot = bootstrap_i_inh$estiboot[3,1],
                                                                 Ki.boot.SE = bootstrap_i_inh$estiboot[3,2],
                                                                 Ki.boot.median = bootstrap_i_inh$bootCI[3,1],
                                                                 Ki.boot.CIlow = bootstrap_i_inh$bootCI[3,2],# 2.5% quantile
                                                                 Ki.boot.CIhigh = bootstrap_i_inh$bootCI[3,3]))# 97.5% quantile
      }
    }
    
    # Graphical check ----
    par(mfrow = c(2,2), oma = c(0,0,1,0)) # allow 2x2 subplots
    # ... MM model behaviour ----
    plot(y = df_i$R, x=df_i$S, main = "Model fit",
         ylab = paste("Reaction rate",enzyme_ddd),
         xlab = paste("Substrate concentration ",enzyme_ddd,"-MUF", sep=""))
    abline(v = max_pre_decrease_S, col = "lightgrey", lty = 2, lwd = 0.5)
    nd <- data.frame(S = seq(from = min(df_i$S), to = max(df_i$S), length.out = 101)) # complete S range
    lines(nd$S, predict(drm_i, newdata = nd), col = "red", lwd = 2) # complete S range model
    title(sub = paste("MAPE = ", signif(mpstats$MAPE[i], 4)), line = -1, col.sub = "red") # complete S range
    nd_pds <- data.frame(S = seq(from = min(df_i$S), to = max_pre_decrease_S, length.out = 101)) # low S range
    lines(nd_pds$S, predict(drm_i_pdS, newdata = nd_pds), col = "orange", lwd = 2) # low S range model
    title(sub = paste("MAPE = ", signif(mpstats$MAPEp1pdS[i], 4)), line = -2, col.sub = "orange") # low S range model stats
    # ... MM model residuals ----
    plot(y =  residuals(drm_i, type = "pearson"), x = fitted(drm_i), col = "red",
         main = "Model residuals", xlab = "Fitted values", ylab = "Residuals") # second subplot: model residuals
    points(y =  residuals(drm_i_pdS, type = "pearson"), x = fitted(drm_i_pdS), col = "orange") # residuals of fit to low S range
    abline(h = 0, lty = 3)
    lines(y = predict(loess(y ~ x, data = data.frame(x = fitted(drm_i), y = residuals(drm_i, type = "pearson"))),
                      data.frame(x = fitted(drm_i))), x = fitted(drm_i), col = "#FF000075")
    title(sub = paste("RMSE = ", signif(mpstats$RMSE[i], 4)), line = -1, col.sub = "red")
    title(sub = paste("RMSE = ", signif(mpstats$RMSEp1pdS[i], 4)), line = -2, col.sub = "orange")
    # ... INH model behaviour ----
    plot(y = df_i$R, x=df_i$S, main = "Model fit",
         ylab = paste("Reaction rate",enzyme_ddd),
         xlab = paste("Substrate concentration ",enzyme_ddd,"-MUF", sep=""))
    abline(v = max_pre_decrease_S, col = "lightgrey", lty = 2, lwd = 0.5)
    if(!is.null(drm_i_inh)){
      lines(nd$S, predict(drm_i_inh, newdata = nd), col = "blue", lwd = 2)# first subplot: model
      title(sub = paste("MAPE = ", signif(mpstats_inh$MAPE[i], 4)), line = -1, col.sub = "blue")
      lines(nd_pds$S, predict(drm_i_inh, newdata = nd_pds), col = "skyblue", lwd = 2, lty = 3)# first subplot: model
      title(sub = paste("MAPE = ", signif(mpstats_inh$MAPEp1[i], 4)), line = -2, col.sub = "skyblue")
      # ... INH model residuals ----
      plot(y =  residuals(drm_i_inh, type = "pearson"), x = fitted(drm_i_inh), main = "Model residuals", xlab = "Fitted values", ylab = "Residuals") # second subplot: model residuals
      abline(h = 0, lty = 3)
      lines(y = predict(loess(y ~ x, data = data.frame(x = fitted(drm_i_inh), y = residuals(drm_i_inh, type = "pearson"))),
                        data.frame(x = fitted(drm_i_inh))), x = fitted(drm_i_inh), col = "#FF000075")
      title(sub = paste("RMSE = ", signif(mpstats_inh$RMSE[i], 4)), line = -1, col.sub = "blue")
      title(sub = paste("RMSE = ", signif(mpstats_inh$RMSEp1[i], 4)), line = -2, col.sub = "skyblue")
    }else{plot(1,1,col="#00000000",xaxt="n",yaxt="n",xlab="",ylab="");text(1,1,"3P model not converged")}
    title(main = df_i$ID[1], outer = TRUE, line = -1)
  }
  dev.off() # to close the device and "save" the figures put to the "model_fit.pdf" (see working directory)
  
  # Write results to csv ----
  write.csv(params, file = paste("./output/params_high_",enzyme_ddd,"_nls.csv",sep=""), row.names = FALSE) 
  write.csv(mpstats, file = paste("./output/mpstats_high_",enzyme_ddd,"_nls.csv",sep=""), row.names = FALSE)
  write.csv(params_inh, file = paste("./output/params_high_",enzyme_ddd,"_nls_inh.csv",sep=""), row.names = FALSE) 
  write.csv(mpstats_inh, file = paste("./output/mpstats_high_",enzyme_ddd,"_nls_inh.csv",sep=""), row.names = FALSE)
  write.csv(bootstrap_est, file = paste("./output/bootstrap_high_",enzyme_ddd,"_nls.csv",sep=""), row.names = FALSE)
  if (bootstrap_i){ write.csv(bootstrap_inh_est, file = paste("./output/bootstrap_high_",enzyme_ddd,"_nls_inh.csv",sep=""), row.names = FALSE) }
  write.csv(aov_comp, file = paste("./output/aov_high_",enzyme_ddd,".csv",sep=""), row.names = FALSE)
}

