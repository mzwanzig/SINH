library(drc)
library(readxl)
library(tidyr)
library(nlme)
library(nlstools)
library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)
library(scales)

inh_f <- function(S, Vmax, Km, Ki){
  R = Vmax * S / (Km + S *(1+S/Ki))  # 3 Parameter Model
  return(R)
}

mm_f <- function(S, Vmax, Km){
  R = (Vmax*S)/(Km+S)  # Michaelis-Menten function
  return(R)
}

# self defined functions to evaluate model performance

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

# functions to calculate model tipping points of the 3P model (inh)
Rmax <- function(Vmax, Km, Ki) { Vmax / (1+2*sqrt(Km/Ki)) }
SRmax <- function(Km, Ki){ sqrt(Km * Ki) }


drm_i <- nls(R ~ ((Vmax*S)/(Km+S)), # 2P MODEL
             data = df_i, start = c(Vmax = max(df_i$R), Km = max(df_i$S/3)),
             control = nls.control(maxiter = 50, tol = 1e-05))
x1 <- Vmax_i <- coef(drm_i)[1] # Vmax estimate
x2 <- Km_i <- coef(drm_i)[2] # Km estimate
# 2P MODEL, restricted to ranges of S showing no decrease in R
agg <- aggregate(df_i, list(df_i$S), mean)
max_pre_decrease_S <- agg$S[agg$R == max(agg$R)]
drm_i_pdS <- nls(R ~ ((Vmax*S)/(Km+S)),
                 data = df_i[df_i$S <= max_pre_decrease_S,], # restricted to ranges of S showing no decrease in R
                 start = c(Vmax = max(df_i$R), Km = max(df_i$S/3)),
                 control = nls.control(maxiter = 50, tol = 1e-05))
# 3P MODEL
drm_i_inh <- nls(R ~ (Vmax * S / (Km + S + (S^2) / Ki)), 
                 data = df_i, start = list(Vmax = x1, Km = x2, Ki = max(df_i$S)), # using 2P Vmax and Km estimates as initial guesses
                 trace = FALSE, control = nls.control(maxiter = 100, tol = 1e-04, warnOnly = TRUE))
if(!drm_i_inh$convInfo$isConv){ # not converged? -> try another algorithm!
  drm_i_inh <- nls(R ~ (Vmax * S / (Km + S + (S^2) / Ki)), data=df_i,
                   start = list(Vmax = x1, Km = x2, Ki = max(df_i$S)),
                   upper = c(Inf,Inf,1e7), # definition of upper parameter bounds ... CHECK upper bound
                   lower = c(1e-16,1e-16,1e-16), # definition of lower parameter bounds
                   trace = FALSE, # set trace = TRUE for more control
                   control = nls.control(maxiter = 100, tol = 1e-4, warnOnly = TRUE),
                   algorithm = "port") # ‘nl2sol’ algorithm from the Port library
  if(!drm_i_inh$convInfo$isConv){ # if INH did not converge -> another algorithm is used
    drm_i_inh <- nls(R ~ (Vmax * S / (Km + S + (S^2) / Ki)), data=df_i,
                     start = list(Vmax=x1, Km=x2, Ki=1e7),
                     upper = c(Inf,Inf,1e7), # definition of upper parameter bounds ... CHECK upper bound
                     lower = c(1e-16,1e-16,1e-16), # definition of lower parameter bounds
                     trace = FALSE, # set trace = TRUE for more control
                     control = nls.control(maxiter = 100, tol = 1e-4),
                     algorithm = "port") # ‘nl2sol’ algorithm from the Port library
  }
}

if(!drm_i_inh$convInfo$isConv){drm_i_inh <- NULL}
if(!is.null(drm_i_inh)){
  Rmax_i <- Rmax(Vmax = coef(drm_i_inh)[1],
                 Km = coef(drm_i_inh)[2],
                 Ki = coef(drm_i_inh)[3])
  SRmax_i <- SRmax(Km = coef(drm_i_inh)[2],
                   Ki = coef(drm_i_inh)[3])
  Sseq <- seq(0, SRmax_i, length.out = 10000)
  SRmax_i_05 <- Sseq[match(TRUE, predict(drm_i_inh, newdata = data.frame(S = Sseq)) > Rmax_i/2)]
}

maxSpred <- max(df_i$S) * 1.5

nd <- data.frame(S = seq(from = min(df_i$S), to = max(df_i$S), length.out = 101)) # complete S range
pc_MM <- data.frame(S = nd$S, R = predict(drm_i, newdata = nd), Model = "MM")
if(!is.null(drm_i_inh)){pc_INH <- data.frame(S = nd$S, R = predict(drm_i_inh, newdata = nd), Model = "INH")}else{pc_INH = NULL}
nd_pds <- data.frame(S = seq(from = min(df_i$S), to = max_pre_decrease_S, length.out = 101)) # low S range
pc_MMpds <- data.frame(S = nd_pds$S, R = predict(drm_i_pdS, newdata = nd_pds), Model = "MM (low S fit)")
pc <- rbind(pc_MM, pc_INH, pc_MMpds) # , pc_INHpds
if(!is.null(drm_i_inh)){
  pc$Model <- factor(pc$Model, levels = c("MM", "MM (low S fit)", "INH"))
}else{
  pc$Model <- factor(pc$Model, levels = c("MM", "MM (low S fit)"))
}

ndx <- data.frame(S = seq(from = max(df_i$S), to = maxSpred, length.out = 101)) # complete S range
pc_MMx <- data.frame(S = ndx$S, R = predict(drm_i, newdata = ndx), Model = "MM")
if(!is.null(drm_i_inh)){pc_INHx <- data.frame(S = ndx$S, R = predict(drm_i_inh, newdata = ndx), Model = "INH")}else{pc_INHx = NULL}
nd_pdsx <- data.frame(S = seq(from = max_pre_decrease_S, to = maxSpred, length.out = 101)) # low S range
pc_MMpdsx <- data.frame(S = nd_pdsx$S, R = predict(drm_i_pdS, newdata = nd_pdsx), Model = "MM (low S fit)")
pcx <- rbind(pc_MMx, pc_INHx, pc_MMpdsx) # , pc_INHpds
if(!is.null(drm_i_inh)){
  pcx$Model <- factor(pcx$Model, levels = c("MM", "MM (low S fit)", "INH"))
}else{
  pcx$Model <- factor(pcx$Model, levels = c("MM", "MM (low S fit)"))
}

mmbp <- mmbpmm <- inhbp <- inhbpmm <- NULL
boot = TRUE
if (boot) {
  srange = seq(from = 0, to = maxSpred, length.out = 101)
  bootstrap_i <- nlsBoot(drm_i, niter = 999) # bootstrapped parameter estimation of 2P model
  for (bb in 1:999){
    mmbp <- rbind(mmbp, data.frame(S = srange,
                                   R = mm_f(S = srange,
                                            Vmax = bootstrap_i$coefboot[bb,"Vmax"],
                                            Km = bootstrap_i$coefboot[bb,"Km"]),
                                   Model = "MM",
                                   Boot = paste("B",bb,sep="")))
    mmbp$Boot <- factor(mmbp$Boot)
  }
  if(!is.null(drm_i_inh)){
    bootstrap_i_inh <-  nlsBoot(drm_i_inh, niter = 999) # bootstrapped parameter estimation of 3P model
    for (bb in 1:999){
      inhbp <- rbind(inhbp, data.frame(S = srange,
                                       R= inh_f(S = srange,
                                                Vmax = bootstrap_i_inh$coefboot[bb,"Vmax"],
                                                Km = bootstrap_i_inh$coefboot[bb,"Km"],
                                                Ki = bootstrap_i_inh$coefboot[bb,"Ki"]),
                                       Model = "INH",
                                       Boot = paste("B",bb,sep="")))
      inhbp$Boot <- factor(inhbp$Boot)
    }
  }
}

mcol <- c(viridis::plasma(n=5, alpha = 0.66)[2],
          viridis::plasma(n=5, alpha = 0.33)[2],
          viridis::plasma(n=5, alpha = 0.66)[4])
if(is.null(drm_i_inh)){mcol <- mcol[1:2]}

# Create residual plots
prMM <- ggplot(data = data.frame(y = residuals(drm_i, type = "pearson")[-c(1:3)],
                                 x = fitted(drm_i)[-c(1:3)]),
               aes(x = x, y = y)) +
  geom_hline(yintercept = 0) +
  geom_smooth(method = "loess", se = TRUE, formula = y ~ x, lwd = 0.5, col = mcol[1]) +
  geom_point(col = mcol[1], stroke = 0, size = 2) +
  scale_x_continuous(n.breaks = 3) + scale_y_continuous(n.breaks = 5) +
  xlab("Fitted") + ylab("Pearson residuals") +
  theme_bw() + #theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  annotation_custom(grobTree(textGrob(paste("RMSE =",round(RMSE(fit = fitted(drm_i), obs = df_i$R), 0)),
                                      x=0.1,  y=0.9, hjust=0,
                                      gp=gpar(col=viridis::plasma(n=5, alpha = 0.9)[2],
                                              fontsize=9, fontface="bold.italic"))))
prpdS <- ggplot(data = data.frame(y = residuals(drm_i_pdS, type = "pearson")[-c(1:3)],
                                  x = fitted(drm_i_pdS)[-c(1:3)]),
                aes(x = x, y = y)) +
  geom_hline(yintercept = 0) +
  geom_smooth(method = "loess", se = TRUE, formula = y ~ x, lwd = 0.5, col = mcol[2]) +
  geom_point(col = mcol[2], stroke = 0, size = 2) +
  scale_x_continuous(n.breaks = 3) + scale_y_continuous(n.breaks = 5) +
  xlab("Fitted") + ylab("Pearson residuals") +
  theme_bw() + #theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  annotation_custom(grobTree(textGrob(paste("RMSE =",round(RMSE(fit = fitted(drm_i_pdS), obs = df_i$R[df_i$S <= max_pre_decrease_S]), 0)),
                                      x=0.1,  y=0.9, hjust=0, gp=gpar(col=mcol[2], fontsize=9, fontface="bold.italic"))))

prINH <- ggplot(data = data.frame(y = residuals(drm_i_inh, type = "pearson")[-c(1:3)],
                                  x = fitted(drm_i_inh)[-c(1:3)]),
                aes(x = x, y = y)) + 
  geom_hline(yintercept = 0) +
  geom_smooth(method = "loess", se = TRUE, formula = y ~ x, lwd = 0.5, col = mcol[3]) +
  geom_point(col = mcol[3], stroke = 0, size = 2) +
  scale_x_continuous(n.breaks = 3) + scale_y_continuous(n.breaks = 5) +
  xlab("Fitted") + ylab("Pearson residuals") +
  theme_bw() + #theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  annotation_custom(grobTree(textGrob(paste("RMSE =",round(RMSE(fit = fitted(drm_i_inh), obs = df_i$R), 0)),
                                      x=0.1,  y=0.9, hjust=0,
                                      gp=gpar(col=viridis::plasma(n=5, alpha = 0.9)[4],
                                              fontsize=9, fontface="bold.italic"))))

# Plot data and model predictions
p0 <- ggplot(data = df_i, aes(x = S, y = R)) +
  #ggtitle(paste(Enzyme,df_i[1,"ID"])) +
  xlab(expression("Substrate concentration"~(µmol~L^-1))) +
  ylab(expression("Reaction rate"~(nmol~{g[dm]}^-1~h^-1))) +
  geom_point(alpha = 0.33, stroke = 0, col = "black", size = 2) +
  ylim(c(0, max(max(df_i$R), Vmax_i, coef(drm_i_pdS)[1]))) + xlim(c(0, maxSpred)) +
  theme_bw() #+ theme(axis.title.x = element_blank(), axis.title.y = element_blank())
if (boot) {
  p0 <- p0 +
    geom_line(data = mmbp, aes(group = Boot), col = viridis::plasma(n=5, alpha = 0.05)[2], lwd = 0.1) +
    geom_line(data = inhbp, aes(group = Boot), col = viridis::plasma(n=5, alpha = 0.05)[4], lwd = 0.1)
}
p1 <- p0 + 
  geom_line(data = pc, aes(col = Model), lwd = 0.75) +
  geom_line(data = pcx, aes(col = Model), lwd = 0.75, lty = 2) +
  scale_color_manual(values = mcol) +
  geom_vline(xintercept = max_pre_decrease_S, col = "darkgrey", lwd = 0.2) +
  geom_hline(yintercept = Vmax_i, col = mcol[1], lty = 4, lwd = 0.5) + # Vmax
  geom_line(data = data.frame(x = rep(Km_i, 2), y = c(0, Vmax_i/2)), aes(x = x, y = y), # Km
            col = mcol[1], lty = 3, lwd = 0.5) +
  geom_line(data = data.frame(x = c(0, Km_i), y = rep(Vmax_i/2, 2)), aes(x = x, y = y), # Vmax/2
            col = mcol[1], lty = 3, lwd = 0.5) +
  geom_hline(yintercept = coef(drm_i_pdS)[1], col = mcol[2], lty = 4, lwd = 0.5) + # Vmax - low S fit
  geom_line(data = data.frame(x = rep(coef(drm_i_pdS)[2], 2), y = c(0, coef(drm_i_pdS)[1]/2)), aes(x = x, y = y), # Km - low S fit
            col = mcol[2], lty = 3, lwd = 0.5) +
  geom_line(data = data.frame(x = c(0, coef(drm_i_pdS)[2]), y = rep(coef(drm_i_pdS)[1]/2, 2)), aes(x = x, y = y), # Vmax/2 - low S fit
            col = mcol[2], lty = 3, lwd = 0.5) +
  geom_point(data = df_i, aes(x = S, y = R), alpha = 0.33, stroke = 0, col = "black", size = 2)

p2 <- p1 + geom_line(data = data.frame(x = rep(SRmax_i, 2), y = c(0, Rmax_i)), aes(x = x, y = y), # SRmax
                     col = mcol[3], lty = 4, lwd = 0.5) +
  geom_hline(yintercept = Rmax_i, col = mcol[3], lty = 4, lwd = 0.5) + # Rmax
  # geom_line(data = data.frame(x = c(0, SRmax_i), y = rep(Rmax_i, 2)), aes(x = x, y = y), # Rmax
  #           col = mcol[3], lty = 2, lwd = 0.5) +
  geom_line(data = data.frame(x = rep(SRmax_i_05, 2), y = c(0, Rmax_i/2)), aes(x = x, y = y), # SRmax_i_05
            col = mcol[3], lty = 3, lwd = 0.5) +
  geom_line(data = data.frame(x = c(0, SRmax_i_05), y = rep(Rmax_i/2,2)), aes(x = x, y = y), # Rmax/2
            col = mcol[3], lty = 3, lwd = 0.5) +
  geom_point(data = df_i, aes(x = S, y = R), alpha = 0.33, stroke = 0, col = "black", size = 2)

((p2 + theme(legend.position = "none") | (prMM / prINH)) +
    plot_layout(width = c(2, 1)))
#}