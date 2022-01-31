# Exemplary model fitting and visualization
# with bootstrapping and residual plots

# The data frame has to contain
# R (velocity of enzymatic reaction)
# S (substrate concentrations):
df_i <- data.frame(S = rep(c(0,5,10,25,50,100,500,1000), 3),
                   R = c(0, 174, 221, 520, 955, 1588, 3366, 3243,
                         0, 110, 199, 496, 934, 1336, 2808, 2618,
                         0, 157, 276, 553, 871, 1438, 2847, 2516))
df_i

# Fit the 2-parametric model
mMM <- nls(R ~ (Vmax*S)/(Km+S), # Michaelis-Menten function
           data = df_i,
           start = c(Vmax = max(df_i$R), Km = max(df_i$S)/3),
           control = nls.control(maxiter = 50, tol = 1e-04))
mMM # MM parameters and stats

# Fit the 3-parametric model
mSINH <- nls(R ~ Vmax * S / (Km + S *(1+S/Ki)), # Substrate inhibition model SINH
             data = df_i,
             start=list(Vmax = max(df_i$R), # initial guess
                        Km = max(df_i$S)/3, # initial guess
                        Ki = max(df_i$S)), # initial guess
             control = nls.control(maxiter = 50, tol = 1e-04))
summary(mSINH) # SINH parameters and stats
# functions to calculate SINH tipping points
# maximal enzyme reaction rate: Rmax
Rmax <- function(Vmax, Km, Ki) { data.frame(Rmax = Vmax / (1+2*sqrt(Km/Ki)), row.names = NULL) }
Rmax_i <- Rmax(Vmax = coef(mSINH)["Vmax"], Km = coef(mSINH)["Km"], Ki = coef(mSINH)["Ki"]); Rmax_i
# substrate concentration for maximum enzyme reaction rate: S(Rmax)
SRmax <- function(Km, Ki){ data.frame(SRmax = sqrt(Km * Ki), row.names = NULL) }
SRmax_i <- SRmax(Km = coef(mSINH)["Km"], Ki = coef(mSINH)["Ki"]); SRmax_i
# substrate concentration for half maximum reaction rate: S(Rmax/2)
Sseq <- seq(0, as.numeric(SRmax_i), length.out = 10000)
SRmax_i_05 <- Sseq[match(TRUE, predict(mSINH, newdata = data.frame(S = Sseq)) > as.numeric(Rmax_i)/2)]; SRmax_i_05

# Compares both models by residual sum of squares
anova(mMM, mSINH)

# Visualize model fits including bootstrapped predictions,
# root mean square error (RMSE) and residual plots:
source("plot_MMandSINH.R")
(p2 | (prMM / prINH)) +
  plot_layout(width = c(2, 1), guides = "collect")
ggsave("example.pdf", width = 7, height = 4)
