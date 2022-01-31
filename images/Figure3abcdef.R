# Figure 3

library(ggplot2)
library(patchwork)

# Vmax, Km, Rmax, S(Rmax/2) relations ----

# Load and combine datasets to one long table:
params_BG <- cbind(read.csv(file = "./output/params_high_BG_nls.csv"), Model = "MM")
params_NAG <- cbind(read.csv(file = "./output/params_high_NAG_nls.csv"), Model = "MM")
params_AP <- cbind(read.csv(file = "./output/params_high_AP_nls.csv"), Model = "MM")

params_inh_BG <- cbind(read.csv(file = "./output/params_high_BG_nls_inh.csv"), Model = "INH")
params_inh_NAG <- cbind(read.csv(file = "./output/params_high_NAG_nls_inh.csv"), Model = "INH")
params_inh_AP <- cbind(read.csv(file = "./output/params_high_AP_nls_inh.csv"), Model = "INH")

# put subsets together:
params <- rbind(params_BG, params_NAG, params_AP)
params_inh <- rbind(params_inh_BG, params_inh_NAG, params_inh_AP)

params$max <- params$Vmax
params$half <- params$Km
params_inh$max <- params_inh$Rmax
params_inh$half <- params_inh$SRmax05
md <- rbind(params[,c("sample","Enzyme","Model","max","half")],
            params_inh[,c("sample","Enzyme","Model","max","half")])
md$Model <- factor(md$Model, levels = c("MM","INH"))
md$Enzyme <- factor(md$Enzyme, levels = c("BG","NAG","AP"))

sum(params[,c("sample")] - params_inh[,c("sample")])
pai <- params_inh
pai$inhVmax <- params_inh$Vmax
pai$inhKm <- params_inh$Km
pai$inhKi <- params_inh$Ki
pd <- cbind(params[,c("sample","Enzyme","Vmax","Km")],
            pai[,c("inhVmax","inhKm","inhKi","Rmax","SRmax","SRmax05","SRmax05highS","Sinhibition","INHbetter_aov_pvalue")])
pd$Enzyme <- factor(pd$Enzyme, levels = c("BG","NAG","AP"))

# INH model parameters


# Scatterplots for relation of parameters

RmaxVmax <- ggplot(data = pd, aes(x = Rmax, y = Vmax, col = Enzyme)) +
  geom_point(aes(shape = Sinhibition), alpha = 0.5, stroke = 0) +
  geom_abline(slope = 1, intercept = 0) +
  scale_y_log10() + scale_x_log10() + coord_fixed() +
  xlab(expression(R[max]~"(INH model)")) + ylab(expression(V[max]~"(MM model)")) +
  #facet_grid(~INHbetter_aov_pvalue<0.05) +
  scale_color_viridis_d(option = "D", begin = 0.1, end = 0.9) + theme_bw()
VmaxVmax <- ggplot(data = pd, aes(x = inhVmax, y = Vmax, col = Enzyme)) +
  geom_point(aes(shape = Sinhibition), alpha = 0.5, stroke = 0) +
  geom_abline(slope = 1, intercept = 0) +
  scale_y_log10() + scale_x_log10() + coord_fixed() +
  xlab(expression(V[max]~"(INH model)")) + ylab(expression(V[max]~"(MM model)")) +
  #facet_grid(~INHbetter_aov_pvalue<0.05) +
  scale_color_viridis_d(option = "D", begin = 0.1, end = 0.9) + theme_bw()
KmKm <- ggplot(data = pd, aes(x = inhKm, y = Km, col = Enzyme)) +
  geom_point(aes(shape = Sinhibition), alpha = 0.5, stroke = 0) +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed() + scale_y_log10() + scale_x_log10() + 
  xlab(expression(K[m]~"(INH model)")) + ylab(expression(K[m]~"(MM model)")) +
  #facet_grid(~INHbetter_aov_pvalue<0.05) +
  scale_color_viridis_d(option = "D", begin = 0.1, end = 0.9) + theme_bw()
SRmax05Km <- ggplot(data = pd, aes(x = SRmax05, y = Km, col = Enzyme)) +
  geom_point(aes(shape = Sinhibition), alpha = 0.5, stroke = 0) +
  #geom_smooth(method = "lm", alpha = 0.1, lwd = 0.2) +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed() + scale_y_log10() + scale_x_log10() + 
  xlab(expression(S(R[max]/2)~"(INH model)")) + ylab(expression(K[m]~"(MM model)")) +
  #facet_grid(~INHbetter_aov_pvalue<0.05) +
  scale_color_viridis_d(option = "D", begin = 0.1, end = 0.9) + theme_bw()
SRmax05Int_Ki <- ggplot(data = pd,#subset(pd, inhKi < 9e+06),
                        aes(y = ((SRmax05highS - SRmax) / SRmax05),
                      x = inhKi, col = Enzyme)) +
  geom_point(aes(shape = Sinhibition), alpha = 0.5, stroke = 0) +
  xlab(expression(K[i]~"(INH model)")) +
  ylab(expression((S(R[max]/2)[high] - S(R[max])) / S(R[max]/2)[low] ~"(INH model)")) +
  scale_x_log10() + scale_y_log10() +
  scale_color_viridis_d(option = "D", begin = 0.1, end = 0.9) + theme_bw()

# PRSE ----

# Load and combine datasets to one long table:
mpstats_BG <- cbind(read.csv(file = "./output/mpstats_high_BG_nls.csv"), PRSE_Ki = NA, Model = "MM")
mpstats_NAG <- cbind(read.csv(file = "./output/mpstats_high_NAG_nls.csv"), PRSE_Ki = NA, Model = "MM")
mpstats_AP <- cbind(read.csv(file = "./output/mpstats_high_AP_nls.csv"), PRSE_Ki = NA, Model = "MM")

mpstats_inh_BG <- cbind(read.csv(file = "./output/mpstats_high_BG_nls_inh.csv"), Model = "INH")
mpstats_inh_NAG <- cbind(read.csv(file = "./output/mpstats_high_NAG_nls_inh.csv"), Model = "INH")
mpstats_inh_AP <- cbind(read.csv(file = "./output/mpstats_high_AP_nls_inh.csv"), Model = "INH")

# put subsets together:
mpstats <- rbind(mpstats_BG, mpstats_NAG, mpstats_AP) # for AP, MM model fits are reduced to cases where the INH model also converged 
mpstats_inh <- rbind(mpstats_inh_BG, mpstats_inh_NAG, mpstats_inh_AP)

mp <- rbind(mpstats, mpstats_inh)
mp$Model <- factor(mp$Model, levels = c("MM","INH"))
mp$Enzyme <- factor(mp$Enzyme, levels = c("BG","NAG","AP"))

mcol <- c(viridis::plasma(n=5, alpha = 0.66)[2],
          viridis::plasma(n=5, alpha = 0.33)[2],
          viridis::plasma(n=5, alpha = 0.66)[4])

pPRSEmax_all <- ggplot(data = mp, #subset(mp, PRSEmax < 100),
                       aes(x = Enzyme, y = PRSEmax, fill = Model)) + 
  geom_boxplot(outlier.size = 0.25) + #scale_y_log10() +
  ylab(expression(PRSE[max]~"(%)")) +
  ylim(0, 100) +
  scale_fill_manual(values = mcol[c(1,3)]) +
  #facet_grid(INHbetter_aov_pvalue<0.05 ~ Sinhibition) +
  #facet_grid(~ Sinhibition, labeller = "label_both") +
  theme_bw() + theme(axis.title.x = element_blank())

# (VmaxVmax / KmKm) |
#   (RmaxVmax / SRmax05Km) |
#   (pPRSEmax_Sin / pPRSEmax_Sno) +
VmaxVmax + KmKm + pPRSEmax_all +
  RmaxVmax + SRmax05Km +
  SRmax05Int_Ki + 
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect", ncol = 3)
ggsave("images/Figure3abcdef.pdf", width = 8, height = 5)

# Report number of observations not included in illustrated y-ranges
tapply(mp$PRSEmax, INDEX = list(mp$Model, mp$Enzyme, mp$Sinhibition),
       FUN = function(x){sum(na.omit(x)>100)})
