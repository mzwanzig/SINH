# Figure 2

library(ggplot2)
library(patchwork)

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

# First row: model performance over the full S range
pMAPE_all <- ggplot(data = mp, aes(x = Enzyme, y = MAPE, fill = Model)) +
  geom_boxplot(outlier.size = 0.5) + scale_y_log10() + ylab("MAPE (%)") +
  scale_fill_manual(values = mcol[c(1,3)]) +
  theme_bw() + theme(axis.title.x = element_blank())
pbiasPERCENT_all <- ggplot(data = mp, aes(x = Enzyme, y = biasPERCENT, fill = Model)) +
  geom_hline(yintercept = 0) + geom_boxplot(outlier.size = 0.5) +
  ylab("Bias (%)") + ylim(-20, 50) +
  scale_fill_manual(values = mcol[c(1,3)]) +
  theme_bw() + theme(axis.title.x = element_blank())

# lower lines: model performance at low and high S levels only

mp2 <- subset(mp, Model == "MM")
mp2$biasPERCENTp1 <- mp2$biasPERCENTp1pdS
mp2$MAPEp1 <- mp2$MAPEp1pdS
mp2$Model <- "MM (low S fit)"

mp3 <- rbind(mp, mp2)
mp3$Model <- factor(mp3$Model, levels = c("MM", "MM (low S fit)", "INH"))

MAPEp1 <- ggplot(data = mp3, aes(x = Enzyme, y = MAPEp1, fill = Model)) +
  geom_boxplot(outlier.size = 0.5) +
  ylab("MAPE (%)") + ylim(0, 50) +
  scale_fill_manual(values = mcol) +
  theme_bw() + theme(axis.title.x = element_blank())
MAPEp2 <- ggplot(data = mp, aes(x = Enzyme, y = MAPEp2, fill = Model)) +
  geom_boxplot(outlier.size = 0.5) +
  ylab("MAPE (%)") + ylim(0, 50) +
  scale_fill_manual(values = mcol[c(1,3)]) +
  theme_bw() + theme(axis.title.x = element_blank())
biasPERCENTp1 <- ggplot(data = mp3, aes(x = Enzyme, y = biasPERCENTp1, fill = Model)) +
  geom_hline(yintercept = 0) + geom_boxplot(outlier.size = 0.5) + 
  ylab("Bias (%)") + ylim(-20, 50) +
  scale_fill_manual(values = mcol) +
  theme_bw() + theme(axis.title.x = element_blank())
biasPERCENTp2 <- ggplot(data = mp, aes(x = Enzyme, y = biasPERCENTp2, fill = Model)) +
  geom_hline(yintercept = 0) + geom_boxplot(outlier.size = 0.5) + 
  ylab("Bias (%)") + ylim(-20, 50) +
  scale_fill_manual(values = mcol[c(1,3)]) +
  theme_bw() + theme(axis.title.x = element_blank())

pMAPE_all + pbiasPERCENT_all +
  MAPEp1 + biasPERCENTp1 +
  MAPEp2 + biasPERCENTp2 +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect", ncol = 2)

ggsave("images/Figure2bcdefg.pdf", width = 6, height = 7)

# Report number of observations not included in illustrated y-ranges
tapply(mp$MAPE, INDEX = mp$Model, FUN = function(x){sum(na.omit(x)>50)})
tapply(mp$biasPERCENT, INDEX = mp$Model, FUN = function(x){sum(na.omit(x)>50)})
tapply(mp3$MAPEp1, INDEX = mp3$Model, FUN = function(x){sum(na.omit(x)>50)})
tapply(mp3$biasPERCENTp1, INDEX = mp3$Model, FUN = function(x){sum(na.omit(x)>50)})
tapply(mp3$biasPERCENTp1, INDEX = mp3$Model, FUN = function(x){sum(na.omit(x)<(-20))})
tapply(mp$MAPEp2, INDEX = mp$Model, FUN = function(x){sum(na.omit(x)>50)})
tapply(mp$biasPERCENTp2, INDEX = mp$Model, FUN = function(x){sum(na.omit(x)>50)})
tapply(mp$biasPERCENTp2, INDEX = mp$Model, FUN = function(x){sum(na.omit(x)<(-20))})


# Check condition of substrate inhibition vs. INH-MM model suitability
df <- subset(mp, Model == "INH") # avoids repetitions when using model comparison results only
tapply(df$INHbetter_aov_pvalue < 0.05, INDEX = df$Enzyme, FUN = sum)
tapply(df$INHbetter_aov_pvalue > 0, INDEX = df$Enzyme, FUN = sum)

tapply(df$INHbetter_aov_pvalue < 0.05, INDEX = list(df$Enzyme, df$Sinhibition), FUN = sum)
tapply(df$INHbetter_aov_pvalue > 0, INDEX = list(df$Enzyme, df$Sinhibition), FUN = sum)

t0 <- tapply(df$INHbetter_aov_pvalue > 0, INDEX = list(df$Enzyme, df$Sinhibition), FUN = sum)
Sinh_df <- data.frame(Count = c(t0[,1], t0[,2]),
                      Enzyme = rep(dimnames(t0)[[1]],2),
                      SInhibition = rep(dimnames(t0)[[2]], each = 3))
Sinh_df$Enzyme <- factor(Sinh_df$Enzyme, levels = c("BG","NAG","AP"))

t1 <- tapply(df$INHbetter_aov_pvalue < 0.05, INDEX = list(df$Enzyme, df$Sinhibition), FUN = sum)
INHb_df <- data.frame(Count = c(t1[,1], t1[,2]),
                      Enzyme = rep(dimnames(t1)[[1]],2),
                      SInhibition = rep(dimnames(t1)[[2]], each = 3))
INHb_df$Enzyme <- factor(INHb_df$Enzyme, levels = c("BG","NAG","AP"))

ggplot() +
  geom_col(data = Sinh_df, aes(y = Count, x = Enzyme),
           fill =  viridis::plasma(n=5)[2], position = 'dodge') +
  geom_col(data = INHb_df, aes(y = Count, x = Enzyme),
           fill =  viridis::plasma(n=5)[4], position = 'dodge') +
  facet_grid(~SInhibition) +
  theme_bw()
ggsave("images/Figure2a.pdf", width = 2.5, height = 2)
