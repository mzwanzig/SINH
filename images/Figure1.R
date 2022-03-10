library(readxl)

# Load data
dr <- read_xlsx("./data/pk1_hydrolases_sep.xlsx", sheet = "pk1_hydrolases_high") # ___________________CHECK !
dr <- dr[order(dr$Sub_conc),]

# PATCHWORK - PLOT
Enzyme <- c("BG","NAG","AP")[1] # define Enzyme
sampleID <- 80 # define sample ID
df <- dr[,c("Sample_ID", "Sub_conc", Enzyme)]; names(df)[1] <- "ID"; names(df)[2] <- "S" ; names(df)[3] <- "R"; df$ID <- as.factor(df$ID)
df_i <- subset(df, ID == sampleID)
source("plot_MMandSINH.R")
# singular gradient error? -> repeat previous line until it once works
p11a <- p2
p11b <- prMM
p11c <- prINH

Enzyme <- c("BG","NAG","AP")[1] # define Enzyme
sampleID <- 75 # define sample ID
df <- dr[,c("Sample_ID", "Sub_conc", Enzyme)]; names(df)[1] <- "ID"; names(df)[2] <- "S" ; names(df)[3] <- "R"; df$ID <- as.factor(df$ID)
df_i <- subset(df, ID == sampleID)
source("plot_MMandSINH.R")
p21a <- p2
p21b <- prMM
p21c <- prINH

Enzyme <- c("BG","NAG","AP")[2] # define Enzyme
sampleID <- 54 # define sample ID
df <- dr[,c("Sample_ID", "Sub_conc", Enzyme)]; names(df)[1] <- "ID"; names(df)[2] <- "S" ; names(df)[3] <- "R"; df$ID <- as.factor(df$ID)
df_i <- subset(df, ID == sampleID)
source("plot_MMandSINH.R")
p12a <- p2
p12b <- prMM
p12c <- prINH

Enzyme <- c("BG","NAG","AP")[2] # define Enzyme
sampleID <- 20 # define sample ID
df <- dr[,c("Sample_ID", "Sub_conc", Enzyme)]; names(df)[1] <- "ID"; names(df)[2] <- "S" ; names(df)[3] <- "R"; df$ID <- as.factor(df$ID)
df_i <- subset(df, ID == sampleID)
source("plot_MMandSINH.R")
p22a <- p2
p22b <- prMM
p22c <- prINH

Enzyme <- c("BG","NAG","AP")[3] # define Enzyme
sampleID <- 78 # define sample ID
df <- dr[,c("Sample_ID", "Sub_conc", Enzyme)]; names(df)[1] <- "ID"; names(df)[2] <- "S" ; names(df)[3] <- "R"; df$ID <- as.factor(df$ID)
df_i <- subset(df, ID == sampleID)
df_i$R <- df_i$R
source("plot_MMandSINH.R")
p13a <- p2
p13b <- prMM
p13c <- prINH

Enzyme <- c("BG","NAG","AP")[3] # define Enzyme
sampleID <- 114 # define sample ID
df <- dr[,c("Sample_ID", "Sub_conc", Enzyme)]; names(df)[1] <- "ID"; names(df)[2] <- "S" ; names(df)[3] <- "R"; df$ID <- as.factor(df$ID)
df_i <- subset(df, ID == sampleID)
df_i$R <- df_i$R
source("plot_MMandSINH.R")
p23a <- p2
p23b <- prMM
p23c <- prINH

layout <-'
AABDDEGGH
AACDDFGGI
JJKMMNPPQ
JJLMMOPPR
'

wrap_plots(A = p11a + theme(legend.position = "none") + ggtitle("a"), B = p11b, C = p11c,
           D = p12a + theme(legend.position = "none") + ggtitle("b"), E = p12b, F = p12c,
           G = p13a + theme(legend.position = "none") + ggtitle("c"), H = p13b, I = p13c,
           J = p21a + theme(legend.position = "none") + ggtitle("d"), K = p21b, L = p21c,
           M = p22a + theme(legend.position = "none") + ggtitle("e"), N = p22b, O = p22c,
           P = p23a + theme(legend.position = "none") + ggtitle("f"), Q = p23b, R = p23c,
           design = layout)

ggsave("images/Figure1.pdf", width = 14, height = 7)
