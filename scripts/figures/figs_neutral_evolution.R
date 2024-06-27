# author: Valeria Añorve-Garibay
# Neutral Evolution: Figures

library(ggplot2)
library(dplyr)
library(ggpubr)
options(scipen = 999)

fname = "tPheno_aPS"

# aPS accuracy
NULL.corr = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/neutral_evolution/metrics/", fname, ".cor", sep = ""), header = TRUE)
NULL.corr$generation = factor(NULL.corr$generation, levels = c(400, 300, 200, 100, 0))
NULL.corr$hsq = factor(NULL.corr$hsq, levels = c("total", "mid"))
NULL.corr$r.squared = (NULL.corr$corr)^2
NULL.corr.plt = ggplot(NULL.corr, aes(x = generation, y = r.squared, fill = hsq)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8, position = position_dodge(0.8)) +
  theme_linedraw() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12)) +
  labs(title = "aPS accuracy", x = expression("generations before the present, " ~ tau), y = expression(italic(r)^2(Y[j], hat(Y)[j]), list(j == "")), fill = expression(h^2)) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442"), labels = c("1.0", "0.5")) + ylim(0,1)

# aPS MSE
NULL.mse = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/neutral_evolution/metrics/", fname, ".mse", sep = ""), header = TRUE)
NULL.mse$generation = factor(NULL.mse$generation, levels = c(400, 300, 200, 100, 0))
NULL.mse$hsq = factor(NULL.mse$hsq, levels = c("total", "mid"))
NULL.mse.plt = ggplot(NULL.mse, aes(x = generation, y = MSE, fill = hsq)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8, position = position_dodge(0.8)) +
  theme_linedraw() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12)) + 
  labs(title = "aPS MSE", x = expression("generations before the present, " ~ tau), y = expression(italic(MSE)(Y[j], hat(Y)[j]), list(j == "")), fill = expression(h^2)) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442"), labels = c("1.0", "0.5")) + ylim(0,30)

Fig_2 = ggarrange(NULL.corr.plt, NULL.mse.plt, common.legend = TRUE, ncol = 2, legend = "left", labels = c("A", "B"))
ggsave("/Users/valeriagby/desktop/ancient-pheno-prediction/figures/Figure_2.png", Fig_2, width = 9, height = 4, dpi = 300)

# QTL: effect sizes
fname = "tPheno_aPS"
NULL.QTL = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/neutral_evolution/metrics/", fname, "_0.8-0.8-0.1.QTLs", sep = ""), header = TRUE)
NULL.QTL$bin = factor(NULL.QTL$bin, levels = unique(NULL.QTL$bin))
NULL.QTL$hsq = factor(NULL.QTL$hsq, levels = c("total", "mid"))
NULL.QTL.plt = NULL.QTL %>%
  mutate(hsq_math = factor(hsq,
                           levels = c("total", "mid"),
                           labels = c(expression(h^2 == 1), expression(h^2 == 0.5))))

NULL.QTL.EZ = ggplot(NULL.QTL.plt, aes(x = bin, y = nQTLs, fill = group)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() + 
  facet_wrap(~hsq_math, labeller = label_parsed) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "effect size (s)", y = "number of QTL mutations") +
  scale_fill_manual(values = c("#D55E00", "#56B4E9"), labels = c("Conserved", "Lost"))


NULL.QTL.plt.C = NULL.QTL.plt[which(NULL.QTL.plt$group == "conserved"),]
RATIO = ggplot(NULL.QTL.plt.C, aes(x = bin, y = coeff, fill = hsq)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() +
  facet_wrap(~hsq_math, labeller = label_parsed) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "effect size (s)",  y = "Conserved QTLs / Conserved QTLs + Lost QTLs", fill = expression(h^2)) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442"), labels = c("1.0", "0.5"))
  
Fig_2_1 = ggarrange(NULL.QTL.EZ, RATIO, ncol = 1, align = "v", labels = c("A", "B"))
ggsave("/Users/valeriagby/desktop/ancient-pheno-prediction/figures/Figure_2_1.png", Fig_2_1, width = 12, height = 10, dpi = 300)

# QTL: frequencies
fname = "tPheno_aPS"
NULL.QTL.frq = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/neutral_evolution/metrics/", fname, "_0-1-0.001.QTLs.frq", sep = ""), header = TRUE)
NULL.QTL.frq.L = NULL.QTL.frq[which(NULL.QTL.frq$hsq == "total" & NULL.QTL.frq$group == "lost"),]
NULL.QTL.frq.L$bin = factor(NULL.QTL.frq.L$bin, levels = unique(NULL.QTL.frq.L$bin))
percent_bin = c()
max_n = c()
for (bin in levels(NULL.QTL.frq.L$bin)) {
  percent_bin = c(percent_bin, (sum(NULL.QTL.frq.L[which(NULL.QTL.frq.L$bin == bin),"nQTLs"]) / sum(NULL.QTL.frq.L$nQTLs)) * 100)
  max_n = c(max_n, max(NULL.QTL.frq.L[which(NULL.QTL.frq.L$bin == bin),"nQTLs"]))
}
annot_text = data.frame(bin = levels(NULL.QTL.frq.L$bin)[0:20], label = cumsum(percent_bin)[0:20], y = max_n[0:20])
NULL.QTL.frq.L.tmp = NULL.QTL.frq.L[1:2000,]
NULL.QTL.frq.L.tmp$bin = factor(NULL.QTL.frq.L.tmp$bin, levels = unique(NULL.QTL.frq.L.tmp$bin))
frq.plt = ggplot(NULL.QTL.frq.L.tmp, aes(x = bin, y = nQTLs)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8, fill = "#56B4E9") +
  geom_hline(yintercept = 75, linetype = "dashed", color = "gray70") +
  theme_linedraw() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "frequency", y = "number of QTL mutations") + ylim(0,100) +
  geom_text(data = data.frame(), aes(x = annot_text$bin , y = 78, 
                                     label = paste(round(annot_text$label, 1), "%", sep = "")), size = 3.2)

# DTWF: prob of transitioning from f = X to f = 0
DTWF = read.table("/Users/valeriagby/desktop/ancient-pheno-prediction/scripts/fastDTWF/DTWF_seq_0_02_OUT.txt", header = FALSE)
probs = seq(0, 20000, 1) / 20000
colnames(DTWF) = c("freq", "generation", paste("fg", probs, sep = ""))
DTWF$rowSums = rowSums(DTWF[,-which(names(DTWF) %in% c("freq", "generation"))])
DTWF$probSeg = rowSums(DTWF[,-which(names(DTWF) %in% c("freq", "generation", "fg0", "rowSums"))])

# probability of transitioning from f% to 0%
# probability of lossing allele A
DTWF$generation = factor(DTWF$generation)
DTWF.plt = ggplot(DTWF, aes(x = freq, y = fg0, color = generation)) +
  geom_line() + geom_point() +
  theme_linedraw() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "", x = "frequency", main = "probability of transitioning from f freq to f = 0 freq", color = "", y = "probability") +
  scale_color_manual(values = c("#CC6677", "#DDCC77", "#6699CC", "#009E73"), labels = c(100, 200, 300, 400)) + ylim(0,1) +
  scale_x_continuous(breaks = seq(0,0.02,0.001))

Fig_2_2 = ggarrange(frq.plt, DTWF.plt, ncol = 1, labels = c("A", "B"), align = "hv")
ggsave("/Users/valeriagby/desktop/ancient-pheno-prediction/figures/Figure_2_2.png", Fig_2_2, width = 10, height = 8, dpi = 300)


################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################
# Ratio of Conserved QTLs / Total QTLs: all evolutionary dynamics
# Neutral Evolution
fname = "tPheno_aPS"
NULL.QTL = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/neutral_evolution/metrics/", fname, "_0.8-0.8-0.1.QTLs", sep = ""), header = TRUE)
NULL.QTL = NULL.QTL[which(NULL.QTL$hsq == "total"),]
NULL.QTL$bin = factor(NULL.QTL$bin, levels = unique(NULL.QTL$bin))
NULL.QTL.C = NULL.QTL[which(NULL.QTL$group == "conserved"),]
NULL.coeff.plt = ggplot(NULL.QTL.C, aes(x = bin, y = coeff)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "effect size (s)",  y = "Conserved QTLs / Conserved QTLs + Lost QTLs", title = "Neutral Evolution")

# Stabilizing Selection
# Effect Sizes
STABS.QTL = read.table("/Users/valeriagby/desktop/ancient-pheno-prediction/output/stabilizing_selection/metrics/QTLs_summ.txt", header = TRUE)
STABS.QTL = STABS.QTL[which(STABS.QTL$hsq == "total" & STABS.QTL$w == 1),]
STABS.QTL$bin = factor(STABS.QTL$bin, levels = unique(STABS.QTL$bin))
STABS.QTL.C = STABS.QTL[which(STABS.QTL$group == "conserved"),]
STAB.coeff.plt = ggplot(STABS.QTL.C, aes(x = bin, y = coeff)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "effect size (s)", y = "Conserved QTLs / Conserved QTLs + Lost QTLs", title = "Stabilizing Selection")

# Directional Selection
DIR.QTL.0025 = data.frame()
for (genbp in c(400,300,200,100,0)) {
  tmp.0025 = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/directional_selection/metrics/QTLs_summ_0025_0_", genbp, ".txt", sep = ""), header = TRUE, colClasses = c("integer", "character", "integer", "character", "character", "integer", "character", "numeric"))
  tmp.0025$gbp = genbp
  DIR.QTL.0025 = rbind(DIR.QTL.0025, tmp.0025)
}
DIR.QTL.0025 = DIR.QTL.0025[which(DIR.QTL.0025$hsq == "total" & DIR.QTL.0025$gbp == "400"),]
DIR.QTL.0025$bin = factor(DIR.QTL.0025$bin, levels = unique(DIR.QTL.0025$bin))
DIR.QTL.0025.C = DIR.QTL.0025[which(DIR.QTL.0025$group == "conserved"),]
DIR.coeff.plt = ggplot(DIR.QTL.0025.C, aes(x = bin, y = coeff)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "effect size (s)", y = "Conserved QTLs / Conserved QTLs + Lost QTLs", title = "Directional Selection")

Fig_7 = ggarrange(NULL.coeff.plt, STAB.coeff.plt, DIR.coeff.plt, ncol = 3,
          align = "hv", labels = c("A", "B", "C"))
ggsave("/Users/valeriagby/desktop/ancient-pheno-prediction/figures/Figure_7.png", Fig_7, width = 15, height = 6, dpi = 300)





