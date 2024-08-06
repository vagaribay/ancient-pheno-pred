# author: Valeria Añorve-Garibay
# Stabilizing Selection: Figures

library(ggplot2)
library(dplyr)
library(ggpubr)
options(scipen = 999)

fname = "tPheno_aPS"

# aPS correlation
STABS.corr = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/stabilizing_selection/metrics/", fname, ".cor", sep = ""), header = TRUE)
STABS.corr$generation = factor(STABS.corr$generation, levels = c(400, 300, 200, 100, 0))
STABS.corr$hsq = factor(STABS.corr$hsq, levels = c("total", "mid"))
STABS.corr$w = factor(STABS.corr$w)
STABS.corr$r.squared = (STABS.corr$corr)^2
STABS.corr = STABS.corr %>%
  mutate(hsq_math = factor(hsq,
                           levels = c("total", "mid"),
                           labels = c(expression(h^2 == 1), expression(h^2 == 0.5))))
STABS.corr.plt = ggplot(STABS.corr, aes(x = generation, y = r.squared, fill = w)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() +
  facet_wrap(. ~ hsq_math, labeller = label_parsed) +
  theme_linedraw() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  labs(title = "aPS accuracy", x = expression("generations before the present, " ~ tau), y = expression(italic(r)^2(Y[j], hat(Y)[j]), list(j == "")), fill = expression(italic(w))) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")) + ylim(0,1)

# aPS MSE
STABS.mse = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/stabilizing_selection/metrics/", fname, ".mse", sep = ""), header = TRUE)
STABS.mse$generation = factor(STABS.mse$generation, levels = c(400, 300, 200, 100, 0))
STABS.mse$hsq = factor(STABS.mse$hsq, levels = c("total", "mid"))
STABS.mse$w = factor(STABS.mse$w)
STABS.mse = STABS.mse %>%
  mutate(hsq_math = factor(hsq,
                           levels = c("total", "mid"),
                           labels = c(expression(h^2 == 1), expression(h^2 == 0.5))))
STABS.mse.plt = ggplot(STABS.mse, aes(x = generation, y = round(MSE,2), fill = w)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() +
  facet_wrap(. ~ hsq_math, labeller = label_parsed) +
  theme_linedraw() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  labs(title = "aPS MSE", x = expression("generations before the present, " ~ tau), y = expression(italic(MSE)(Y[j], hat(Y)[j]), list(j == "")), fill = expression(italic(w))) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")) + ylim(0,1.5)

Fig_3 = ggarrange(STABS.corr.plt, STABS.mse.plt, common.legend = TRUE, ncol = 1, legend = "left", labels = c("A", "B"))
#ggsave("/Users/valeriagby/desktop/ancient-pheno-prediction/figures/Figure_3.png", Fig_3, width = 8, height = 6.5, dpi = 300)

# lm for tpheno & aPS 
data = data.frame()
lm.1.slope = c()
lm.1.rsquared = c()
lm.5.slope = c()
lm.5.rsquared = c()
for (gen in c(100400, 100300, 100200, 100100, 100000)) {
  tpheno.1.tmp = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/stabilizing_selection/total/1/tpheno_", gen, "_33.txt", sep = ""), header = FALSE, col.names = "tpheno")
  aPS.1.tmp = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/stabilizing_selection/total/1/aPs_", gen, "_33.txt", sep = ""), header = FALSE, col.names = "aPS")
  
  tpheno.5.tmp = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/stabilizing_selection/total/5/tpheno_", gen, "_33.txt", sep = ""), header = FALSE, col.names = "tpheno")
  aPS.5.tmp = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/stabilizing_selection/total/5/aPs_", gen, "_33.txt", sep = ""), header = FALSE, col.names = "aPS")
  
  tmp.1 = cbind(tpheno.1.tmp, aPS.1.tmp)
  tmp.1$w = "1"
  tmp.5 = cbind(tpheno.5.tmp, aPS.5.tmp)
  tmp.5$w = "5"
  
  lm.tmp.1 = lm(aPS ~ tpheno, data = tmp.1)
  lm.1.slope = c(lm.1.slope, coef(lm.tmp.1)[2])
  lm.1.rsquared = c(lm.1.rsquared, summary(lm.tmp.1)$r.squared)
  
  lm.tmp.5 = lm(aPS ~ tpheno, data = tmp.5)
  lm.5.slope = c(lm.5.slope, coef(lm.tmp.5)[2])
  lm.5.rsquared = c(lm.5.rsquared, summary(lm.tmp.5)$r.squared)
  
  tmp = rbind(tmp.1, tmp.5)
  pgen = abs(gen - 100400)
  tmp$gen = pgen
  data = rbind(data, tmp)
}

data = data %>%
  mutate(w_math = factor(w,
                         levels = c(1, 5),
                         labels = c(expression(italic(w) == 1), expression(italic(w) == 5))))
lmTphenoAPS = ggplot(data, aes(x = tpheno, y = aPS)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red", lwd = 1, se = FALSE) +
  facet_grid(rows = vars(gen), cols = vars(w_math), scales="free_y", switch = 'y', labeller = label_parsed) +
  theme_linedraw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "true phenotype", y = "ancient polygenic score (aPS)") + ylim(-2.2,2.2)

ggsave("/Users/valeriagby/desktop/ancient-pheno-prediction/output/figures/Figure_3_Supp_1.png", lmTphenoAPS, width = 8, height = 5, dpi = 300)

# QTL: effect sizes
STABS.QTL = read.table("/Users/valeriagby/desktop/ancient-pheno-prediction/output/stabilizing_selection/metrics/QTLs_summ.txt", header = TRUE)
STABS.QTL$bin = factor(STABS.QTL$bin, levels = unique(STABS.QTL$bin))
STABS.QTL$hsq = factor(STABS.QTL$hsq, levels = c("total", "mid"))
STABS.QTL$w = factor(STABS.QTL$w)
STABS.QTL = STABS.QTL %>%
  mutate(hsq_math = factor(hsq,
                           levels = c("total", "mid"),
                           labels = c(expression(h^2 == 1), expression(h^2 == 0.5))))
STABS.QTL = STABS.QTL %>%
  mutate(w_math = factor(w,
                         levels = 1:5,
                         labels = c(expression(italic(w) == 1), expression(italic(w) == 2), expression(italic(w) == 3), expression(italic(w) == 4), expression(italic(w) == 5))))

STABS.QTL.plt = ggplot(STABS.QTL, aes(x = bin, y = nQTLs, fill = group)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() + 
  facet_grid(rows = vars(hsq_math), cols = vars(w_math), scales="free_y", switch = 'y',
             labeller = label_parsed) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "effect size (s)", y = "number of QTL mutations") +
  scale_fill_manual(values = c("#D55E00", "#56B4E9"), labels = c("Conserved", "Lost")) + ylim(0,150)

STABS.QTL.C = STABS.QTL[which(STABS.QTL$group == "conserved"),]
RATIO = ggplot(STABS.QTL.C, aes(x = bin, y = coeff, fill = w)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() +
  facet_grid(rows = vars(hsq_math), cols = vars(w_math), scales="free_y", switch = 'y',
             labeller = label_parsed) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "effect size (s)", y = "Conserved QTLs / (Conserved QTLs + Lost QTLs)") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7"))

Fig_3_2 = ggarrange(STABS.QTL.plt, RATIO, ncol = 1, align = "v", labels = c("A", "B"))
ggsave("/Users/valeriagby/desktop/ancient-pheno-prediction/figures/Figure_3_2.png", Fig_3_2, width = 15, height = 12, dpi = 300)

####################################################################################

library(ggplot2)
library(dplyr)
library(ggpubr)

options(scipen = 999)

# Stabilizing Selection with real parameters by Sanjak et al
fname = "tPheno_aPS"

# aPS correlation
STABS.corr = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/real_estimates/metrics/", fname, ".cor", sep = ""), header = TRUE)
STABS.corr$generation = factor(STABS.corr$generation, levels = c(400, 300, 200, 100, 0))
STABS.corr$pheno = factor(STABS.corr$hsq, levels = c("height", "bmi"))
STABS.corr$r.squared = (STABS.corr$corr)^2
corr.real = ggplot(STABS.corr, aes(x = generation, y = r.squared, fill = pheno)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  labs(title = "aPS accuracy", x = expression("generations before the present, " ~ tau), y = expression(italic(r)^2(Y[j], hat(Y)[j]), list(j == "")), fill = "") +
  scale_fill_manual(values = c("#DDCC77", "#44AA99"), labels = c("Height", "BMI")) + ylim(0,1)

# aPS MSE
STABS.mse = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/real_estimates/metrics/", fname, ".mse", sep = ""), header = TRUE)
STABS.mse$generation = factor(STABS.mse$generation, levels = c(400, 300, 200, 100, 0))
STABS.mse$pheno = factor(STABS.mse$hsq, levels = c("height", "bmi"))

mse.real = ggplot(STABS.mse, aes(x = generation, y = round(MSE,2), fill = pheno)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  labs(title = "aPS MSE", x = expression("generations before the present, " ~ tau), y = expression(italic(MSE)(Y[j], hat(Y)[j]), list(j == "")), fill = "") +
  scale_fill_manual(values = c("#DDCC77", "#44AA99"), labels = c("Height", "BMI")) + ylim(0,1.5)

Fig_6 = ggarrange(corr.real, mse.real, common.legend = TRUE, ncol = 2, legend = "left", labels = c("A", "B"))
ggsave("/Users/valeriagby/desktop/ancient-pheno-prediction/figures/Figure_6.png", real.plt, width = 9, height = 4, dpi = 300)
