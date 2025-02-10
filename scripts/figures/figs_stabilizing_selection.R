# author: Valeria Añorve-Garibay
# Stabilizing Selection: Figures

library(ggplot2)
library(dplyr)
library(ggpubr)
options(scipen = 999)

# aPS correlation
STABS.corr = read.table('/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/metrics/tPheno_aPS.cor', header = TRUE)
STABS.corr$generation = factor(STABS.corr$generation, levels = c(200000, 200100, 200200, 200300, 200400))
STABS.corr$hsq = factor(STABS.corr$hsq, levels = c("total", "mid"))
STABS.corr$w = factor(STABS.corr$w)
STABS.corr$r.squared = (STABS.corr$corr)^2
STABS.corr = STABS.corr %>%
  mutate(hsq_math = factor(hsq,
                           levels = c("total", "mid"),
                           labels = c(expression(h^2 == 1), expression(h^2 == 0.5))))

STABS.corr.100 = STABS.corr[which(STABS.corr$present_sample_size == 100),]

#median(STABS.corr.100[which(STABS.corr.100$hsq == "total" & STABS.corr.100$generation == 0),"r.squared"])
#median(STABS.corr.100[which(STABS.corr.100$hsq == "total" & STABS.corr.100$generation == 100 & STABS.corr.100$w == 1),"r.squared"])
#median(STABS.corr.100[which(STABS.corr.100$hsq == "total" & STABS.corr.100$generation == 200 & STABS.corr.100$w == 1),"r.squared"])
#median(STABS.corr.100[which(STABS.corr.100$hsq == "total" & STABS.corr.100$generation == 300 & STABS.corr.100$w == 1),"r.squared"])
#median(STABS.corr.100[which(STABS.corr.100$hsq == "total" & STABS.corr.100$generation == 400 & STABS.corr.100$w == 1),"r.squared"])

#median(STABS.corr.100[which(STABS.corr.100$hsq == "mid" & STABS.corr.100$generation == 0 & STABS.corr.100$w == 1),"r.squared"])
#median(STABS.corr.100[which(STABS.corr.100$hsq == "mid" & STABS.corr.100$generation == 100 & STABS.corr.100$w == 1),"r.squared"])
#median(STABS.corr.100[which(STABS.corr.100$hsq == "mid" & STABS.corr.100$generation == 200 & STABS.corr.100$w == 1),"r.squared"])
#median(STABS.corr.100[which(STABS.corr.100$hsq == "mid" & STABS.corr.100$generation == 300 & STABS.corr.100$w == 1),"r.squared"])
#median(STABS.corr.100[which(STABS.corr.100$hsq == "mid" & STABS.corr.100$generation == 400 & STABS.corr.100$w == 1),"r.squared"])


STABS.corr.plt = ggplot(STABS.corr.100, aes(x = generation, y = r.squared, fill = w)) +
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
  labs(title = "aPGS accuracy", x = expression("generations before the present, " ~ tau), y = expression(italic(r)^2(Y[i], hat(Y)[i]), list(j == "")), fill = expression(italic(w))) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")) + ylim(0,1) +
  scale_x_discrete(breaks = c(200000, 200100, 200200, 200300, 200400), labels = c(400, 300, 200, 100, 0))

# aPS MSE
STABS.mse = read.table('/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/metrics/tPheno_aPS_scaled.mse', header = TRUE)
STABS.mse$generation = factor(STABS.mse$generation, levels = c(200000, 200100, 200200, 200300, 200400))
STABS.mse$hsq = factor(STABS.mse$hsq, levels = c("total", "mid"))
STABS.mse$w = factor(STABS.mse$w)
STABS.mse = STABS.mse %>%
  mutate(hsq_math = factor(hsq,
                           levels = c("total", "mid"),
                           labels = c(expression(h^2 == 1), expression(h^2 == 0.5))))


STABS.mse.100 = STABS.mse[which(STABS.mse$present_sample_size == 100),]
STABS.mse.plt = ggplot(STABS.mse.100, aes(x = generation, y = round(MSE,2), fill = w)) +
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
  labs(title = "aPGS MSE", x = expression("generations before the present, " ~ tau), y = expression(italic(MSE)(Y[i], hat(Y)[i]), list(j == "")), fill = expression(italic(w))) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")) + ylim(0,3.5) +
  scale_x_discrete(breaks = c(200000, 200100, 200200, 200300, 200400), labels = c(400, 300, 200, 100, 0))

# aPS bias
STABS.bias = read.table('/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/metrics/tPheno_aPS.bias', header = TRUE)
STABS.bias$generation = factor(STABS.bias$generation, levels = c(200000, 200100, 200200, 200300, 200400))
STABS.bias$hsq = factor(STABS.bias$hsq, levels = c("total", "mid"))
STABS.bias$w = factor(STABS.bias$w)
STABS.bias = STABS.bias %>%
  mutate(hsq_math = factor(hsq,
                           levels = c("total", "mid"),
                           labels = c(expression(h^2 == 1), expression(h^2 == 0.5))))
STABS.bias.100 = STABS.bias[which(STABS.bias$present_sample_size == 100),]
STABS.bias.plt = ggplot(STABS.bias.100, aes(x = generation, y = bias, fill = w)) +
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
  labs(title = "aPGS Bias", x = expression("generations before the present, " ~ tau), y = expression(italic(Bias)(Y[i], hat(Y)[i]), list(j == "")), fill = expression(italic(w))) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")) +  ylim(-1,1) +
  scale_x_discrete(breaks = c(200000, 200100, 200200, 200300, 200400), labels = c(400, 300, 200, 100, 0))
Fig_3 = ggarrange(STABS.corr.plt, STABS.mse.plt, STABS.bias.plt, common.legend = TRUE, ncol = 1, legend = "left", labels = c("A", "B", "C"))
ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_3.png", Fig_3, width = 8, height = 8, dpi = 300)


# lm for tpheno & aPS 
data = data.frame()
lm.1.slope = c()
lm.1.rsquared = c()
lm.5.slope = c()
lm.5.rsquared = c()
for (gen in c(200400, 200300, 200200, 200100, 200000)) {
  tpheno.1.tmp = read.table(paste("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/total/1/100_tpheno_", gen, "_33.txt", sep = ""), header = FALSE, col.names = "tpheno")
  aPS.1.tmp = read.table(paste("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/total/1/100_aPs_", gen, "_33.txt", sep = ""), header = FALSE, col.names = "aPS")
  
  tpheno.5.tmp = read.table(paste("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/total/5/100_tpheno_", gen, "_33.txt", sep = ""), header = FALSE, col.names = "tpheno")
  aPS.5.tmp = read.table(paste("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/total/5/100_aPs_", gen, "_33.txt", sep = ""), header = FALSE, col.names = "aPS")
  
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
  tmp$gen = gen
  data = rbind(data, tmp)
}

data = data %>%
  mutate(w_math = factor(w,
                         levels = c(1, 5),
                         labels = c(expression(italic(w) == 1), expression(italic(w) == 5))))

data$pgen = abs(data$gen - 200400)
data$pgen = factor(data$pgen)
lmTphenoAPS = ggplot(data, aes(x = tpheno, y = aPS)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red", lwd = 1, se = FALSE) +
  facet_grid(rows = vars(pgen), cols = vars(w_math), scales="free_y", switch = 'y', labeller = label_parsed) +
  theme_linedraw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "true phenotype", y = "ancient polygenic score (aPS)") + ylim(-3.5,3.5)
ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_3_1.png", lmTphenoAPS, width = 8, height = 5, dpi = 300)

regression_stats = data %>%
  group_by(gen, w_math) %>%
  summarise(
    slope = coef(lm(aPS ~ tpheno, data = cur_data()))[2],
    r_squared = summary(lm(aPS ~ tpheno, data = cur_data()))$r.squared,
    .groups = "drop"
  )

regression_stats$slope = round(regression_stats$slope, 2)
regression_stats$r_squared = round(regression_stats$r_squared, 2)

# QTL: effect sizes
STABS.QTL = read.table("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/metrics/QTLs_summ.txt", header = TRUE)
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
  scale_fill_manual(values = c("#D55E00", "#56B4E9"), labels = c("Conserved", "Lost"))


STABS.QTL.C = STABS.QTL[which(STABS.QTL$group == "conserved"),]
STABS.QTL.C[is.na(STABS.QTL.C$coeff),"coeff"] = 0
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
ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_3_2.png", Fig_3_2, width = 15, height = 12, dpi = 300)

# metrics w changing sample size in present
STABS.corr.mid = STABS.corr[which(STABS.corr$hsq == "mid" & (STABS.corr$w == 1 | STABS.corr$w == 5)),]
STABS.corr.mid$present_sample_size = as.factor(STABS.corr.mid$present_sample_size)
STABS.corr.mid = STABS.corr.mid %>%
  mutate(hsq_math = factor(w,
                           levels = c("1", "5"),
                           labels = c(expression(w == 1), expression(w == 5))))
STABS.corr.plt = ggplot(STABS.corr.mid, aes(x = generation, y = r.squared, fill = present_sample_size)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() +
  facet_wrap(. ~ hsq_math, labeller = label_parsed)  + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  labs(title = "aPGS accuracy", x = expression("generations before the present, " ~ tau), y = expression(italic(r)^2(Y[i], hat(Y)[i]), list(j == "")), fill = expression("sample size at" ~ tau == 0)) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442", "#E69F00")) + ylim(0,1) +
  scale_x_discrete(breaks = c(200000, 200100, 200200, 200300, 200400), labels = c(400, 300, 200, 100, 0))


STABS.mse.mid = STABS.mse[which(STABS.mse$hsq == "mid" & (STABS.mse$w == 1 | STABS.mse$w == 5)),]
STABS.mse.mid$present_sample_size = as.factor(STABS.mse.mid$present_sample_size)
STABS.mse.mid = STABS.mse.mid %>%
  mutate(hsq_math = factor(w,
                           levels = c("1", "5"),
                           labels = c(expression(w == 1), expression(w == 5))))

STABS.mse.plt = ggplot(STABS.mse.mid, aes(x = generation, y = round(MSE,2), fill = present_sample_size)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() +
  facet_wrap(. ~ hsq_math, labeller = label_parsed) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  labs(title = "aPGS MSE", x = expression("generations before the present, " ~ tau), y = expression(italic(MSE)(Y[i], hat(Y)[i]), list(j == "")), fill = expression("sample size at" ~ tau == 0)) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442", "#E69F00")) + ylim(0,3) +
  scale_x_discrete(breaks = c(200000, 200100, 200200, 200300, 200400), labels = c(400, 300, 200, 100, 0))


STABS.bias.mid = STABS.bias[which(STABS.bias$hsq == "mid" & (STABS.bias$w == 1 | STABS.bias$w == 5)),]
STABS.bias.mid$present_sample_size = as.factor(STABS.bias.mid$present_sample_size)
STABS.bias.mid = STABS.bias.mid %>%
  mutate(hsq_math = factor(w,
                           levels = c("1", "5"),
                           labels = c(expression(w == 1), expression(w == 5))))

STABS.bias.plt = ggplot(STABS.bias.mid, aes(x = generation, y = bias, fill = present_sample_size)) +
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
  labs(title = "aPGS Bias", x = expression("generations before the present, " ~ tau), y = expression(italic(Bias)(Y[i], hat(Y)[i]), list(j == "")), fill = expression("sample size at" ~ tau == 0)) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442", "#E69F00")) +
  scale_x_discrete(breaks = c(200000, 200100, 200200, 200300, 200400), labels = c(400, 300, 200, 100, 0))

Fig_3_3 = ggarrange(STABS.corr.plt, STABS.mse.plt, STABS.bias.plt, common.legend = TRUE, ncol = 1, legend = "left", labels = c("A", "B", "C"))
#ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_3_3.png", Fig_3_3, width = 8, height = 10, dpi = 300)


#################################################################################
#################################################################################
# Sims with bigger Genome Size
# r^2
# 500 kB / 20 regions
STABS.corr = read.table('/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/metrics/tPheno_aPS.cor', header = TRUE)
STABS.corr.mid.100 = STABS.corr[which(STABS.corr$hsq == "mid" & STABS.corr$present_sample_size == 100 & (STABS.corr$w == 1 | STABS.corr$w == 5)),]
STABS.corr.mid.100$genome_size = '20'
# 1.25 Mb / 50 regions
STABS.corr.big = read.table('/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/sims_w_bigger_genome/metrics/tPheno_aPS.cor', header = TRUE)
STABS.corr.big$genome_size = '50'
# merge
STABS.corr.gsize = rbind(STABS.corr.mid.100, STABS.corr.big)
STABS.corr.gsize = STABS.corr.gsize %>%
  mutate(hsq_math = factor(w,
                           levels = c("1", "5"),
                           labels = c(expression(w == 1), expression(w == 5))))
STABS.corr.gsize$generation = factor(STABS.corr.gsize$generation, levels = c(200000, 200100, 200200, 200300, 200400))
STABS.corr.gsize$r.squared = (STABS.corr.gsize$corr)^2
STABS.corr.plt = ggplot(STABS.corr.gsize, aes(x = generation, y = r.squared, fill = genome_size)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() +
  facet_wrap(. ~ hsq_math, labeller = label_parsed)  + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  labs(title = "aPGS accuracy", x = expression("generations before the present, " ~ tau), y = expression(italic(r)^2(Y[i], hat(Y)[i]), list(j == "")), fill = "genome size") +
  scale_fill_manual(values = c("#CC79A7", "#F0E442"), labels = c("500 kb", "1.25 Mb")) + ylim(0,1) +
  scale_x_discrete(breaks = c(200000, 200100, 200200, 200300, 200400), labels = c(400, 300, 200, 100, 0))
# MSE
# 500 kB / 20 regions
STABS.MSE = read.table('/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/metrics/tPheno_aPS_scaled.mse', header = TRUE)
STABS.MSE.mid.100 = STABS.MSE[which(STABS.MSE$hsq == "mid" & STABS.MSE$present_sample_size == 100 & (STABS.MSE$w == 1 | STABS.MSE$w == 5)),]
STABS.MSE.mid.100$genome_size = '20'
# 1.25 Mb / 50 regions
STABS.MSE.big = read.table('/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/sims_w_bigger_genome/metrics/tPheno_aPS_scaled.mse', header = TRUE)
STABS.MSE.big$genome_size = '50'
# merge
STABS.MSE.gsize = rbind(STABS.MSE.mid.100, STABS.MSE.big)
STABS.MSE.gsize = STABS.MSE.gsize %>%
  mutate(hsq_math = factor(w,
                           levels = c("1", "5"),
                           labels = c(expression(w == 1), expression(w == 5))))
STABS.MSE.gsize$generation = factor(STABS.MSE.gsize$generation, levels = c(200000, 200100, 200200, 200300, 200400))
STABS.mse.plt = ggplot(STABS.MSE.gsize, aes(x = generation, y = MSE, fill = genome_size)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() +
  facet_wrap(. ~ hsq_math, labeller = label_parsed)  + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14)) + 
  labs(title = "aPGS MSE", x = expression("generations before the present, " ~ tau), y = expression(italic(MSE)(Y[i], hat(Y)[i]), list(j == "")), fill = "genome size") +
  scale_fill_manual(values = c("#CC79A7", "#F0E442"), labels = c("500 kb", "1.25 Mb")) + ylim(0,2) +
  scale_x_discrete(breaks = c(200000, 200100, 200200, 200300, 200400), labels = c(400, 300, 200, 100, 0))
# Bias
# 500 kB / 20 regions
STABS.bias = read.table('/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/metrics/tPheno_aPS.bias', header = TRUE)
STABS.bias.mid.100 = STABS.bias[which(STABS.bias$hsq == "mid" & STABS.bias$present_sample_size == 100 & (STABS.bias$w == 1 | STABS.bias$w == 5)),]
STABS.bias.mid.100$genome_size = '20'
# 1.25 Mb / 50 regions
STABS.bias.big = read.table('/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/sims_w_bigger_genome/metrics/tPheno_aPS.bias', header = TRUE)
STABS.bias.big$genome_size = '50'
# merge
STABS.bias.gsize = rbind(STABS.bias.mid.100, STABS.bias.big)
STABS.bias.gsize = STABS.bias.gsize %>%
  mutate(hsq_math = factor(w,
                           levels = c("1", "5"),
                           labels = c(expression(w == 1), expression(w == 5))))
STABS.bias.gsize$generation = factor(STABS.bias.gsize$generation, levels = c(200000, 200100, 200200, 200300, 200400))
STABS.bias.plt = ggplot(STABS.bias.gsize, aes(x = generation, y = bias, fill = genome_size)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() +
  facet_wrap(. ~ hsq_math, labeller = label_parsed)  + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14)) + 
  labs(title = "aPGS Bias", x = expression("generations before the present, " ~ tau), y = expression(italic(Bias)(Y[i], hat(Y)[i]), list(j == "")), fill = "genome size") +
  scale_fill_manual(values = c("#CC79A7", "#F0E442"), labels = c("500 kb", "1.25 Mb")) + ylim(-1.5,1.5) +
  scale_x_discrete(breaks = c(200000, 200100, 200200, 200300, 200400), labels = c(400, 300, 200, 100, 0))
S_4 = ggarrange(STABS.corr.plt, STABS.mse.plt, STABS.bias.plt, common.legend = TRUE, ncol = 1, legend = "left", labels = c("A", "B", "C"))
ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_3_4.png", S_4, width = 7, height = 8, dpi = 300)

####################################################################################

library(ggplot2)
library(dplyr)
library(ggpubr)

options(scipen = 999)
# aPGS correlation
STABS.corr = read.table("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/otuput/stabilizing_selection/metrics/HT_BMI.cor", header = TRUE)
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
  labs(title = "aPGS accuracy", x = expression("generations before the present, " ~ tau), y = expression(italic(r)^2(Y[i], hat(Y)[i]), list(j == "")), fill = "") +
  scale_fill_manual(values = c("#DDCC77", "#44AA99"), labels = c("Height", "BMI")) + ylim(0,1)

# aPS MSE
STABS.mse = read.table("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/otuput/stabilizing_selection/metrics/HT_BMI_scaled.MSE", header = TRUE)
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
  labs(title = "aPGS MSE", x = expression("generations before the present, " ~ tau), y = expression(italic(MSE)(Y[i], hat(Y)[i]), list(j == "")), fill = "") +
  scale_fill_manual(values = c("#DDCC77", "#44AA99"), labels = c("Height", "BMI")) + ylim(0,1.5)

# aPS bias
STABS.bias = read.table('/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/otuput/stabilizing_selection/metrics/HT_BMI.bias', header = TRUE)
STABS.bias$generation = factor(STABS.bias$generation, levels = c(400, 300, 200, 100, 0))
STABS.bias$pheno = factor(STABS.bias$hsq, levels = c("height", "bmi"))

bias.real = ggplot(STABS.bias, aes(x = generation, y = bias, fill = pheno)) +
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
  labs(title = "aPGS Bias", x = expression("generations before the present, " ~ tau), y = expression(italic(bias)(Y[i], hat(Y)[i]), list(j == ""))) +
  scale_fill_manual(values = c("#DDCC77", "#44AA99"), labels = c("Height", "BMI")) + ylim(-1.5,1.5)


Fig_6 = ggarrange(corr.real, mse.real, bias.real, common.legend = TRUE, ncol = 1, legend = "left", labels = c("A", "B", "C"))
ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/otuput/new_figures/Figure_6.png", Fig_6, width = 8, height = 8, dpi = 300)


################################################################################################
################################################################################################
STABS.corr = read.table("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/realistic_params/metrics/tPheno_aPS.cor", header = TRUE)
STABS.corr$generation = factor(STABS.corr$generation, levels = c(200000, 200100, 200200, 200300, 200400))
STABS.corr$pheno = factor(STABS.corr$hsq, levels = c("HT", "BMI"))
STABS.corr$r.squared = (STABS.corr$corr)^2
corr = ggplot(STABS.corr, aes(x = generation, y = r.squared, fill = pheno)) +
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
  labs(title = "aPGS accuracy", x = expression("generations before the present, " ~ tau), y = expression(italic(r)^2(Y[i], hat(Y)[i]), list(j == "")), fill = "") +
  scale_fill_manual(values = c("#DDCC77", "#44AA99"), labels = c("Height", "BMI")) + ylim(0,1) +
  scale_x_discrete(breaks = c(200000, 200100, 200200, 200300, 200400), labels = c(400, 300, 200, 100, 0))


STABS.mse = read.table("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/realistic_params/metrics/tPheno_aPS_scaled.mse", header = TRUE)
STABS.mse$generation = factor(STABS.mse$generation, levels = c(200000, 200100, 200200, 200300, 200400))
STABS.mse$pheno = factor(STABS.mse$hsq, levels = c("HT", "BMI"))
mse = ggplot(STABS.mse, aes(x = generation, y = round(MSE,2), fill = pheno)) +
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
  labs(title = "aPGS MSE", x = expression("generations before the present, " ~ tau), y = expression(italic(MSE)(Y[i], hat(Y)[i]), list(j == "")), fill = "") +
  scale_fill_manual(values = c("#DDCC77", "#44AA99"), labels = c("Height", "BMI")) + 
  scale_x_discrete(breaks = c(200000, 200100, 200200, 200300, 200400), labels = c(400, 300, 200, 100, 0))

STABS.bias = read.table("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/realistic_params/metrics/tPheno_aPS.bias", header = TRUE)
STABS.bias$generation = factor(STABS.bias$generation, levels = c(200000, 200100, 200200, 200300, 200400))
STABS.bias$pheno = factor(STABS.bias$hsq, levels = c("HT", "BMI"))
bias = ggplot(STABS.bias, aes(x = generation, y = bias, fill = pheno)) +
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
  labs(title = "aPGS Bias", x = expression("generations before the present, " ~ tau), y = expression(italic(Bias)(Y[i], hat(Y)[i]), list(j == ""))) +
  scale_fill_manual(values = c("#DDCC77", "#44AA99"), labels = c("Height", "BMI")) + ylim(-1.5,1.5) + 
  scale_x_discrete(breaks = c(200000, 200100, 200200, 200300, 200400), labels = c(400, 300, 200, 100, 0))

Fig_5 = ggarrange(corr, mse, bias, ncol = 1, align = "v", labels = c("A", "B", "C"), common.legend = TRUE, legend = "left")
ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_5.png", Fig_5, width = 8, height = 8, dpi = 300)

# derived allele freq vs effect sizes

process_QTLs = function(path, type, case, seed = 41, sample_size = 3) {
  set.seed(seed)
  chosen_reps = sample(1:100, 3)
  QTLs_data = data.frame()
  for (rep in chosen_reps) {
    file_path = paste0(path, type, "/", case, "/QTLs_200400_", rep, ".txt")
    tmp = read.table(file_path, header = TRUE)
    tmp$rep = rep
    QTLs_data = rbind(QTLs_data, tmp)
  }
  plot_list = list()
  index = 1
  for (rep in chosen_reps) {
    rep_data = QTLs_data[which(QTLs_data$rep == rep),]
    plt = ggplot(rep_data, aes(x = freq, y = selcoeff)) +
      geom_point() +
      theme_linedraw() + 
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "left",
            legend.title = element_text(hjust = 0.5, size = 12),
            legend.text = element_text(size = 12),
            strip.text.x = element_text(size = 14)) +
      labs(title = "", x = "derived allele frequency", y = "effect size")
    plot_list[[index]] <- plt
    index = index + 1
  }
  return(plot_list)
  print(chosen_reps)
}

path = "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/"
plot_list_w1 = process_QTLs(path, "stabilizing_selection", "mid/1")
plot_list_w5 = process_QTLs(path, "stabilizing_selection", "mid/5")

fplot = ggarrange(ggarrange(plot_list_w1[[1]], plot_list_w1[[2]], plot_list_w1[[3]], ncol = 1),
                  ggarrange(plot_list_w5[[1]], plot_list_w5[[2]], plot_list_w5[[3]], ncol = 1), 
                  ncol = 2)

ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_X_1.png", fplot, width = 8, height = 8, dpi = 300)

path = "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/"
plot_list_d025 = process_QTLs(path, "directional_selection", "mid/025")
plot_list_d00025 = process_QTLs(path, "directional_selection", "mid/00025")

fplot = ggarrange(ggarrange(plot_list_d025[[1]], plot_list_d025[[2]], plot_list_d025[[3]], ncol = 1),
                  ggarrange(plot_list_d00025[[1]], plot_list_d00025[[2]], plot_list_d00025[[3]], ncol = 1), 
                  ncol = 2)

ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_X_2.png", fplot, width = 8, height = 8, dpi = 300)



