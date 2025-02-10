# author: Valeria Añorve-Garibay
# Directional Selection: Figures

library(ggplot2)
library(dplyr)
library(ggpubr)
options(scipen = 999)

hsq.values = c("total", "mid")
QTLsSD.values = c("025", "0025", "00025")

# average phenotype mean and variance after optimum shift
data.plt = data.frame()
for (hsq in hsq.values) {
  for (QTLsSD in QTLsSD.values){
    data = data.frame()
    for (rep in 1:100) {
      OUT.n = readLines(paste("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/directional_selection/", hsq, "/", QTLsSD, "/OUT.", rep, sep = ""))
      OUT.n = OUT.n[!grepl("COALESCED", OUT.n)]
      OUT.n = OUT.n[35:length(OUT.n)]
      OUT.n = read.csv(text = OUT.n, header = FALSE, sep = " ")
      OUT.n = data.frame(gen = as.integer(gsub(":", "", OUT.n$V1)), mean_pheno = as.numeric(gsub(",", "", OUT.n$V4)), var_pheno = as.numeric(OUT.n$V7))
      OUT.n$rep = rep
      
      data = rbind(data, OUT.n)
    }
    data.summ = data %>%
      group_by(gen) %>%
      summarize(
        mean = mean(mean_pheno),
        var = mean(var_pheno)
      )
    data.summ$QTLsSD = QTLsSD
    data.summ$hsq = hsq
    data.plt = rbind(data.plt, data.summ)
  }
}
data.plt = data.plt %>%
  mutate(hsq_math = factor(hsq,
                           levels = c("total", "mid"),
                           labels = c(expression(h^2 == 1), expression(h^2 == 0.5))))
data.plt$QTLsSD = factor(data.plt$QTLsSD, levels = c("025", "0025", "00025"))
data.plt$pgen = abs(data.plt$gen - 200000)
mean.plt = ggplot(data.plt, aes(x = pgen, y = mean, color = QTLsSD)) +
  geom_point() + theme_linedraw() + 
  facet_wrap(~hsq_math, labeller = label_parsed) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  labs(title = "", x = "generations since optimum shift", y = expression(mean(bar(Y)[j])), color = expression("QTLs " * sigma * "")) +
  scale_color_manual(values = c("#E69F00", "#F0E442", "#009E73"), labels = c("0.25", "0.025", "0.0025")) + ylim(0,1) +
  scale_y_continuous(breaks = seq(0,1,0.2))

var.plt = ggplot(data.plt, aes(x = pgen, y = var, color = QTLsSD)) +
  geom_point() +
  theme_linedraw() + 
  facet_wrap(~hsq_math, labeller = label_parsed) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  labs(title = "", x = "generations since optimum shift", y = expression(var(Y[j])), color = expression("QTLs " * sigma * "")) +
  scale_color_manual(values = c("#E69F00", "#F0E442", "#009E73"), labels = c("0.25", "0.025", "0.0025")) + ylim(0,0.5)

Fig_4 = ggarrange(mean.plt, var.plt, nrow = 2, common.legend = TRUE, labels = c("A", "B"), legend = "left", align = "v")
ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_4.png", Fig_4, width = 8, height = 6.5, dpi = 300)

# metrics
# aPS correlation
DIR.corr = read.table('/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/directional_selection/metrics/tPheno_aPS.cor', header = TRUE, colClasses = c("integer", "integer", "numeric", "character", "character", "character"))
generations = c(200400, 200300, 200200, 200100, 200090, 200080, 200070, 200060, 200050, 200040, 200030, 200020, 200010, 200000)
DIR.corr$hsq = factor(DIR.corr$hsq, levels = c("total", "mid"))
DIR.corr$QTLsd = factor(DIR.corr$QTLsd, levels = c("025", "0025", "00025"))
DIR.corr$r.squared = (DIR.corr$corr)^2
DIR.corr = DIR.corr %>%
  mutate(hsq_math = factor(hsq,
                           levels = c("total", "mid"),
                           labels = c(expression(h^2 == 1), expression(h^2 == 0.5))))
# Figure 4_2
DIR.corr.100 = DIR.corr[which(DIR.corr$present_sample_size == 100),]
DIR.corr.100.wZoom.plt = ggplot(DIR.corr.100[which(DIR.corr.100$QTLsd == '025'),], aes(x = factor(generation), y = r.squared, fill = QTLsd)) +
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
  labs(title = "aPGS accuracy", x = expression("generations before the present, " ~ tau), y = expression(italic(r)^2(Y[j], hat(Y)[j]), list(j == "")), fill = expression("QTLs " * sigma * "")) +
  scale_fill_manual(values = c("#E69F00", "#F0E442", "#009E73"), labels = c("0.25", "0.025", "0.0025")) + ylim(0,1) +
  scale_x_discrete(breaks = generations, labels = abs(generations-200400))

# Figure 4
DIR.corr.100.noZoom = DIR.corr.100[DIR.corr.100$generation %in% c(200400, 200300, 200200, 200100, 200000), ]
DIR.corr.100.noZoom$generation = factor(DIR.corr.100.noZoom$generation, levels = c(200000, 200100, 200200, 200300, 200400))
DIR.corr.100.noZoom.plt = ggplot(DIR.corr.100.noZoom, aes(x = generation, y = r.squared, fill = QTLsd)) +
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
  labs(title = "aPGS accuracy", x = expression("generations before the present, " ~ tau), y = expression(italic(r)^2(Y[j], hat(Y)[j]), list(j == "")), fill = expression("QTLs " * sigma * "")) +
  scale_fill_manual(values = c("#E69F00", "#F0E442", "#009E73"), labels = c("0.25", "0.025", "0.0025")) + ylim(0,1) +
  scale_x_discrete(breaks = c(200000, 200100, 200200, 200300, 200400), labels = c(400, 300, 200, 100, 0))

# MSE
DIR.MSE = read.table('/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/directional_selection/metrics/tPheno_aPS_scaled.mse', header = TRUE, colClasses = c("integer", "integer", "numeric", "character", "character", "character"))
DIR.MSE$hsq = factor(DIR.MSE$hsq, levels = c("total", "mid"))
DIR.MSE$QTLsd = factor(DIR.MSE$QTLsd, levels = c("025", "0025", "00025"))
DIR.MSE = DIR.MSE %>%
  mutate(hsq_math = factor(hsq,
                           levels = c("total", "mid"),
                           labels = c(expression(h^2 == 1), expression(h^2 == 0.5))))
# Figure 4_2
DIR.MSE.100 = DIR.MSE[which(DIR.MSE$present_sample_size == 100),]
DIR.MSE.100.wZoom.plt = ggplot(DIR.MSE.100[which(DIR.MSE.100$QTLsd == '025'),], aes(x = factor(generation), y = round(MSE,2), fill = QTLsd)) +
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
  labs(title = "aPGS MSE", x = expression("generations before the present, " ~ tau), y = expression(italic(MSE)(Y[j], hat(Y)[j]), list(j == "")), fill = expression("QTLs " * sigma * "")) +
  scale_fill_manual(values = c("#E69F00", "#F0E442", "#009E73"), labels = c("0.25", "0.025", "0.0025")) +
  scale_x_discrete(breaks = generations, labels = abs(generations-200400))

# Figure 4
DIR.MSE.100.noZoom = DIR.MSE.100[DIR.MSE.100$generation %in% c(200000, 200100, 200200, 200300, 200400), ]
DIR.MSE.100.noZoom$generation = factor(DIR.MSE.100.noZoom$generation, levels = c(200000, 200100, 200200, 200300, 200400))
DIR.MSE.100.noZoom.plt = ggplot(DIR.MSE.100.noZoom, aes(x = generation, y = round(MSE,2), fill = QTLsd)) +
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
  labs(title = "aPGS MSE", x = expression("generations before the present, " ~ tau), y = expression(italic(MSE)(Y[j], hat(Y)[j]), list(j == "")), fill = expression("QTLs " * sigma * "")) +
  scale_fill_manual(values = c("#E69F00", "#F0E442", "#009E73"), labels = c("0.25", "0.025", "0.0025")) +
  scale_x_discrete(breaks = c(200000, 200100, 200200, 200300, 200400), labels = c(400, 300, 200, 100, 0))
# Bias
DIR.bias = read.table('/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/directional_selection/metrics/tPheno_aPS.bias', header = TRUE, colClasses = c("integer", "integer", "numeric", "character", "character", "character"))
DIR.bias$hsq = factor(DIR.bias$hsq, levels = c("total", "mid"))
DIR.bias$QTLsd = factor(DIR.bias$QTLsd, levels = c("025", "0025", "00025"))
DIR.bias = DIR.bias %>%
  mutate(hsq_math = factor(hsq,
                           levels = c("total", "mid"),
                           labels = c(expression(h^2 == 1), expression(h^2 == 0.5))))
# Figure 4_2
DIR.bias.100 = DIR.bias[which(DIR.bias$present_sample_size == 100),]
DIR.bias.100.wZoom.plt = ggplot(DIR.bias.100[which(DIR.bias.100$QTLsd == '025'),], aes(x = factor(generation), y = bias, fill = QTLsd)) +
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
  labs(title = "aPGS Bias", x = expression("generations before the present, " ~ tau), y = expression(italic(Bias)(Y[i], hat(Y)[i]), list(j == "")), fill = expression("QTLs " * sigma * "")) +
  scale_fill_manual(values = c("#E69F00", "#F0E442", "#009E73"), labels = c("0.25", "0.025", "0.0025")) +
  scale_x_discrete(breaks = generations, labels = abs(generations-200400))

DIR.bias.100.noZoom = DIR.bias.100[DIR.bias.100$generation %in% c(200000, 200100, 200200, 200300, 200400), ]
DIR.bias.100.noZoom$generation = factor(DIR.bias.100.noZoom$generation, levels = c(200000, 200100, 200200, 200300, 200400))
DIR.bias.100.noZoom.plt = ggplot(DIR.bias.100.noZoom, aes(x = generation, y = bias, fill = QTLsd)) +
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
  labs(title = "aPGS Bias", x = expression("generations before the present, " ~ tau), y = expression(italic(Bias)(Y[i], hat(Y)[i]), list(j == "")), fill = expression("QTLs " * sigma * "")) +
  scale_fill_manual(values = c("#E69F00", "#F0E442", "#009E73"), labels = c("0.25", "0.025", "0.0025")) +
  scale_x_discrete(breaks = c(200000, 200100, 200200, 200300, 200400), labels = c(400, 300, 200, 100, 0))

Fig_4 = ggarrange(DIR.corr.100.noZoom.plt, DIR.MSE.100.noZoom.plt, DIR.bias.100.noZoom.plt, common.legend = TRUE, ncol = 1, legend = "left", labels = c("A", "B", "C"))
#ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_4.png", Fig_4, width = 8, height = 8, dpi = 300)

Fig_4_2 = ggarrange(DIR.corr.100.wZoom.plt, DIR.MSE.100.wZoom.plt, DIR.bias.100.wZoom.plt, common.legend = TRUE, ncol = 1, legend = "left", labels = c("A", "B", "C"))
#ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_4_2.png", Fig_4_2, width = 12, height = 10, dpi = 300)


median(DIR.corr.100.noZoom[which(DIR.corr.100.noZoom$hsq == "total" & DIR.corr.100.noZoom$generation == 200400 & DIR.corr.100.noZoom$QTLsd == "025"),"r.squared"])
median(DIR.corr.100.noZoom[which(DIR.corr.100.noZoom$hsq == "total" & DIR.corr.100.noZoom$generation == 200300 & DIR.corr.100.noZoom$QTLsd == "025"),"r.squared"])
median(DIR.corr.100.noZoom[which(DIR.corr.100.noZoom$hsq == "total" & DIR.corr.100.noZoom$generation == 200200 & DIR.corr.100.noZoom$QTLsd == "025"),"r.squared"])
median(DIR.corr.100.noZoom[which(DIR.corr.100.noZoom$hsq == "total" & DIR.corr.100.noZoom$generation == 200100 & DIR.corr.100.noZoom$QTLsd == "025"),"r.squared"])
median(DIR.corr.100.noZoom[which(DIR.corr.100.noZoom$hsq == "total" & DIR.corr.100.noZoom$generation == 200000 & DIR.corr.100.noZoom$QTLsd == "025"),"r.squared"])

# metrics w bigger sample size
DIR.corr.025 = DIR.corr[which(DIR.corr$hsq == "mid" & DIR.corr$QTLsd == "025"),]
DIR.corr.025$present_sample_size = factor(DIR.corr.025$present_sample_size, levels = c(100, 500, 1000))
DIR.corr.025.plt = ggplot(DIR.corr.025, aes(x = factor(generation), y = r.squared, fill = present_sample_size)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8, position = position_dodge(0.8)) +
  theme_linedraw() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12)) +
  labs(title = "aPGS accuracy", x = expression("generations before the present, " ~ tau), y = expression(italic(r)^2(Y[i], hat(Y)[i]), list(j == "")), fill = expression("sample size at " ~ tau * " = 0")) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442", "#E69F00")) + ylim(0,1) +
  scale_x_discrete(breaks = generations, labels = abs(generations-200400))

DIR.MSE.025 = DIR.MSE[which(DIR.MSE$hsq == "mid" & DIR.MSE$QTLsd == "025"),]
DIR.MSE.025$present_sample_size = factor(DIR.MSE.025$present_sample_size, levels = c(100, 500, 1000))
DIR.MSE.025.plt = ggplot(DIR.MSE.025, aes(x = factor(generation), y = MSE, fill = present_sample_size)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8, position = position_dodge(0.8)) +
  theme_linedraw() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12)) + 
  labs(title = "aPGS MSE", x = expression("generations before the present, " ~ tau), y = expression(italic(MSE)(Y[i], hat(Y)[i]), list(j == "")), fill = expression("sample size at " ~ tau * " = 0")) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442", "#E69F00")) +
  scale_x_discrete(breaks = generations, labels = abs(generations-200400))


DIR.bias.025 = DIR.bias[which(DIR.bias$hsq == "mid" & DIR.bias$QTLsd == "025"),]
DIR.bias.025$present_sample_size = factor(DIR.bias.025$present_sample_size, levels = c(100, 500, 1000))
DIR.bias.025.plt = ggplot(DIR.bias.025, aes(x = factor(generation), y = bias, fill = present_sample_size)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8, position = position_dodge(0.8)) +
  theme_linedraw() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12)) + 
  labs(title = "aPGS Bias", x = expression("generations before the present, " ~ tau), y = expression(italic(Bias)(Y[i], hat(Y)[i]), list(j == "")), fill = expression("sample size at " ~ tau * " = 0")) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442", "#E69F00")) +
  scale_x_discrete(breaks = generations, labels = abs(generations-200400))


Fig_4_11 = ggarrange(DIR.corr.025.plt, DIR.MSE.025.plt, DIR.bias.025.plt, common.legend = TRUE, ncol = 1, legend = "left", labels = c("A", "B", "C"))
ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_4_11.png", Fig_4_11, width = 8, height = 8, dpi = 300)



############################################################################
############################################################################
############################################################################
############################################################################

# QTL: effect sizes
plot_QTLs = function(data, QTLsd, color, generations_sorted, hsq, zoom = FALSE){
  data = data[which(data$QTLsd == QTLsd),]
  data = data[which(data$hsq == hsq),]
  
  # subset data if not zoomed
  if (!zoom) {
    data = data[data$gbp %in% c(400, 300, 200, 100, 0), ]
  }
  # boxplot for QTLs: C vs L
  CvsL = ggplot(data, aes(x = factor(bin, levels = unique(bin)), y = nQTLs, fill = group)) +
    geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
    theme_linedraw() + 
    facet_grid(cols = vars(factor(gbp, levels = generations_sorted)), rows = vars(QTLsSD_math), scales="free_y", switch = 'y',
               labeller = label_parsed) +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 14),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "effect size (s)", y = "number of QTL mutations") +
    scale_fill_manual(values = c("#D55E00", "#56B4E9"), labels = c("Conserved", "Lost")) +
    scale_x_discrete(expand = c(0.1, 0.2))
  # boxplot for QTLs ratio
  data = data[which(data$group == "conserved"),]
  data[which(is.na(data$coeff) & data$gbp != 0),"coeff"] = 0
  Ratio = ggplot(data, aes(x = factor(bin, levels = unique(bin)), y = coeff)) +
    geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8, fill = color) +
    theme_linedraw() +
    facet_grid(rows = vars(QTLsSD_math), cols = vars(factor(gbp, levels = generations_sorted)), scales="free_y", switch = 'y',
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
    scale_x_discrete(expand = c(0.1, 0.2))
  Fig = ggarrange(CvsL, Ratio, ncol = 1, align = "v", labels = c("A", "B"))
}
generations_sorted = c(seq(400, 300, -10), 200, 100, 0)
DIR.QTL.effect.sizes = data.frame()
for (genbp in generations_sorted) {
  tmp = read.table(paste("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/directional_selection/metrics/QTLs_summ_0_", genbp, ".txt", sep = ""), header = TRUE, colClasses = c("integer", "character", "integer", "character", "character", "integer", "character", "numeric"))
  DIR.QTL.effect.sizes = rbind(DIR.QTL.effect.sizes, tmp)
}
DIR.QTL.effect.sizes = DIR.QTL.effect.sizes %>%
  mutate(QTLsSD_math = factor(QTLsd,
                              levels = c("00025", "0025", "025"),
                              labels = c(expression("QTLs " * sigma * "" == 0.0025), expression("QTLs " * sigma * "" == 0.025), expression("QTLs " * sigma * "" == 0.25))))
# 0.25
# no Zoom
DIR.QTL.effect.sizes.025.1 = plot_QTLs(DIR.QTL.effect.sizes, "025", "#E69F00", generations_sorted, "total", zoom = FALSE)
DIR.QTL.effect.sizes.025.05 = plot_QTLs(DIR.QTL.effect.sizes, "025", "#E69F00", generations_sorted, "mid", zoom = FALSE)
# w Zoom
DIR.QTL.effect.sizes.025.1.Zoom = plot_QTLs(DIR.QTL.effect.sizes, "025", "#E69F00", generations_sorted, "total", zoom = TRUE)
DIR.QTL.effect.sizes.025.05.Zoom = plot_QTLs(DIR.QTL.effect.sizes, "025", "#E69F00", generations_sorted, "mid", zoom = TRUE)
# 0.025
# no Zoom
DIR.QTL.effect.sizes.0025.1 = plot_QTLs(DIR.QTL.effect.sizes, "0025", "#F0E442", generations_sorted, "total", zoom = FALSE)
DIR.QTL.effect.sizes.0025.05 = plot_QTLs(DIR.QTL.effect.sizes, "0025", "#F0E442", generations_sorted, "mid", zoom = FALSE)
# 0.0025
# no Zoom
DIR.QTL.effect.sizes.00025.1 = plot_QTLs(DIR.QTL.effect.sizes, "00025", "#009E73", generations_sorted, "total", zoom = FALSE)
DIR.QTL.effect.sizes.00025.05 = plot_QTLs(DIR.QTL.effect.sizes, "00025", "#009E73", generations_sorted, "mid", zoom = FALSE)

ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_4_3.png", DIR.QTL.effect.sizes.025.1, width = 15, height = 12, dpi = 300)
ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_4_4.png", DIR.QTL.effect.sizes.025.1.Zoom, width = 35, height = 18, dpi = 600, limitsize = FALSE)
ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_4_5.png", DIR.QTL.effect.sizes.025.05, width = 15, height = 12, dpi = 300)
ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_4_6.png", DIR.QTL.effect.sizes.025.05.Zoom, width = 35, height = 18, dpi = 600, limitsize = FALSE)

ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_4_7.png", DIR.QTL.effect.sizes.0025.1, width = 15, height = 12, dpi = 300)
ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_4_8.png", DIR.QTL.effect.sizes.0025.05, width = 15, height = 12, dpi = 300)

ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_4_9.png", DIR.QTL.effect.sizes.00025.1, width = 15, height = 12, dpi = 300)
ggsave("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/figures/Figure_4_10.png", DIR.QTL.effect.sizes.00025.05, width = 15, height = 12, dpi = 300)
