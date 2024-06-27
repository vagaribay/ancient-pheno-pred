# author: Valeria Añorve-Garibay
# Directional Selection: Figures

library(ggplot2)
library(dplyr)
library(ggpubr)
options(scipen = 999)

hsq.values = c("total", "mid")
QTLsSD.values = c("025_zoom", "0025", "00025")

# Average phenotype mean and variance after optimum shift
data.plt = data.frame()
for (hsq in hsq.values) {
  for (QTLsSD in QTLsSD.values){
    data = data.frame()
    for (rep in 1:100) {
      OUT.n = readLines(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/directional_selection/", hsq, "/", QTLsSD, "/OUT.", rep, sep = ""))
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
data.plt$QTLsSD = factor(data.plt$QTLsSD, levels = c("025_zoom", "0025", "00025"))
mean.plt = ggplot(data.plt, aes(x = gen-100000, y = mean, color = QTLsSD)) +
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
  
var.plt = ggplot(data.plt, aes(x = gen-100000, y = var, color = QTLsSD)) +
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
ggsave("/Users/valeriagby/desktop/ancient-pheno-prediction/figures/Figure_4.png", Fig_4, width = 8, height = 6.5, dpi = 300)

# aPS correlation
fname = "tPheno_aPS"
DIR.corr.tmp = read.table("/Users/valeriagby/desktop/ancient-pheno-prediction/output/directional_selection/metrics/tPheno_aPS_wzoom.cor", header = TRUE, colClasses = c("integer", "integer", "numeric", "character", "character"))
DIR.corr = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/directional_selection/metrics/", fname, ".cor", sep = ""), header = TRUE, colClasses = c("integer", "integer", "numeric", "character", "character"))
DIR.corr = rbind(DIR.corr.tmp, DIR.corr)
DIR.corr = DIR.corr[-which(DIR.corr$QTLsSD == "025"),]

# DIR without zoom
zoomxgen = seq(310, 390, 10)
DIR.corr_noZoom = DIR.corr %>% filter(!generation %in% zoomxgen)
DIR.corr_noZoom$generation = factor(DIR.corr_noZoom$generation, levels = c(400, 300, 200, 100, 0))
DIR.corr_noZoom$hsq = factor(DIR.corr_noZoom$hsq, levels = c("total", "mid"))
DIR.corr_noZoom$QTLsSD = factor(DIR.corr_noZoom$QTLsSD, levels = c("025_zoom", "0025", "00025"))
DIR.corr_noZoom$r.squared = (DIR.corr_noZoom$corr)^2
DIR.corr_noZoom = DIR.corr_noZoom %>%
  mutate(hsq_math = factor(hsq,
                           levels = c("total", "mid"),
                           labels = c(expression(h^2 == 1), expression(h^2 == 0.5))))

DIR.corr.plt = ggplot(DIR.corr_noZoom, aes(x = generation, y = r.squared, fill = QTLsSD)) +
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
  labs(title = "aPS accuracy", x = expression("generations before the present, " ~ tau), y = expression(italic(r)^2(Y[j], hat(Y)[j]), list(j == "")), fill = expression("QTLs " * sigma * "")) +
  scale_fill_manual(values = c("#E69F00", "#F0E442", "#009E73"), labels = c("0.25", "0.025", "0.0025")) + ylim(0,1)
# DIR with zoom at 0.25
#zoomxgen = seq(300, 400, 10)
#DIR.corr_wZoom = DIR.corr %>% filter(generation %in% zoomxgen)
DIR.corr_wZoom = DIR.corr[which(DIR.corr$QTLsSD == "025_zoom"),]
DIR.corr_wZoom$generation = factor(DIR.corr_wZoom$generation, levels = c(seq(400, 300, -10), 200, 100, 0))
DIR.corr_wZoom$hsq = factor(DIR.corr_wZoom$hsq, levels = c("total", "mid"))
DIR.corr_wZoom$QTLsSD = factor(DIR.corr_wZoom$QTLsSD, levels = c("025_zoom"))
DIR.corr_wZoom$r.squared = (DIR.corr_wZoom$corr)^2
DIR.corr_wZoom = DIR.corr_wZoom %>%
  mutate(hsq_math = factor(hsq,
                           levels = c("total", "mid"),
                           labels = c(expression(h^2 == 1), expression(h^2 == 0.5))))

DIR.corr.pltWzoom = ggplot(DIR.corr_wZoom, aes(x = generation, y = r.squared, fill = QTLsSD)) +
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
  labs(title = "aPS accuracy", x = expression("generations before the present, " ~ tau), y = expression(italic(r)^2(Y[j], hat(Y)[j]), list(j == "")), fill = expression("QTLs " * sigma * "")) +
  scale_fill_manual(values = c("#E69F00", "#F0E442", "#009E73"), labels = c("0.25", "0.025", "0.0025")) + ylim(0,1)

# aPS MSE
fname = "tPheno_aPS"
DIR.MSE.tmp = read.table("/Users/valeriagby/desktop/ancient-pheno-prediction/output/directional_selection/metrics/tPheno_aPS_wzoom.mse", header = TRUE, colClasses = c("integer", "integer", "numeric", "character", "character"))
DIR.mse = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/directional_selection/metrics/", fname, ".mse", sep = ""), header = TRUE, colClasses = c("integer", "integer", "numeric", "character", "character"))
DIR.mse = rbind(DIR.MSE.tmp, DIR.mse)
DIR.mse = DIR.mse[-which(DIR.mse$QTLsSD == "025"),]

# DIR without zoom
zoomxgen = seq(310, 390, 10)
DIR.mse_noZoom = DIR.mse %>% filter(!generation %in% zoomxgen)
DIR.mse_noZoom$generation = factor(DIR.mse_noZoom$generation, levels = c(400, 300, 200, 100, 0))
DIR.mse_noZoom$hsq = factor(DIR.mse_noZoom$hsq, levels = c("total", "mid"))
DIR.mse_noZoom$QTLsSD = factor(DIR.mse_noZoom$QTLsSD, levels = c("025_zoom", "0025", "00025"))
DIR.mse_noZoom = DIR.mse_noZoom %>%
  mutate(hsq_math = factor(hsq,
                           levels = c("total", "mid"),
                           labels = c(expression(h^2 == 1), expression(h^2 == 0.5))))
DIR.mse.plt = ggplot(DIR.mse_noZoom, aes(x = generation, y = round(MSE,2), fill = QTLsSD)) +
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
  labs(title = "aPS MSE", x = expression("generations before the present, " ~ tau), y = expression(italic(MSE)(Y[j], hat(Y)[j]), list(j == "")), fill = expression("QTLs " * sigma * "")) +
  scale_fill_manual(values = c("#E69F00", "#F0E442", "#009E73"), labels = c("0.25", "0.025", "0.0025")) + ylim(0, 0.8)

# DIR with zoom at 0.25
#zoomxgen = seq(300, 400, 10)
#DIR.mse_wZoom = DIR.mse %>% filter(generation %in% zoomxgen)
DIR.mse_wZoom = DIR.mse[which(DIR.mse$QTLsSD == "025_zoom"),]
DIR.mse_wZoom$generation = factor(DIR.mse_wZoom$generation, levels = c(seq(400, 300, -10), 200, 100, 0))
DIR.mse_wZoom$hsq = factor(DIR.mse_wZoom$hsq, levels = c("total", "mid"))
DIR.mse_wZoom$QTLsSD = factor(DIR.mse_wZoom$QTLsSD, levels = c("025_zoom"))
DIR.mse_wZoom = DIR.mse_wZoom %>%
  mutate(hsq_math = factor(hsq,
                           levels = c("total", "mid"),
                           labels = c(expression(h^2 == 1), expression(h^2 == 0.5))))
DIR.mse.pltWzoom = ggplot(DIR.mse_wZoom, aes(x = generation, y = round(MSE,2), fill = QTLsSD)) +
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
  labs(title = "aPS MSE", x = expression("generations before the present, " ~ tau), y = expression(italic(MSE)(Y[j], hat(Y)[j]), list(j == "")), fill = expression("QTLs " * sigma * "")) +
  scale_fill_manual(values = c("#E69F00", "#F0E442", "#009E73"), labels = c("0.25", "0.025", "0.0025"))

Fig_5 = ggarrange(DIR.corr.plt, DIR.mse.plt, common.legend = TRUE, ncol = 1, legend = "left", labels = c("A", "B"))
Fig5_1 = ggarrange(DIR.corr.pltWzoom, DIR.mse.pltWzoom, common.legend = TRUE, ncol = 1, legend = "left", labels = c("A", "B"))
ggsave("/Users/valeriagby/desktop/ancient-pheno-prediction/figures/Figure_5.png", Fig_5, width = 8, height = 6.5, dpi = 300)
ggsave("/Users/valeriagby/desktop/ancient-pheno-prediction/figures/Figure_5_1.png", Fig5_1, width = 15, height = 8, dpi = 300)

# QTL: effect sizes
DIR.QTL.025 = data.frame()
for (genbp in c(seq(400, 300, -10), 200, 100, 0)) {
  tmp.025 = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/directional_selection/metrics/QTLs_summ_025wzoom_0_", genbp, ".txt", sep = ""), header = TRUE, colClasses = c("integer", "character", "integer", "character", "character", "integer", "character", "numeric"))
  DIR.QTL.025 = rbind(DIR.QTL.025, tmp.025)
}
# 0.25
DIR.QTL.025 = DIR.QTL.025[which(DIR.QTL.025$hsq == "mid"),]
DIR.QTL.025$bin = factor(DIR.QTL.025$bin, levels = unique(DIR.QTL.025$bin))
DIR.QTL.025 = DIR.QTL.025 %>%
  mutate(QTLsSD_math = factor(QTLsSD,
                              levels = c("025_zoom"),
                              labels = c(expression("QTLs " * sigma * "" == 0.25))))

#DIR.QTL.025 = DIR.QTL.025 %>% filter(!gbp %in% seq(310, 390, 10))
DIR.QTL.025 = DIR.QTL.025 %>% filter(gbp %in% seq(300, 400, 10))
DIR.QTL.025$gbp = factor(DIR.QTL.025$gbp, levels = seq(400, 300, -10))
#DIR.QTL.025$gbp = factor(DIR.QTL.025$gbp, levels = c(400, 300, 200, 100, 0))
DIR.QTL.025.plt = ggplot(DIR.QTL.025, aes(x = bin, y = nQTLs, fill = group)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() + 
  facet_grid(cols = vars(gbp), rows = vars(QTLsSD_math), scales="free_y", switch = 'y',
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
  scale_fill_manual(values = c("#D55E00", "#56B4E9"), labels = c("Conserved", "Lost"))

DIR.QTL.025.C = DIR.QTL.025[which(DIR.QTL.025$group == "conserved"),]
RATIO = ggplot(DIR.QTL.025.C, aes(x = bin, y = coeff)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8, fill = "#E69F00") +
  theme_linedraw() +
  facet_grid(rows = vars(QTLsSD_math), cols = vars(gbp), scales="free_y", switch = 'y',
             labeller = label_parsed) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "effect size (s)", y = "Conserved QTLs / Conserved QTLs + Lost QTLs")

#supp = ggarrange(DIR.QTL.025.plt, RATIO, ncol = 1, align = "v", labels = c("A", "B"))
#ggsave("/Users/valeriagby/desktop/ancient-pheno-prediction/output/figures/Figure5_Supp_X.png", supp, dpi = 300, width = 25, height = 15)

######################################################################
######################################################################
DIR.QTL.0025 = data.frame()
DIR.QTL.00025 = data.frame()
for (genbp in c(400,300,200,100,0)) {
  tmp.0025 = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/directional_selection/metrics/QTLs_summ_0025_0_", genbp, ".txt", sep = ""), header = TRUE, colClasses = c("integer", "character", "integer", "character", "character", "integer", "character", "numeric"))
  tmp.0025$gbp = genbp
  tmp.00025 = read.table(paste("/Users/valeriagby/desktop/ancient-pheno-prediction/output/directional_selection/metrics/QTLs_summ_00025_0_", genbp, ".txt", sep = ""), header = TRUE, colClasses = c("integer", "character", "integer", "character", "character", "integer", "character", "numeric"))
  tmp.00025$gbp = genbp
  DIR.QTL.0025 = rbind(DIR.QTL.0025, tmp.0025)
  DIR.QTL.00025 = rbind(DIR.QTL.00025, tmp.00025)
}
DIR.QTL.0025 = DIR.QTL.0025[which(DIR.QTL.0025$hsq == "mid"),]
DIR.QTL.0025$bin = factor(DIR.QTL.0025$bin, levels = unique(DIR.QTL.0025$bin))
DIR.QTL.0025 = DIR.QTL.0025 %>%
  mutate(QTLsSD_math = factor(QTLsSD,
                              levels = c("0025"),
                              labels = c(expression("QTLs " * sigma * "" == 0.025))))
DIR.QTL.0025$gbp = factor(DIR.QTL.0025$gbp, c(400,300,200,100,0))
DIR.QTL.0025.plt = ggplot(DIR.QTL.0025, aes(x = bin, y = nQTLs, fill = group)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() + 
  facet_grid(cols = vars(gbp), rows = vars(QTLsSD_math), scales="free_y", switch = 'y',
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
  scale_fill_manual(values = c("#D55E00", "#56B4E9"), labels = c("Conserved", "Lost")) + ylim(0,210)

DIR.QTL.0025.C = DIR.QTL.0025[which(DIR.QTL.0025$group == "conserved"),]
RATIO = ggplot(DIR.QTL.0025.C, aes(x = bin, y = coeff)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8, fill = "#F0E442") +
  theme_linedraw() +
  facet_grid(rows = vars(QTLsSD_math), cols = vars(gbp), scales="free_y", switch = 'y',
             labeller = label_parsed) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "effect size (s)", y = "Conserved QTLs / Conserved QTLs + Lost QTLs")

#supp = ggarrange(DIR.QTL.0025.plt, RATIO, ncol = 1, align = "v", labels = c("A", "B"))
#ggsave("/Users/valeriagby/desktop/ancient-pheno-prediction/output/figures/Figure5_Supp_7.png", supp, width = 15, height = 12, dpi = 300)

# 0.0025
DIR.QTL.00025 = DIR.QTL.00025[which(DIR.QTL.00025$hsq == "mid"),]
DIR.QTL.00025$bin = factor(DIR.QTL.00025$bin, levels = unique(DIR.QTL.00025$bin))
DIR.QTL.00025 = DIR.QTL.00025 %>%
  mutate(QTLsSD_math = factor(QTLsSD,
                              levels = c("00025"),
                              labels = c(expression("QTLs " * sigma * "" == 0.0025))))
DIR.QTL.00025$gbp = factor(DIR.QTL.00025$gbp, c(400,300,200,100,0))
DIR.QTL.00025.plt = ggplot(DIR.QTL.00025, aes(x = bin, y = nQTLs, fill = group)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8) +
  theme_linedraw() + 
  facet_grid(cols = vars(gbp), rows = vars(QTLsSD_math), scales="free_y", switch = 'y',
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
  scale_fill_manual(values = c("#D55E00", "#56B4E9"), labels = c("Conserved", "Lost")) + ylim(0,210)


DIR.QTL.00025.C = DIR.QTL.00025[which(DIR.QTL.00025$group == "conserved"),]
RATIO = ggplot(DIR.QTL.00025.C, aes(x = bin, y = coeff)) +
  geom_boxplot(outlier.shape = 3, outlier.color = "red", outlier.size = 0.8, fill = "#009E73") +
  theme_linedraw() +
  facet_grid(rows = vars(QTLsSD_math), cols = vars(gbp), scales="free_y", switch = 'y',
             labeller = label_parsed) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "effect size (s)", y = "Conserved QTLs / Conserved QTLs + Lost QTLs")

#supp = ggarrange(DIR.QTL.00025.plt, RATIO, ncol = 1, align = "v", labels = c("A", "B"))
#ggsave("/Users/valeriagby/desktop/ancient-pheno-prediction/output/figures/Figure5_Supp_9.png", supp, width = 15, height = 12, dpi = 300)
