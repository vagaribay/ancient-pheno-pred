# author: Valeria Añorve-Garibay

library(vcfR)
library(stringr)
library(hash)
library(tidyr)

options(scipen = 999)

source("/Users/valeriagby/desktop/ancient-pheno-prediction/scripts/ancient-pheno-pred-FUNCTIONS.R")
path.new = "/Users/valeriagby/desktop/ancient-pheno-prediction/output/"

# Neutral Evolution
hsq.values = c("total", "mid")
path.tmp = paste(path.new, "neutral_evolution/", sep = "")
# compute ancient-polygenic scores
for (hsq in hsq.values) {
  print(paste("Running hsq = ", hsq, sep = ""))
  path = paste(path.tmp, hsq, "/", sep = "")
  nrep = 100
  compute_aPS(nrep, path)
}

# compute precision metrics: r and MSE
corr = data.frame()
MSE = data.frame()
for (hsq in hsq.values) {
  print(paste("Running hsq = ", hsq, sep = ""))
  path = paste(path.tmp, hsq, "/", sep = "")
  nrep = 100
  metrics = compute_metrics(nrep, path, hsq)
  corr = rbind(corr, metrics$pcorr)
  MSE = rbind(MSE, metrics$MSE)
}
write.table(corr, paste(path.tmp, "metrics/tPheno_aPS.cor", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(MSE, paste(path.tmp, "metrics/tPheno_aPS.mse", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE)

# QTLs: effect sizes per bin
QTLs.s = data.frame()
for (hsq in hsq.values) {
  print(paste("Running hsq = ", hsq, sep = ""))
  path = paste(path.tmp, hsq, "/", sep = "")
  nrep = 100
  tmp.out = compute_QTLsDist(nrep, path, gbp = 400, bins = seq(-0.8, 0.8, 0.1), hsq)
  QTLs.s = rbind(QTLs.s, tmp.out)
}
QTLs.s$coeff = (QTLs.s$nQTLs / QTLs.s$nQTLs.x.bin)
write.table(QTLs.s, paste(path.tmp, "metrics/tPheno_aPS_0.8-0.8-0.1.QTLs", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE)

# QTLs: frequencies per bin
QTLs.frq = data.frame()
for (hsq in hsq.values) {
  print(paste("Running hsq = ", hsq, sep = ""))
  path = paste(path.tmp, hsq, "/", sep = "")
  nrep = 100
  tmp.out = compute_QTLsFreq(nrep, path, gbp = 400, bins = seq(0, 1, 0.001), hsq)
  QTLs.frq = rbind(QTLs.frq, tmp.out)
}
write.table(QTLs.frq, paste(path.tmp, "metrics/tPheno_aPS_0-1-0.001.QTLs.frq", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE)

########################################################################################################################################################################################################
########################################################################################################################################################################################################

# Stabilizing Selection
hsq.values = c("total", "mid")
w.values = c("1", "2", "3", "4", "5")
path.tmp = paste(path.new, "stabilizing_selection/", sep = "")
# compute ancient-polygenic scores
for (hsq in hsq.values) {
  print(paste("Running hsq = ", hsq, sep = ""))
  for (w in w.values) {
    path = paste(path.tmp, hsq, "/", w, "/", sep = "")
    nrep = 100
    compute_aPS(nrep, path)
  }
}

# compute precision metrics: r and MSE
corr = data.frame()
MSE = data.frame()
for (hsq in hsq.values) {
  print(paste("Running hsq = ", hsq, sep = ""))
  corr.tmp = data.frame()
  MSE.tmp = data.frame()
  for (w in w.values) {
    path = paste(path.tmp, hsq, "/", w, "/", sep = "")
    nrep = 100
    metrics = compute_metrics(nrep, path, hsq)
    pcorr.tmp = metrics$pcorr
    pcorr.tmp$w = w

    pMSE.tmp = metrics$MSE
    pMSE.tmp$w = w
    
    corr.tmp = rbind(corr.tmp, pcorr.tmp)
    MSE.tmp = rbind(MSE.tmp, pMSE.tmp)
  }
  corr = rbind(corr, corr.tmp)
  MSE = rbind(MSE, MSE.tmp)
}
write.table(corr, paste(path.tmp, "metrics/tPheno_aPS.cor", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(MSE, paste(path.tmp, "metrics/tPheno_aPS.mse", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE)

# QTLs: effect sizes per bin
QTLs.s = data.frame()
for (hsq in hsq.values) {
  print(paste("Running hsq = ", hsq, sep = ""))
  QTLs.s.tmp = data.frame()
  for (w in w.values) {
   path = paste(path.tmp, hsq, "/", w, "/", sep = "")
    nrep = 100
    tmp.out = compute_QTLsDist(nrep, path, gbp = 400, bins = seq(-0.8, 0.8, 0.1), hsq)
    tmp.out$w = w
    QTLs.s.tmp = rbind(QTLs.s.tmp, tmp.out)
  }
  QTLs.s = rbind(QTLs.s, QTLs.s.tmp)
}
QTLs.s$coeff = (QTLs.s$nQTLs / QTLs.s$nQTLs.x.bin)
write.table(QTLs.s, paste(path.tmp, "metrics/QTLs_summ.txt", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE)

########################################################################################################################################################################################################
########################################################################################################################################################################################################

# Directional Selection
hsq.values = c("total", "mid")
#QTLsSD.values = c("025_zoom") # changes must be done in baseline functions
QTLsSD.values = c("0025", "00025")
path.tmp = paste(path.new, "directional_selection/", sep = "")
# compute ancient-polygenic scores
for (hsq in hsq.values) {
  print(paste("Running hsq = ", hsq, sep = ""))
  for (QTLsSD in QTLsSD.values) {
    path = paste(path.tmp, hsq, "/", QTLsSD, "/", sep = "")
    nrep = 100
    compute_aPS(nrep, path)
  }
}

# compute precision metrics: r and MSE
corr = data.frame()
MSE = data.frame()
for (hsq in hsq.values) {
  print(paste("Running hsq = ", hsq, sep = ""))
  corr.tmp = data.frame()
  MSE.tmp = data.frame()
  for (QTLsSD in QTLsSD.values) {
    path = paste(path.tmp, hsq, "/", QTLsSD, "/", sep = "")
    nrep = 100
    metrics = compute_metrics(nrep, path, hsq)
    pcorr.tmp = metrics$pcorr
    pcorr.tmp$QTLsSD = QTLsSD
    
    pMSE.tmp = metrics$MSE
    pMSE.tmp$QTLsSD = QTLsSD

    corr.tmp = rbind(corr.tmp, pcorr.tmp)
    MSE.tmp = rbind(MSE.tmp, pMSE.tmp)
  }
  corr = rbind(corr, corr.tmp)
  MSE = rbind(MSE, MSE.tmp)
}
write.table(corr, paste(path.tmp, "metrics/tPheno_aPS_wzoom.cor", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(MSE, paste(path.tmp, "metrics/tPheno_aPS_wzoom.mse", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE)

# QTLs: effect sizes per bin
#gbp.values = abs(c(100400, 100300, 100200, 100100, 100090, 100080, 100070, 100060, 100050, 100040, 100030, 100020, 100010, 100000) - 100400)
gbp.values = abs(c(100400, 100300, 100200, 100100, 100000) - 100400)

for (gbp_index in gbp.values) {
  QTLs.s = data.frame()
  for (hsq in hsq.values) {
    print(paste("Running hsq = ", hsq, sep = ""))
    QTLs.s.tmp = data.frame()
    for (QTLsSD in QTLsSD.values) {
      path = paste(path.tmp, hsq, "/", QTLsSD, "/", sep = "")
      nrep = 100
      tmp.out = compute_QTLsDist(nrep, path, gbp = gbp_index, bins = seq(-0.5, 0.5, 0.1), hsq)
      tmp.out$QTLsSD = QTLsSD
      QTLs.s.tmp = rbind(QTLs.s.tmp, tmp.out)
    }
    QTLs.s = rbind(QTLs.s, QTLs.s.tmp)
  }
  QTLs.s$coeff = (QTLs.s$nQTLs / QTLs.s$nQTLs.x.bin)
  QTLs.s$gbp = gbp_index
  write.table(QTLs.s, paste(path.tmp, "metrics/QTLs_summ_025wzoom_short_0_", gbp_index, ".txt", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE)
}

########################################################################################################################################################################################################
########################################################################################################################################################################################################

# Stabilizing Selection: realistic parameters for w and hsq
phenos = c("height", "bmi")
path.tmp = paste(path.new, "/real_estimates/", sep = "")
# compute ancient-polygenic scores
for (pheno in phenos) {
  path = paste(path.tmp, pheno, "/", sep = "")
  nrep = 100
  compute_aPS(nrep, path)
}
# compute precision metrics: r and MSE
corr = data.frame()
MSE = data.frame()
for (pheno in phenos) {
  print(paste("Running ", pheno, sep = ""))
  path = paste(path.tmp, pheno, "/", sep = "")
  nrep = 100
  metrics = compute_metrics(nrep, path, pheno)
  corr = rbind(corr, metrics$pcorr)
  MSE = rbind(MSE, metrics$MSE)
}
write.table(corr, paste(path.tmp, "metrics/tPheno_aPS.cor", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(MSE, paste(path.tmp, "metrics/tPheno_aPS.mse", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE)
