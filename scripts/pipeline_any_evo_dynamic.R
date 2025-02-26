# author: Valeria Añorve-Garibay
# the following commands were used to run functions to compute polygenic scores, precision metrics (mse, r^2, bias) and the number of lost/conserved QTLs per evolutionary scenario
# for details in how metrics are computed see ancient-pheno-pred-FUNCTIONS.R

library(vcfR)
library(stringr)
library(hash)
library(tidyr)

options(scipen = 999)

# load functions
source("/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/scripts/ancient-pheno-pred-FUNCTIONS_new.R")

source.path = "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/"

################################################################################
################################################################################

# Neutral Evolution
path.tmp = paste(source.path, "neutral_evolution/", sep = "")
hsq.values = c("total", "mid")
# aPS
#for (hsq in hsq.values) {
#  message("Running aPS for hsq = ", hsq)
#  path = paste0(path.tmp, hsq, "/")
#  nrep = 100
#  # Run compute_aPS with base parameters
#  compute_aPS(nrep, path, 0, 100)
#  if (hsq == "mid") {
#    compute_aPS(nrep, path, 0, 500)
#    compute_aPS(nrep, path, 0, 1000)
#  }
#}

# compute precision metrics: r^2, MSE and bias
corr = data.frame()
MSE = data.frame()
bias = data.frame()
sample.sizes = c(100, 500, 1000)
for (hsq in hsq.values) {
  print(paste("Running hsq = ", hsq, sep = ""))
  path = paste(path.tmp, hsq, "/", sep = "")
  nrep = 100
  
  if (hsq == "mid") {
    for (sample_size in sample.sizes) {
      message("Running metrics for sample_size = ", sample_size)
      metrics = compute_metrics(nrep, path, hsq, 0, sample_size)
      
      # Append metrics for each sample_size
      corr = rbind(corr, metrics$pcorr)
      MSE = rbind(MSE, metrics$MSE)
      bias = rbind(bias, metrics$bias)
    }
  } else {
    # For other hsq values
    metrics = compute_metrics(nrep, path, hsq, 0, sample.sizes[1])
    corr = rbind(corr, metrics$pcorr)
    MSE = rbind(MSE, metrics$MSE)
    bias = rbind(bias, metrics$bias)
  }
}
write.table(corr, "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/neutral_evolution/metrics/tPheno_aPS.cor", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(MSE, "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/neutral_evolution/metrics/tPheno_aPS_scaled.mse", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(bias, "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/neutral_evolution/metrics/tPheno_aPS.bias", row.names = FALSE, col.names = TRUE, quote = FALSE)

# QTLs: effect sizes per bin
QTLs.s = data.frame()
for (hsq in hsq.values) {
  message("Running QTLs: effect sizes per bin for hsq = ", hsq)
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
  message("Running QTLs: frequencies per bin for hsq = ", hsq)
  path = paste(path.tmp, hsq, "/", sep = "")
  nrep = 100
  tmp.out = compute_QTLsFreq(nrep, path, gbp = 400, bins = seq(0, 1, 0.001), hsq)
  QTLs.frq = rbind(QTLs.frq, tmp.out)
}
write.table(QTLs.frq, paste(path.tmp, "metrics/tPheno_aPS_0-1-0.001.QTLs.frq", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE)

################################################################################
################################################################################
# Stabilizing Selection
path.tmp = paste(source.path, "stabilizing_selection/", sep = "")
hsq.values = c("total", "mid")
w.values = c("1", "2", "3", "4", "5")
#run_compute_aPS = function(hsq_values, path_tmp, w) {
#  for (hsq in hsq_values) {
#    message("Running hsq = ", hsq, " and w = ", w)
#    path = paste(path_tmp, hsq, "/", w, "/", sep = "")
#    nrep = 100
#    if (hsq == "mid" & (w == 1 | w == 5)) {
#      message("sample size is 100")
#      compute_aPS(nrep, path, 0, 100)
      #message("sample size is 500")
      #compute_aPS(nrep, path, 0, 500)
      #message("sample size is 1000")
      #compute_aPS(nrep, path, 0, 1000)
#    }
#    else{
#      message("sample size is 100")
#      compute_aPS(nrep, path, 0, 100)
#    }
#  }
#}
# aPS
#for (w in w.values) {
#  run_compute_aPS(hsq.values, path.tmp, w)
#}

# compute metrics
corr = data.frame()
MSE = data.frame()
bias = data.frame()
sample.sizes = c(100, 500, 1000)
for (hsq in hsq.values) {
  print(paste("Running hsq = ", hsq, sep = ""))
  corr.tmp = data.frame()
  MSE.tmp = data.frame()
  Bias.tmp = data.frame()
  for (w in w.values) {
    path = paste(path.tmp, hsq, "/", w, "/", sep = "")
    nrep = 100
    if (hsq == "mid" & (w == 1 | w == 5)) {
      for (sample_size in sample.sizes) {
        print(paste("Running sample_size = ", sample_size, sep = ""))
        metrics = compute_metrics(nrep, path, hsq, 0, sample_size)
        
        pcorr.tmp = metrics$pcorr
        pcorr.tmp$w = w
        
        pMSE.tmp = metrics$MSE
        pMSE.tmp$w = w
        
        pBias.tmp = metrics$bias
        pBias.tmp$w = w
        
        corr.tmp = rbind(corr.tmp, pcorr.tmp)
        MSE.tmp = rbind(MSE.tmp, pMSE.tmp)
        Bias.tmp = rbind(Bias.tmp, pBias.tmp)
      }
    } else {
      metrics = compute_metrics(nrep, path, hsq, 0, sample.sizes[1])
      pcorr.tmp = metrics$pcorr
      pcorr.tmp$w = w
      
      pMSE.tmp = metrics$MSE
      pMSE.tmp$w = w
      
      pBias.tmp = metrics$bias
      pBias.tmp$w = w
      
      corr.tmp = rbind(corr.tmp, pcorr.tmp)
      MSE.tmp = rbind(MSE.tmp, pMSE.tmp)
      Bias.tmp = rbind(Bias.tmp, pBias.tmp)
    }
  }
  # append metrics for each sample_size
  corr = rbind(corr, corr.tmp)
  MSE = rbind(MSE, MSE.tmp)
  bias = rbind(bias, Bias.tmp)
}
write.table(corr, "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/metrics/tPheno_aPS.cor", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(MSE, "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/metrics/tPheno_aPS_scaled.mse", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(bias, "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/metrics/tPheno_aPS.bias", row.names = FALSE, col.names = TRUE, quote = FALSE)

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

# for sim replicates w bigger genome size
# aPS
#for (w in c("1", "5")) {
#  run_compute_aPS(c("mid"), 
#                  paste(source.path,
#                        "stabilizing_selection/sims_w_bigger_genome/", sep = ""),
#                  w)
#}
# compute metrics
corr = data.frame()
MSE = data.frame()
bias = data.frame()
sample.sizes = c(100)
for (hsq in c("mid")) {
  print(paste("Running hsq = ", hsq, sep = ""))
  corr.tmp = data.frame()
  MSE.tmp = data.frame()
  Bias.tmp = data.frame()
  for (w in c("1", "5")) {
    path = paste(paste(source.path,
                       "stabilizing_selection/sims_w_bigger_genome/", sep = ""), hsq, "/", w, "/", sep = "")
    nrep = 100
    if (hsq == "mid" & (w == 1 | w == 5)) {
      for (sample_size in sample.sizes) {
        print(paste("Running sample_size = ", sample_size, sep = ""))
        metrics = compute_metrics(nrep, path, hsq, 0, sample_size)
        
        pcorr.tmp = metrics$pcorr
        pcorr.tmp$w = w
        
        pMSE.tmp = metrics$MSE
        pMSE.tmp$w = w
        
        pBias.tmp = metrics$bias
        pBias.tmp$w = w
        
        corr.tmp = rbind(corr.tmp, pcorr.tmp)
        MSE.tmp = rbind(MSE.tmp, pMSE.tmp)
        Bias.tmp = rbind(Bias.tmp, pBias.tmp)
      }
    } else {
      metrics = compute_metrics(nrep, path, hsq, 0, sample.sizes[1])
      pcorr.tmp = metrics$pcorr
      pcorr.tmp$w = w
      
      pMSE.tmp = metrics$MSE
      pMSE.tmp$w = w
      
      pBias.tmp = metrics$bias
      pBias.tmp$w = w
      
      corr.tmp = rbind(corr.tmp, pcorr.tmp)
      MSE.tmp = rbind(MSE.tmp, pMSE.tmp)
      Bias.tmp = rbind(Bias.tmp, pBias.tmp)
    }
  }
  # append metrics for each sample_size
  corr = rbind(corr, corr.tmp)
  MSE = rbind(MSE, MSE.tmp)
  bias = rbind(bias, Bias.tmp)
}
write.table(corr, "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/sims_w_bigger_genome/metrics/tPheno_aPS.cor", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(MSE, "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/sims_w_bigger_genome/metrics/tPheno_aPS_scaled.mse", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(bias, "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/sims_w_bigger_genome/metrics/tPheno_aPS.bias", row.names = FALSE, col.names = TRUE, quote = FALSE)

################################################################################
################################################################################
# Directional Selection
path.tmp = paste(source.path, "directional_selection/", sep = "")
hsq.values = c("total", "mid")
QTLsd.values = c("00025", "0025", "025")
#run_compute_aPS = function(hsq_values, path_tmp, QTLsd) {
# for (hsq in hsq.values) {
#    message("Running hsq = ", hsq, " and QTL sd = ", QTLsd)
#    path = paste(path.tmp, hsq, "/", QTLsd, "/", sep = "")
#    nrep = 100
#    if (hsq == "mid" & QTLsd == "025") {
#      message("sample size is 100")
#      compute_aPS(nrep, path, 1, 100)
#      message("sample size is 500")
#      compute_aPS(nrep, path, 1, 500)
#      message("sample size is 1000")
#      compute_aPS(nrep, path, 1, 1000)
#    }
#    else{
#      message("sample size is 100")
#      compute_aPS(nrep, path, 1, 100)
#    }
#  }
#}
# aPS
#for (QTLsd in QTLsd.values) {
#  run_compute_aPS(hsq.values, path.tmp, QTLsd)
#}

# compute metrics
corr = data.frame()
MSE = data.frame()
bias = data.frame()
sample_sizes = c(100, 500, 1000)
for (hsq in hsq.values) {
  print(paste("Running hsq = ", hsq, sep = ""))
  corr.tmp = data.frame()
  MSE.tmp = data.frame()
  Bias.tmp = data.frame()
  for (QTLsd in QTLsd.values) {
    path = paste(path.tmp, hsq, "/", QTLsd, "/", sep = "")
    nrep = 100
    if (hsq == "mid" & QTLsd == "025") {
      for (sample_size in sample_sizes) {
        print(paste("Running sample_size = ", sample_size, sep = ""))
        metrics = compute_metrics(nrep, path, hsq, 1, sample_size)
        
        pcorr.tmp = metrics$pcorr
        pcorr.tmp$QTLsd = QTLsd
        
        pMSE.tmp = metrics$MSE
        pMSE.tmp$QTLsd = QTLsd
        
        pBias.tmp = metrics$bias
        pBias.tmp$QTLsd = QTLsd
        
        corr.tmp = rbind(corr.tmp, pcorr.tmp)
        MSE.tmp = rbind(MSE.tmp, pMSE.tmp)
        Bias.tmp = rbind(Bias.tmp, pBias.tmp)
      }
    } else {
      metrics = compute_metrics(nrep, path, hsq, 1, sample_sizes[1])
      pcorr.tmp = metrics$pcorr
      pcorr.tmp$QTLsd = QTLsd
      
      pMSE.tmp = metrics$MSE
      pMSE.tmp$QTLsd = QTLsd
      
      pBias.tmp = metrics$bias
      pBias.tmp$QTLsd = QTLsd
      
      corr.tmp = rbind(corr.tmp, pcorr.tmp)
      MSE.tmp = rbind(MSE.tmp, pMSE.tmp)
      Bias.tmp = rbind(Bias.tmp, pBias.tmp)
    }
  }
  # append metrics for each sample_size
  corr = rbind(corr, corr.tmp)
  MSE = rbind(MSE, MSE.tmp)
  bias = rbind(bias, Bias.tmp)
}
write.table(corr, "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/directional_selection/metrics/tPheno_aPS.cor", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(MSE, "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/directional_selection/metrics/tPheno_aPS_scaled.mse", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(bias, "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/directional_selection/metrics/tPheno_aPS.bias", row.names = FALSE, col.names = TRUE, quote = FALSE)

# QTLs: effect sizes per bin
gbp.values = abs(c(200400, 200300, 200200, 200100, 200090, 200080, 200070, 200060, 200050, 200040, 200030, 200020, 200010, 200000) - 200400)
#gbp.values = abs(c(100400, 100300, 100200, 100100, 100000) - 100400)
bins_x_value = data.frame(seq(-0.008, 0.008, 0.001), seq(-0.08, 0.08, 0.01), seq(-0.8, 0.8, 0.1))
colnames(bins_x_value) = QTLsd.values
for (gbp_index in gbp.values) {
  QTLs.s = data.frame()
  for (hsq in hsq.values) {
    print(paste("Running hsq = ", hsq, sep = ""))
    QTLs.s.tmp = data.frame()
    for (QTLsd in QTLsd.values) {
      path = paste(path.tmp, hsq, "/", QTLsd, "/", sep = "")
      nrep = 100
      tmp.out = compute_QTLsDist(nrep, path, gbp = gbp_index, bins = bins_x_value[,QTLsd], hsq)
      tmp.out$QTLsd = QTLsd
      QTLs.s.tmp = rbind(QTLs.s.tmp, tmp.out)
    }
    QTLs.s = rbind(QTLs.s, QTLs.s.tmp)
  }
  QTLs.s$coeff = (QTLs.s$nQTLs / QTLs.s$nQTLs.x.bin)
  QTLs.s$gbp = gbp_index
  write.table(QTLs.s, paste(path.tmp, "metrics/QTLs_summ_0_", gbp_index, ".txt", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE)
}

# stabilizing Selection: realistic parameters for w and hsq
path.tmp = paste(source.path, "/stabilizing_selection/realistic_params/", sep = "")
phenos = c("HT", "BMI")
# aPS
for (pheno in phenos) {
  message("Running aPS for pheno = ", pheno)
  path = paste0(path.tmp, pheno, "/")
  nrep = 100
  # run compute_aPS with base parameters
  compute_aPS(nrep, path, 0, 100)
}

# compute precision metrics: r^2, MSE and bias
corr = data.frame()
MSE = data.frame()
bias = data.frame()
for (pheno in phenos) {
  print(paste("Running pheno = ", pheno, sep = ""))
  path = paste(path.tmp, pheno, "/", sep = "")
  nrep = 100
  metrics = compute_metrics(nrep, path, pheno, 0, 100)
  corr = rbind(corr, metrics$pcorr)
  MSE = rbind(MSE, metrics$MSE)
  bias = rbind(bias, metrics$bias)
}
write.table(corr, "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/realistic_params/metrics/tPheno_aPS.cor", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(MSE, "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/realistic_params/metrics/tPheno_aPS_scaled.mse", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(bias, "/Users/valeriaanorve-garibay/desktop/projects/ancient-pheno-prediction/revisions/output/coalesced_sims/stabilizing_selection/realistic_params/metrics/tPheno_aPS.bias", row.names = FALSE, col.names = TRUE, quote = FALSE)
