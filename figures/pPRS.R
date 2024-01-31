library(vcfR)
library(stringr)
library(hash)
library(tidyr)
library(ggplot2)
options(scipen = 999)

# compute polygenic risk scores using ancient and modern data
# pPRS for individuals of generation t

path.tmp = "/Users/valeriagby/desktop/ancient-pheno-prediction/output/neutral/"
hsq.values = c("low", "mid", "total")
null_model = data.frame()
for (hsq in hsq.values) {
  path = paste(path.tmp, hsq, "/", sep = "")
  nrep = 100
  data = matrix(data = NA, nrow = 100, ncol = 5)
  for (rep in 1:nrep) {
    corrs = c()
    for (gen in c(100400, 100300, 100200, 100100, 100000)) {
      # true pheno for generation t
      tpheno = read.table(paste(path, "tpheno_", gen, "_", rep, ".txt", sep = ""), header = FALSE)
      colnames(tpheno) = "tpheno"
      
      # VCF File for generation t (past)
      pastVCF = read.vcfR(paste(path, "outVCF_", gen, "_", rep, ".vcf", sep = ""), verbose = FALSE)
      # positions and effect sizes for each mutation
      tmp.past = data.frame(do.call(rbind, strsplit(as.character(pastVCF@fix[,"INFO"]), ";")))[2]
      QTLs.past = data.frame("position" = as.integer(pastVCF@fix[,"POS"]), "s" = as.numeric(sub("S=", "", tmp.past$X2)))
      
      # VCF File for generation t = 0 (present-day)
      presentVCF = read.vcfR(paste(path, "outVCF_100400_", rep, ".vcf", sep = ""), verbose = FALSE)
      # QTLs: positions and effect sizes
      tmp.present = data.frame(do.call(rbind, strsplit(as.character(presentVCF@fix[,"INFO"]), ";")))[2]
      QTLs.present = data.frame("position" = as.integer(presentVCF@fix[,"POS"]), "s" = as.numeric(sub("S=", "", tmp.present$X2)))
      # create a dictionarie containing present (t = 0) QTLs effect sizes as values and positions as keys
      QTLs.present.dic = hash(values=as.list(QTLs.present$s), keys = QTLs.present$position)
      
      # compute predicted polygenic risk score (pPRS)
      pPRS = c()
      ind.ids = colnames(pastVCF@gt)[-1]
      for (ind in ind.ids) {
        shared_positions = c()
        shared_positions_index = c()
        for (snp in 1:dim(QTLs.past)[1]) {
          position = QTLs.past$position[snp]
          if(!is.null(QTLs.present.dic[[as.character(position)]])){
            shared_positions = c(shared_positions, position)
            shared_positions_index = c(shared_positions_index, snp)
          }
        }
        geno = pastVCF@gt[shared_positions_index,ind]
        # number of derived alleles
        geno[which(geno == "0|0")] = 0
        geno[which(geno == "1|0")] = 1
        geno[which(geno == "0|1")] = 1
        geno[which(geno == "1|1")] = 2
        # Polygenic Score Prediction
        summatory = 0
        for (snp in 1:length(geno)) {
          summatory = summatory + ((as.integer(geno[snp]))*(QTLs.present.dic[[as.character(shared_positions[snp])]]))
        }
        pPRS = c(pPRS, summatory)
      }
      corrs = c(corrs, cor(pPRS, tpheno$tpheno))
    }
    data[rep,] = corrs
  }
  
  data.new = data.frame(data)
  colnames(data.new) = c("0", "100", "200", "300", "400")
  data.new$rep = c(1:nrep)
  
  data.new.plt = gather(data.new, key = "Variable", value = "Value", -rep)
  data.new.plt = data.new.plt[order(data.new.plt$rep),]
  data.new.plt$Variable = as.numeric(data.new.plt$Variable)
  data.new.plt$hsq = hsq
  
  null_model = rbind(null_model, data.new.plt)
}
