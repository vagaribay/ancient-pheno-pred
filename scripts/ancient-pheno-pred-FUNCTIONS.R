# author: Valeria Añorve-Garibay
# compute ancient-polygenic scores
compute_aPS = function(nrep, path, gen_zoom, present_sample_size){
  if(gen_zoom == 1){
    generations = c(200400, 200300, 200200, 200100,
                    200090, 200080, 200070, 200060, 200050,
                    200040, 200030, 200020, 200010, 200000)
  } else{
    generations = c(200400, 200300, 200200, 200100, 200000)
  }
  
  for (rep in 1:nrep) {
    print(rep)
    for (gen in generations) {
      # true pheno for generation t
      #tpheno = read.table(paste(path, "tpheno_", gen, "_", rep, ".txt", sep = ""), header = FALSE)
      #colnames(tpheno) = "tpheno"
      # VCF File for generation t (past)
      pastVCF = read.vcfR(paste(path, 100, "_noMultiallelicsVCF_", gen, "_", rep, ".vcf", sep = ""), verbose = FALSE)
      #pastVCF = read.vcfR(paste(path, "noMultiallelicsVCF_", gen, "_", rep, ".vcf", sep = ""), verbose = FALSE)
      # positions and effect sizes for each mutation
      tmp.past = data.frame(do.call(rbind, strsplit(as.character(pastVCF@fix[,"INFO"]), ";")))[2]
      QTLs.past = data.frame("position" = as.integer(pastVCF@fix[,"POS"]), "s" = as.numeric(sub("S=", "", tmp.past$X2)))
      
      # VCF File for generation t = 0 (present-day)
      presentVCF = read.vcfR(paste(path, present_sample_size, "_noMultiallelicsVCF_200400_", rep, ".vcf", sep = ""), verbose = FALSE)
      #presentVCF = read.vcfR(paste(path, "noMultiallelicsVCF_100400_", rep, ".vcf", sep = ""), verbose = FALSE)
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
      # aPS
      out.path = paste(path, present_sample_size, "_aPS_", gen, "_", rep, ".txt", sep = "")
      write.table(pPRS, out.path, row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
}
# QTL effect sizes distribution
compute_QTLsDist = function(nrep, path, gbp, bins, hsq){
  pGen = 200400 - gbp
  print(pGen)
  data.conserved = matrix(data = NA, nrow = nrep, ncol = length(bins))
  data.lost = matrix(data = NA, nrow = nrep, ncol = length(bins))
  nQTLs.x.bin.out = c()
  for (rep in 1:nrep) {
    # VCF File for generation 200000 (400 generations ago)
    pastVCF = read.vcfR(paste(path, "100_noMultiallelicsVCF_", pGen, "_", rep, ".vcf", sep = ""), verbose = FALSE)
    # positions and effect sizes for each mutation
    tmp.past = data.frame(do.call(rbind, strsplit(as.character(pastVCF@fix[,"INFO"]), ";")))[2]
    QTLs.past = data.frame("position" = as.integer(pastVCF@fix[,"POS"]), "s" = as.numeric(sub("S=", "", tmp.past$X2)))
    
    # VCF File for generation 200400 (present-day)
    presentVCF = read.vcfR(paste(path, "100_noMultiallelicsVCF_200400_", rep, ".vcf", sep = ""), verbose = FALSE)
    # QTLs: positions and effect sizes
    tmp.present = data.frame(do.call(rbind, strsplit(as.character(presentVCF@fix[,"INFO"]), ";")))[2]
    QTLs.present = data.frame("position" = as.integer(presentVCF@fix[,"POS"]), "s" = as.numeric(sub("S=", "", tmp.present$X2)))
    # create a dictionarie containing present (t = 0) QTLs effect sizes as values and positions as keys
    QTLs.present.dic = hash(values=as.list(QTLs.present$s), keys = QTLs.present$position)
    
    shared_positions_index = c()
    lost_positions_index = c()
    for (snp in 1:dim(QTLs.past)[1]) {
      position = QTLs.past$position[snp]
      if(!is.null(QTLs.present.dic[[as.character(position)]])){
        shared_positions_index = c(shared_positions_index, snp)
      } else{
        lost_positions_index = c(lost_positions_index, snp)
      }
    }
    
    QTLs.past$status = NA
    QTLs.past[shared_positions_index, "status"] = "conserved"
    QTLs.past[lost_positions_index, "status"] = "lost"
    write.table(QTLs.past, paste(path, "QTLs_C_L_", rep, ".txt", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE)

    # effect sizes
    # conserved QTLs
    QTLs.conserved.bins = cut(QTLs.past[which(QTLs.past$status == "conserved"),"s"], breaks = bins)
    nConservedQTLs.x.bin = as.vector(unname(table(QTLs.conserved.bins)))
    data.conserved[rep,] = c(rep, nConservedQTLs.x.bin)
    
    # lost QTLs
    QTLs.lost.bins = cut(QTLs.past[which(QTLs.past$status == "lost"),"s"], breaks = bins)
    nLostQTLs.x.bin = as.vector(unname(table(QTLs.lost.bins)))
    data.lost[rep,] = c(rep, nLostQTLs.x.bin)
    
    nQTLs.x.bin = nConservedQTLs.x.bin + nLostQTLs.x.bin
    nQTLs.x.bin.out = c(nQTLs.x.bin.out, rep(nQTLs.x.bin, 2))
  }
  # out data frame
  data.conserved.new = data.frame(data.conserved)
  data.lost.new = data.frame(data.lost)
  colnames(data.conserved.new) = c("rep", levels(QTLs.conserved.bins))
  colnames(data.lost.new) = c("rep", levels(QTLs.lost.bins))
  data.conserved.new = gather(data.conserved.new, key = "bin", value = "nQTLs", -rep)
  data.lost.new = gather(data.lost.new, key = "bin", value = "nQTLs", -rep)
  data.conserved.new$group = "conserved"
  data.lost.new$group = "lost"
  
  data = rbind(data.conserved.new, data.lost.new)
  data$hsq = hsq
  
  data = data[order(data$rep),]
  data$nQTLs.x.bin = nQTLs.x.bin.out
  
  return(data)
  
}

# QTLs frequencies
compute_QTLsFreq = function(nrep, path, gbp, bins, hsq){
  pGen = 200400 - gbp
  print(pGen)
  data.C = matrix(data = NA, nrow = nrep, ncol = length(bins))
  data.L = matrix(data = NA, nrow = nrep, ncol = length(bins))
  for (rep in 1:nrep) {
    fname = paste(path, "QTLs_C_L_", rep, ".txt", sep = "")
    QTLs = read.table(fname, header = TRUE)
    fname = paste(path, "QTLs_", pGen, "_", rep, ".txt", sep = "")
    QTLs.wpop = read.table(fname, header = TRUE)
    QTLs.wpop$position = QTLs.wpop$position+1
    mdata = merge(QTLs, QTLs.wpop[,c("position", "freq")], by = "position",  all.x = TRUE, all.y = FALSE)
    mdata = mdata[!duplicated(mdata$position),]
    
    C.bins = cut(mdata[which(mdata$status == "conserved"),"freq"], breaks = bins)
    nQTLs.C.x.bin = as.vector(unname(table(C.bins)))
    
    L.bins = cut(mdata[which(mdata$status == "lost"),"freq"], breaks = bins)
    nQTLs.L.x.bin = as.vector(unname(table(L.bins)))
    
    data.C[rep,] = c(rep, nQTLs.C.x.bin)
    data.L[rep,] = c(rep, nQTLs.L.x.bin)
    
  }
  data.C.new = data.frame(data.C)
  data.L.new = data.frame(data.L)
  
  colnames(data.C.new) = c("rep", levels(C.bins))
  colnames(data.L.new) = c("rep", levels(L.bins))
  
  data.C.new.m = gather(data.C.new, key = "bin", value = "nQTLs", -rep)
  data.L.new.m = gather(data.L.new, key = "bin", value = "nQTLs", -rep)
  
  data.C.new.m$group = "conserved"
  data.L.new.m$group = "lost"
  
  data = rbind(data.C.new.m, data.L.new.m)
  data$hsq = hsq
  
  return(data)
}

# compute r^2, MSE and bias between true pheno and ancient polygenic score
compute_metrics = function(nrep, path, hsq, gen_zoom, present_sample_size){
  if(gen_zoom == 1){
    generations = c(200400, 200300, 200200, 200100,
                    200090, 200080, 200070, 200060, 200050,
                    200040, 200030, 200020, 200010, 200000)
  } else{
    generations = c(200400, 200300, 200200, 200100, 200000)
  }
  # correlation between true pheno and aPS
  corr = matrix(data = NA, nrow = nrep, ncol = length(generations))
  # MSE between true pheno and aPS
  MSE = matrix(data = NA, nrow = nrep, ncol = length(generations))
  # bias between true pheno and aPS
  bias = matrix(data = NA, nrow = nrep, ncol = length(generations))
  for (rep in 1:nrep) {
    print(rep)
    corrs = c()
    MSEs = c()
    biases = c()
    for (gen in generations) {
      # true pheno for modern population (generation t = 0)
      #tpheno_modern = read.table(paste(path, "tpheno_", 100400, "_", rep, ".txt", sep = ""), header = FALSE)
      tpheno_modern = read.table(paste(path, 100, "_tpheno_", 200400, "_", rep,".txt", sep = ""), header = FALSE)
      colnames(tpheno_modern) = "tpheno"
      # true pheno for generation t
      tpheno = read.table(paste(path, 100, "_tpheno_", gen, "_", rep, ".txt", sep = ""), header = FALSE)
      colnames(tpheno) = "tpheno"
      # aPS pheno for generation t
      aPS = read.table(paste(path, present_sample_size, "_aPS_", gen, "_", rep, ".txt", sep = ""), header = FALSE)
      colnames(aPS) = "aPS"
      corrs = c(corrs, cor(tpheno$tpheno, aPS$aPS, method = "pearson"))
      MSE_value = mean((tpheno$tpheno - aPS$aPS)^2)
      MSE_scaled = (MSE_value)/(var(tpheno_modern$tpheno))
      MSEs = c(MSEs, MSE_scaled)
      bias_value = mean(aPS$aPS - tpheno$tpheno)
      biases = c(biases, bias_value)
    }
    corr[rep,] = corrs
    MSE[rep,] = MSEs
    bias[rep,] = biases
  }
  corr.new = data.frame(corr)
  colnames(corr.new) = as.character(generations)
  corr.new$rep = c(1:nrep)
  corr.new.out = NA
  corr.new.out = gather(corr.new, key = "generation", value = "corr", -rep)
  corr.new.out = corr.new.out[order(corr.new.out$rep),]
  corr.new.out$hsq = hsq
  corr.new.out$present_sample_size = present_sample_size
  
  MSE.new = data.frame(MSE)
  colnames(MSE.new) = as.character(generations)
  MSE.new$rep = c(1:nrep)
  MSE.new.out = NA
  MSE.new.out = gather(MSE.new, key = "generation", value = "MSE", -rep)
  MSE.new.out = MSE.new.out[order(MSE.new.out$rep),]
  MSE.new.out$hsq = hsq
  MSE.new.out$present_sample_size = present_sample_size
  
  bias.new = data.frame(bias)
  colnames(bias.new) = as.character(generations)
  bias.new$rep = c(1:nrep)
  bias.new.out = NA
  bias.new.out = gather(bias.new, key = "generation", value = "bias", -rep)
  bias.new.out = bias.new.out[order(bias.new.out$rep),]
  bias.new.out$hsq = hsq
  bias.new.out$present_sample_size = present_sample_size
  
  return(list(pcorr = corr.new.out, MSE = MSE.new.out, bias = bias.new.out))
}
