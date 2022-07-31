library(ggplot2)
library(ggpubr)
library(gridExtra)
library(forcats)
library(org.Hs.eg.db)
library(reshape2)

setwd("/Users/xu-wenwang/Dropbox/Projects/K25/version2.0/code/figures")

gene_wise = c('Asthma_combined.csv','Breast_cancer_combined.csv','Lung_cancer_combined.csv',"Colorectal_cancer_combined.csv",
              "Gastric_cancer_combined.csv", "Ovarian_cancer_combined.csv", "Prostate_cancer_combined.csv",
              "COPD-pvalue-new.txt", "CVD-pvalue-new.txt","Diabetes-pvalue-new.txt")
phenotype = c('Asthma','Breast_cancer','Lung_cancer','Colorectal_cancer','Gastric_cancer','Ovarian_cancer',
              'Prostate_cancer','COPD_new','CVD_new','Diabetes_new')
method = c('RFIM','DIAMOnD','ModuleDiscoverer','DOMINO','ROBUST','Hierarchical HotNet')
# String
summary_string = data.frame()
for (k in 1:10){
  if (k<8){p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ',')}
  if (k>=8){p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ' ')}
  # load the module of each method
  # RFIM
  module_RFIM = read.csv(paste("../../results/RFIM/network-",phenotype[k],'_String_f_0.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
  # DIAMOND
  module_DIAMOND = read.csv(paste("../../results/DIAMOND/DIAMOND_",phenotype[k],'_0_string.csv',sep = ""),header = FALSE,sep = ',')
  # ModuleDicoverer
  module_ModuleDiscoverer = read.csv(paste("../../results/ModuleDiscoverer/ModuleDiscoverer_",phenotype[k],'_0_string.csv',sep = ""),header = FALSE,sep = ',')
  # DOMINO
  if ((file.size(paste("../../results/DOMINO/String/0/seed_",phenotype[k],'/modules.out',sep = "")))>0){
    module_DOMINO = read.csv(paste("../../results/DOMINO/String/0/seed_",phenotype[k],'/modules.out',sep = ""),header = FALSE,sep = ',',na.strings=c("","NA"))
    module_DOMINO = melt(as.matrix(module_DOMINO))
    module_DOMINO = module_DOMINO[complete.cases(module_DOMINO),]
    module_DOMINO = data.frame(V1 = module_DOMINO$value)
    module_DOMINO$V1 = gsub("\\[","",module_DOMINO$V1)
    module_DOMINO$V1 = gsub("\\]","",module_DOMINO$V1)
    module_DOMINO$V1 = gsub(" ","",module_DOMINO$V1)
    symbols <- (mapIds(org.Hs.eg.db, keys = module_DOMINO$V1, keytype = "ENSEMBL", column="SYMBOL"))
    module_DOMINO$V1 = unname(symbols)
  }
  if ((file.size(paste("../../results/DOMINO/String/0/seed_",phenotype[k],'/modules.out',sep = "")))==0){
    module_DOMINO = p_value[p_value$V2<1e-4,]
  }
  # robust
  module_robust = read.csv(paste("../../results/robust/",phenotype[k],'_0_string.csv',sep = ""),header = T,sep = ',')
  module_robust = data.frame(V1 = module_robust$vertex)
  
  # hierarchical-hotnet
  module_hierarchical_hotnet = read.csv(paste("../../results/hierarchical-hotnet/clusters_",phenotype[k],'_0_score_string.tsv',sep = ""),header = F,sep = '\t',nrows = 1,skip = 7)
  module_sub = c(module_RFIM$V1,module_DIAMOND$V1,module_ModuleDiscoverer$V1,module_DOMINO$V1,module_robust$V1,as.character(module_hierarchical_hotnet[1,]))
  module_sub = as.data.frame(table(module_sub))
  module_sub['phenotype'] = phenotype[k]
  module_sub['p-value'] = p_value$V2[match(module_sub$module_sub,p_value$V1)]
  summary_string = rbind(summary_string,module_sub)
}

# iRefIndex
summary_iRefIndex = data.frame()
for (k in 1:10){
  if (k<8){p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ',')}
  if (k>=8){p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ' ')}
  # load the module of each method
  # RFIM
  module_RFIM = read.csv(paste("../../results/RFIM/network-",phenotype[k],'_iRefIndex_f_0.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
  # DIAMOND
  module_DIAMOND = read.csv(paste("../../results/DIAMOND/DIAMOND_",phenotype[k],'_0_iRefIndex.csv',sep = ""),header = FALSE,sep = ',')
  # ModuleDicoverer
  module_ModuleDiscoverer = read.csv(paste("../../results/ModuleDiscoverer/ModuleDiscoverer_",phenotype[k],'_0_iRefIndex.csv',sep = ""),header = FALSE,sep = ',')
  # DOMINO
  if ((file.size(paste("../../results/DOMINO/iRefIndex/0/seed_",phenotype[k],'/modules.out',sep = "")))>0){
    module_DOMINO = read.csv(paste("../../results/DOMINO/iRefIndex/0/seed_",phenotype[k],'/modules.out',sep = ""),header = FALSE,sep = ',',na.strings=c("","NA"))
    module_DOMINO = melt(as.matrix(module_DOMINO))
    module_DOMINO = module_DOMINO[complete.cases(module_DOMINO),]
    module_DOMINO = data.frame(V1 = module_DOMINO$value)
    module_DOMINO$V1 = gsub("\\[","",module_DOMINO$V1)
    module_DOMINO$V1 = gsub("\\]","",module_DOMINO$V1)
    module_DOMINO$V1 = gsub(" ","",module_DOMINO$V1)
    symbols <- (mapIds(org.Hs.eg.db, keys = module_DOMINO$V1, keytype = "ENSEMBL", column="SYMBOL"))
    module_DOMINO$V1 = unname(symbols)
  }
  if ((file.size(paste("../../results/DOMINO/iRefIndex/0/seed_",phenotype[k],'/modules.out',sep = "")))==0){
    module_DOMINO = p_value[p_value$V2<1e-4,]
  }
  # robust
  module_robust = read.csv(paste("../../results/robust/",phenotype[k],'_0_iRefIndex.csv',sep = ""),header = T,sep = ',')
  module_robust = data.frame(V1 = module_robust$vertex)
  
  # hierarchical-hotnet
  module_hierarchical_hotnet = read.csv(paste("../../results/hierarchical-hotnet/clusters_",phenotype[k],'_0_score_iRefIndex.tsv',sep = ""),header = F,sep = '\t',nrows = 1,skip = 7)
  module_sub = c(module_RFIM$V1,module_DIAMOND$V1,module_ModuleDiscoverer$V1,module_DOMINO$V1,module_robust$V1,as.character(module_hierarchical_hotnet[1,]))
  module_sub = as.data.frame(table(module_sub))
  module_sub['phenotype'] = phenotype[k]
  module_sub['p-value'] = p_value$V2[match(module_sub$module_sub,p_value$V1)]
  summary_iRefIndex = rbind(summary_iRefIndex,module_sub)
}
summary_string_4 = summary_string[summary_string$Freq>=4,]
summary_iRefIndex_4 = summary_iRefIndex[summary_iRefIndex$Freq>=4,]
summary_merged_4 = rbind(summary_string_4,summary_iRefIndex_4)
summary_merged_4 = summary_merged_4[!duplicated(summary_merged_4[c(1,3)]),]
summary_merged_4 = summary_merged_4[order(summary_merged_4$phenotype,summary_merged_4$`p-value`),]
write.table(summary_merged_4,file = "../../results/consensus.csv",sep = ",")
