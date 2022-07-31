library(ggplot2)
library(ggpubr)
library(gridExtra)
library(plyr)
library(org.Hs.eg.db)
library(reshape2)
library(enrichR)
library(ReactomePA)
library(caret)

setwd("/Users/xu-wenwang/Dropbox/Projects/K25/version2.0/code/figures")

gene_wise = c('Asthma_combined.csv','Breast_cancer_combined.csv','Lung_cancer_combined.csv',"Colorectal_cancer_combined.csv",
              "Gastric_cancer_combined.csv", "Ovarian_cancer_combined.csv", "Prostate_cancer_combined.csv",
              "COPD-pvalue-new.txt", "CVD-pvalue-new.txt","Diabetes-pvalue-new.txt")
phenotype = c('Asthma','Breast_cancer','Lung_cancer','Colorectal_cancer','Gastric_cancer','Ovarian_cancer',
              'Prostate_cancer','COPD_new','CVD_new','Diabetes_new')
method = c('RFIM','DIAMOnD','ModuleDiscoverer','DOMINO','ROBUST','Hierarchical HotNet')

query_term = c("../../data/Disgenet/C0004096_disease_gda_summary_Asthma.tsv",
               "../../data/Disgenet/C2938924_disease_gda_summary_Breast_cancer.tsv",
               "../../data/Disgenet/C0278504_disease_gda_summary_lung_merged.tsv",
               "../../data/Disgenet/C0009402_disease_gda_summary_Colorectal_Carcinoma.tsv",
               "../../data/Disgenet/C1708349_disease_gda_summary_Hereditary_Diffuse_Gastric_Cancer.tsv",
               "../../data/Disgenet/C0029925_disease_gda_summary_Ovarian_Carcinoma.tsv",
               "../../data/Disgenet/C0600139_disease_gda_summary_Prostate_carcinoma.tsv",
               "../../data/Disgenet/C0024117_disease_gda_summary_COPD.tsv",
               "../../data/Disgenet/C0007222_disease_gda_summary_CVD.tsv",
               "../../data/Disgenet/C0011847_disease_gda_summary_Diabetes.tsv")

# metrics: (1) Number of pathways; (2) Accuracy; (3) Specificity; (4) Sensitivity

# String
ppi_temp = read.csv("../../data/interactome/String_proceed.txt",header = TRUE,sep = '\t')
String_comparison = data.frame(matrix(ncol = 4, nrow = 0))
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
  module = list(RFIM=module_RFIM$V1,DIAMOND=module_DIAMOND$V1,ModuleDiscoverer=module_ModuleDiscoverer$V1,
                DOMINO=module_DOMINO$V1,robust=module_robust$V1,hierarchical_hotnet=as.character(module_hierarchical_hotnet[1,]))
  ##################################### evaluations #######################
  for (i in 1:6){
    # number of pathways
    module$entrez = mapIds(org.Hs.eg.db,
                           keys=module[[i]],
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
    Path_module = enrichPathway(module$entrez, pvalueCutoff=0.05, organism = 'human',pAdjustMethod = "fdr")
    disgene = read.csv(query_term[k],header = TRUE, check.names = FALSE, sep = '\t')
    if (is.null(Path_module)){
      nn = 0
    } else {
      P_1 = Path_module@result$geneID
      nn_2 = 0
      nn_4 = 0
      for (j1 in 1:length(P_1)){
        s1 = strsplit(P_1[j1],'/')
        s1 = as.numeric(unlist(s1))
        s2 = intersect(disgene$Gene_id,s1)
        if (length(s2)>1){nn_2=nn_2+1}
        if (length(s2)>3){nn_4=nn_4+1}
      }
      nn_2 = nn_2 / length(P_1)
    }
    String_comparison = rbind(String_comparison,c(method[i],'Precision (2)',nn_2,phenotype[k]))
    String_comparison = rbind(String_comparison,c(method[i],'Precision (4)',nn_4,phenotype[k]))
    # classification performance
    true_genes = unique(c(ppi_temp$source,ppi_temp$target))
    true_genes[!true_genes%in%disgene$Gene]=0
    true_genes[true_genes%in%disgene$Gene]=1
    predicted_genes = unique(c(ppi_temp$source,ppi_temp$target))
    predicted_genes[!predicted_genes%in%module[[i]]]=0
    predicted_genes[predicted_genes%in%module[[i]]]=1
    
    cm <- confusionMatrix(as.factor(predicted_genes),as.factor(true_genes),positive="1")
    String_comparison = rbind(String_comparison,c(method[i],'Accuracy',as.numeric(cm$overall[1]),phenotype[k]))
    String_comparison = rbind(String_comparison,c(method[i],'Sensitivity',as.numeric(cm$byClass[1]),phenotype[k]))
    String_comparison = rbind(String_comparison,c(method[i],'Specificity',as.numeric(cm$byClass[2]),phenotype[k]))
  }
}
colnames(String_comparison) = c("Method", "Metric", "Performance", "Phenotype")
String_comparison$Method <- factor(String_comparison$Method,levels = c("DIAMOnD","DOMINO","Hierarchical HotNet",'ModuleDiscoverer',"ROBUST","RFIM"))
String_comparison$Performance = as.numeric(String_comparison$Performance)
g1 = ggplot(String_comparison, aes(x=Method, y=Performance,fill=Method)) + geom_bar(stat="identity")+
  facet_grid(Metric~Phenotype,scales = "free")+theme_bw()+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#66a61e","#e6ab02", "#e7298a"))+
  theme(
    line = element_line(size = 0.5), # Thickness of all lines, mm
    rect = element_rect(size = 0.5), # Thickness of all rectangles, mm
    text = element_text(size = 8), # Size of all text, points
    axis.text.x = element_text(size = 8,color = 'black',angle = 90,hjust=0.95,vjust=0.2), # Angle x-axis text
    axis.text.y = element_text(size = 8,color = 'black'), # Angle x-axis text
    axis.title.x = element_text(size = 10), # Size x-axis title
    axis.title.y = element_text(size = 10), # Size y-axis title
    panel.grid.minor = element_blank(), # No minor grid lines
    panel.grid.major.x = element_blank(), # No major x-axis grid lines
    panel.grid.major.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'none',
    legend.box.background = element_rect(size = 0.2), # Box for legend
    legend.key.size = unit(4, unit = 'mm'),
    legend.text = element_text(size = 8),
    strip.background = element_blank()
  )

# iRefIndex
ppi_temp = read.csv("../../data/interactome/iRefIndex_proceed.txt",header = TRUE,sep = '\t')
iRefIndex_comparison = data.frame(matrix(ncol = 4, nrow = 0))
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
  module = list(RFIM=module_RFIM$V1,DIAMOND=module_DIAMOND$V1,ModuleDiscoverer=module_ModuleDiscoverer$V1,
                DOMINO=module_DOMINO$V1,robust=module_robust$V1,hierarchical_hotnet=as.character(module_hierarchical_hotnet[1,]))
  ##################################### evaluations #######################
  for (i in 1:6){
    # number of pathways
    module$entrez = mapIds(org.Hs.eg.db,
                           keys=module[[i]],
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
    Path_module = enrichPathway(module$entrez, pvalueCutoff=0.05, organism = 'human',pAdjustMethod = "fdr")
    disgene = read.csv(query_term[k],header = TRUE, check.names = FALSE, sep = '\t')
    if (is.null(Path_module)){
      nn = 0
    } else {
      P_1 = Path_module@result$geneID
      nn_2 = 0
      nn_4 = 0
      for (j1 in 1:length(P_1)){
        s1 = strsplit(P_1[j1],'/')
        s1 = as.numeric(unlist(s1))
        s2 = intersect(disgene$Gene_id,s1)
        if (length(s2)>1){nn_2=nn_2+1}
        if (length(s2)>3){nn_4=nn_4+1}
      }
      nn_2 = nn_2 / length(P_1)
      nn_4 = nn_4 / length(P_1)
    }
    iRefIndex_comparison = rbind(iRefIndex_comparison,c(method[i],'Precision (2)',nn_2,phenotype[k]))
    iRefIndex_comparison = rbind(iRefIndex_comparison,c(method[i],'Precision (4)',nn_4,phenotype[k]))
    # classification performance
    true_genes = unique(c(ppi_temp$source,ppi_temp$target))
    true_genes[!true_genes%in%disgene$Gene]=0
    true_genes[true_genes%in%disgene$Gene]=1
    predicted_genes = unique(c(ppi_temp$source,ppi_temp$target))
    predicted_genes[!predicted_genes%in%module[[i]]]=0
    predicted_genes[predicted_genes%in%module[[i]]]=1
    
    cm <- confusionMatrix(as.factor(predicted_genes),as.factor(true_genes),positive="1")
    iRefIndex_comparison = rbind(iRefIndex_comparison,c(method[i],'Accuracy',as.numeric(cm$overall[1]),phenotype[k]))
    iRefIndex_comparison = rbind(iRefIndex_comparison,c(method[i],'Sensitivity',as.numeric(cm$byClass[1]),phenotype[k]))
    iRefIndex_comparison = rbind(iRefIndex_comparison,c(method[i],'Specificity',as.numeric(cm$byClass[2]),phenotype[k]))
  }
}
colnames(iRefIndex_comparison) = c("Method", "Metric", "Performance", "Phenotype")
iRefIndex_comparison$Method <- factor(iRefIndex_comparison$Method,levels = c("DIAMOnD","DOMINO","Hierarchical HotNet",'ModuleDiscoverer',"ROBUST","RFIM"))
iRefIndex_comparison$Performance = as.numeric(iRefIndex_comparison$Performance)
g2 = ggplot(iRefIndex_comparison, aes(x=Method, y=Performance,fill=Method)) + geom_bar(stat="identity")+
  facet_grid(Metric~Phenotype,scales = "free")+theme_bw()+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#66a61e","#e6ab02", "#e7298a"))+
  theme(
    line = element_line(size = 0.5), # Thickness of all lines, mm
    rect = element_rect(size = 0.5), # Thickness of all rectangles, mm
    text = element_text(size = 8), # Size of all text, points
    axis.text.x = element_text(size = 8,color = 'black',angle = 90,hjust=0.95,vjust=0.2), # Angle x-axis text
    axis.text.y = element_text(size = 8,color = 'black'), # Angle x-axis text
    axis.title.x = element_text(size = 10), # Size x-axis title
    axis.title.y = element_text(size = 10), # Size y-axis title
    panel.grid.minor = element_blank(), # No minor grid lines
    panel.grid.major.x = element_blank(), # No major x-axis grid lines
    panel.grid.major.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'none',
    legend.box.background = element_rect(size = 0.2), # Box for legend
    legend.key.size = unit(4, unit = 'mm'),
    legend.text = element_text(size = 8),
    strip.background = element_blank()
  )


p1 = grid.arrange(g1,g2,nrow = 2)
ggsave(p1, file="../../figures/figure4.pdf", width=9, height=11)
