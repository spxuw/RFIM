library(ggplot2)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(org.Hs.eg.db)

setwd("/Users/xu-wenwang/Dropbox/Projects/K25/version2.0/code/figures")

gene_wise = c("COPD-pvalue-new.txt", "CVD-pvalue-new.txt","Diabetes-pvalue-new.txt",
              "COPD-pvalue-old.txt", "CVD-pvalue-old.txt","Diabetes-pvalue-old.txt")
phenotype = c('COPD_old','CVD_old','Diabetes_old','COPD_new','CVD_new','Diabetes_new')
method = c('RFIM','DIAMOnD','ModuleDiscoverer','DOMINO','ROBUST','Hierarchical HotNet')

# String
String_comparison = data.frame(matrix(ncol = 4, nrow = 0))
module = list()
for (k in 1:6){
  p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ' ')
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
  module[[(k-1)*6+1]] = module_RFIM$V1
  module[[(k-1)*6+2]] = module_DIAMOND$V1
  module[[(k-1)*6+3]] = module_ModuleDiscoverer$V1
  module[[(k-1)*6+4]] = module_DOMINO$V1
  module[[(k-1)*6+5]] = module_robust$V1
  module[[(k-1)*6+6]] = as.character(module_hierarchical_hotnet[1,])
}
##################################### evaluations #######################
# overlap of disease module
for (i in 1:6){
  for (j in 1:3){
    overlap = length(intersect(module[[(j-1)*6+i]],module[[(j-1)*6+i+18]]))/min(length(module[[(j-1)*6+i]]),length(module[[(j-1)*6+i+18]]))
    String_comparison = rbind(String_comparison,c(method[i],'Overlap',overlap,phenotype[j]))
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
iRefIndex_comparison = data.frame(matrix(ncol = 4, nrow = 0))
module = list()
for (k in 1:6){
  p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ' ')
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
  module[[(k-1)*6+1]] = module_RFIM$V1
  module[[(k-1)*6+2]] = module_DIAMOND$V1
  module[[(k-1)*6+3]] = module_ModuleDiscoverer$V1
  module[[(k-1)*6+4]] = module_DOMINO$V1
  module[[(k-1)*6+5]] = module_robust$V1
  module[[(k-1)*6+6]] = as.character(module_hierarchical_hotnet[1,])
}
##################################### evaluations #######################
# overlap of disease module
for (i in 1:6){
  for (j in 1:3){
    overlap = length(intersect(module[[(j-1)*6+i]],module[[(j-1)*6+i+18]]))/min(length(module[[(j-1)*6+i]]),length(module[[(j-1)*6+i+18]]))
    iRefIndex_comparison = rbind(iRefIndex_comparison,c(method[i],'Overlap',overlap,phenotype[j]))
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
ggsave(p1, file="../../figures/figure7.pdf", width=8, height=8)
