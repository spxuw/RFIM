library(ggplot2)
library(ggpubr)
library(gridExtra)
library(jcolors)
library(plyr)
library(forcats)
library(igraph)
library(org.Hs.eg.db)
library(NetSci)

setwd("/Users/xu-wenwang/Dropbox/Projects/K25/version2.0/code/figures")

gene_wise = c('Asthma_combined.csv','Breast_cancer_combined.csv','Lung_cancer_combined.csv',"Colorectal_cancer_combined.csv",
              "Gastric_cancer_combined.csv", "Ovarian_cancer_combined.csv", "Prostate_cancer_combined.csv",
              "COPD-pvalue-new.txt", "CVD-pvalue-new.txt","Diabetes-pvalue-new.txt")

phenotype = c('Asthma','Breast_cancer','Lung_cancer','Colorectal_cancer','Gastric_cancer','Ovarian_cancer',
              'Prostate_cancer','COPD_new','CVD_new','Diabetes_new')
method = c('RFIM','DIAMOND','ModuleDicoverer','DOMINO','robust')

ppi_temp = read.csv("../../data/interactome/String_proceed.txt",header = TRUE,sep = '\t')
# metrics: (1) p values of genes in disease module; (2) LCC; (3) ratio

# String
p_all = c()
method_p = c()
phenotype_p = c()

LCC = c()
method_LCC = c()
phenotype_LCC = c()

ratio = c()
method_ratio = c()
phenotype_ratio = c()

for (k in 8:8){
  if (k<8){p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ',')}
  if (k>=8){p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ' ')}
  # load the modue of each method
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
  
  module = list(RFIM=module_RFIM$V1,DIAMOND=module_DIAMOND$V1,ModuleDiscoverer=module_ModuleDiscoverer$V1,DOMINO=module_DOMINO$V1,robust=module_robust$V1)
  # p-values of genes in disease module
  for (i in 1:5){
    p_all = c(p_all, p_value$V2[match(module[[i]],p_value$V1)])
    method_p = c(method_p, rep(method[i],length(module[[i]])))
    phenotype_p = c(phenotype_p, rep(phenotype[k],length(module[[i]])))
  }
  # LCC
  graph_module = graph.data.frame(d = ppi_temp, directed = FALSE)
  for (i in 1:5){
    LCC_i = LCC_Significance(N=1000,Targets = module[[i]],G=graph_module)
    LCC = c(LCC,LCC_i$Z)
    method_LCC = c(method_LCC, method[i])
    phenotype_LCC = c(phenotype_LCC, phenotype[k])
  }
  # ratio
  for (i in 1:5){
    graph_module = graph.data.frame(d = ppi_temp, directed = FALSE)
    graph_module_sub = induced_subgraph(graph_module,module[[i]])
    degree_module_genes = degree(graph_module_sub,module[[i]])
    direct_neighbors = c()
    for (j in 1:length(module[[i]])){direct_neighbors = c(direct_neighbors, names(neighbors(graph_module,module[[i]][j])))}
    graph_module_sub = induced_subgraph(graph_module,unique(direct_neighbors))
    ratio = c(ratio, mean(degree_module_genes)/mean(degree(graph_module_sub,unique(direct_neighbors))))
    method_ratio = c(method_ratio, method[i])
    phenotype_ratio = c(phenotype_ratio,phenotype[k])
  }
}
dat_p = data.frame(p_value=-log10(p_all),method=method_p,phenotype=phenotype_p)
dat_LCC = data.frame(z_score=LCC,method=method_LCC,phenotype=phenotype_LCC)
dat_ratio = data.frame(ratio=ratio,method=method_ratio,phenotype=phenotype_ratio)

dat_LCC$method <- factor(dat_LCC$method, levels=c('DIAMOND','ModuleDicoverer','DOMINO','robust','RFIM'))
dat_ratio$method <- factor(dat_ratio$method, levels=c('DIAMOND','ModuleDicoverer','DOMINO','robust','RFIM'))


g1 = ggplot(dat_p, aes(x=p_value, color=method)) + geom_density(alpha=0.4,size=1)+
  facet_wrap(~phenotype,scales = "free",nrow = 1)+theme_bw()+ylab('Density')+xlab("-log10(p-value)")+
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#66a61e", "#e7298a"))+xlim(0,15)+
  theme(
    line = element_line(size = 0.5), # Thickness of all lines, mm
    rect = element_rect(size = 0.5), # Thickness of all rectangles, mm
    text = element_text(size = 8), # Size of all text, points
    axis.text.x = element_text(size = 8,color = 'black'), # Angle x-axis text
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

dat_LCC$z_score[is.na(dat_LCC$z_score)] = 0
g2 = ggplot(dat_LCC, aes(x=method, y= z_score, fill=method)) + geom_bar(stat="identity",color="white")+
  facet_wrap(~phenotype,scales = "free",nrow = 1)+theme_bw()+ylab('z-score')+xlab("")+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#66a61e", "#e7298a"))+
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

g3 = ggplot(dat_ratio, aes(x=method, y= ratio, fill=method)) + geom_bar(stat="identity",color="white")+
  facet_wrap(~phenotype,scales = "free",nrow = 1)+theme_bw()+ylab('Ratio')+xlab("")+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#66a61e", "#e7298a"))+
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

p_all = c()
method_p = c()
phenotype_p = c()

LCC = c()
method_LCC = c()
phenotype_LCC = c()

ratio = c()
method_ratio = c()
phenotype_ratio = c()

for (k in 8:8){
  if (k<8){p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ',')}
  if (k>=8){p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ' ')}
  # load the modue of each method
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
  
  module = list(RFIM=module_RFIM$V1,DIAMOND=module_DIAMOND$V1,ModuleDiscoverer=module_ModuleDiscoverer$V1,DOMINO=module_DOMINO$V1,robust=module_robust$V1)
  # p-values of genes in disease module
  for (i in 1:5){
    p_all = c(p_all, p_value$V2[match(module[[i]],p_value$V1)])
    method_p = c(method_p, rep(method[i],length(module[[i]])))
    phenotype_p = c(phenotype_p, rep(phenotype[k],length(module[[i]])))
  }
  # LCC
  # LCC
  graph_module = graph.data.frame(d = ppi_temp, directed = FALSE)
  for (i in 1:5){
    LCC_i = LCC_Significance(N=1000,Targets = module[[i]],G=graph_module)
    LCC = c(LCC,LCC_i$Z)
    method_LCC = c(method_LCC, method[i])
    phenotype_LCC = c(phenotype_LCC, phenotype[k])
  }
  # ratio
  for (i in 1:5){
    graph_module = graph.data.frame(d = ppi_temp, directed = FALSE)
    graph_module_sub = induced_subgraph(graph_module,module[[i]])
    degree_module_genes = degree(graph_module_sub,module[[i]])
    direct_neighbors = c()
    for (j in 1:length(module[[i]])){direct_neighbors = c(direct_neighbors, names(neighbors(graph_module,module[[i]][j])))}
    graph_module_sub = induced_subgraph(graph_module,unique(direct_neighbors))
    ratio = c(ratio, mean(degree_module_genes)/mean(degree(graph_module_sub,unique(direct_neighbors))))
    method_ratio = c(method_ratio, method[i])
    phenotype_ratio = c(phenotype_ratio,phenotype[k])
  }
}
dat_p = data.frame(p_value=-log10(p_all),method=method_p,phenotype=phenotype_p)
dat_LCC = data.frame(z_score=LCC,method=method_LCC,phenotype=phenotype_LCC)
dat_ratio = data.frame(ratio=ratio,method=method_ratio,phenotype=phenotype_ratio)

dat_LCC$method <- factor(dat_LCC$method, levels=c('DIAMOND','ModuleDicoverer','DOMINO','robust','RFIM'))
dat_ratio$method <- factor(dat_ratio$method, levels=c('DIAMOND','ModuleDicoverer','DOMINO','robust','RFIM'))

g4 = ggplot(dat_p, aes(x=p_value, color=method)) + geom_density(alpha=0.4,size=1)+
  facet_wrap(~phenotype,scales = "free",nrow = 1)+theme_bw()+ylab('Density')+xlab("-log10(p-value)")+
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#66a61e", "#e7298a"))+xlim(0,15)+
  theme(
    line = element_line(size = 0.5), # Thickness of all lines, mm
    rect = element_rect(size = 0.5), # Thickness of all rectangles, mm
    text = element_text(size = 8), # Size of all text, points
    axis.text.x = element_text(size = 8,color = 'black'), # Angle x-axis text
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

dat_LCC$z_score[is.na(dat_LCC$z_score)] = 0
g5 = ggplot(dat_LCC, aes(x=method, y= z_score, fill=method)) + geom_bar(stat="identity",color="white")+
  facet_wrap(~phenotype,scales = "free",nrow = 1)+theme_bw()+ylab('z-score')+xlab("")+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#66a61e", "#e7298a"))+
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

g6 = ggplot(dat_ratio, aes(x=method, y= ratio, fill=method)) + geom_bar(stat="identity",color="white")+
  facet_wrap(~phenotype,scales = "free",nrow = 1)+theme_bw()+ylab('Ratio')+xlab("")+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#66a61e", "#e7298a"))+
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


######################################################
method = c('RFIM','DIAMOND','ModuleDicoverer','DOMINO','robust')
phenotype = ('COPD')

for (k in 1:1){
  # load the modue of each method
  # RFIM
  module_RFIM_1 = read.csv(paste("../../results/RFIM/network-",'COPD_old','_String_f_0.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
  module_RFIM_2 = read.csv(paste("../../results/RFIM/network-",'COPD_new','_String_f_0.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
  
  # DIAMOND
  module_DIAMOND_1 = read.csv(paste("../../results/DIAMOND/DIAMOND_",'COPD_old','_0_string.csv',sep = ""),header = FALSE,sep = ',')
  module_DIAMOND_2 = read.csv(paste("../../results/DIAMOND/DIAMOND_",'COPD_new','_0_string.csv',sep = ""),header = FALSE,sep = ',')
  
  # ModuleDicoverer
  module_ModuleDiscoverer_1 = read.csv(paste("../../results/ModuleDiscoverer/ModuleDiscoverer_",'COPD_old','_0_string.csv',sep = ""),header = FALSE,sep = ',')
  module_ModuleDiscoverer_2 = read.csv(paste("../../results/ModuleDiscoverer/ModuleDiscoverer_",'COPD_new','_0_string.csv',sep = ""),header = FALSE,sep = ',')
  
  # DOMINO
  if ((file.size(paste("../../results/DOMINO/String/0/seed_",'COPD_old','/modules.out',sep = "")))>0){
    module_DOMINO_1 = read.csv(paste("../../results/DOMINO/String/0/seed_",'COPD_old','/modules.out',sep = ""),header = FALSE,sep = ',',na.strings=c("","NA"))
    module_DOMINO_1 = melt(as.matrix(module_DOMINO_1))
    module_DOMINO_1 = module_DOMINO_1[complete.cases(module_DOMINO_1),]
    module_DOMINO_1 = data.frame(V1 = module_DOMINO_1$value)
    module_DOMINO_1$V1 = gsub("\\[","",module_DOMINO_1$V1)
    module_DOMINO_1$V1 = gsub("\\]","",module_DOMINO_1$V1)
    module_DOMINO_1$V1 = gsub(" ","",module_DOMINO_1$V1)
    symbols <- (mapIds(org.Hs.eg.db, keys = module_DOMINO_1$V1, keytype = "ENSEMBL", column="SYMBOL"))
    module_DOMINO_1$V1 = unname(symbols)
  }
  if ((file.size(paste("../../results/DOMINO/String/0/seed_",'COPD_old','/modules.out',sep = "")))==0){
    p_value = read.csv(paste("../../data/gene_wise_p/",'COPD-pvalue-old.txt',sep = ""),header = FALSE,sep = ' ')
    module_DOMINO_1 = p_value[p_value$V2<1e-4,]
  }
  
  if ((file.size(paste("../../results/DOMINO/String/0/seed_",'COPD_new','/modules.out',sep = "")))>0){
    module_DOMINO_2 = read.csv(paste("../../results/DOMINO/String/0/seed_",'COPD_new','/modules.out',sep = ""),header = FALSE,sep = ',',na.strings=c("","NA"))
    module_DOMINO_2 = melt(as.matrix(module_DOMINO_2))
    module_DOMINO_2 = module_DOMINO_2[complete.cases(module_DOMINO_2),]
    module_DOMINO_2 = data.frame(V1 = module_DOMINO_2$value)
    module_DOMINO_2$V1 = gsub("\\[","",module_DOMINO_2$V1)
    module_DOMINO_2$V1 = gsub("\\]","",module_DOMINO_2$V1)
    module_DOMINO_2$V1 = gsub(" ","",module_DOMINO_2$V1)
    symbols <- (mapIds(org.Hs.eg.db, keys = module_DOMINO_2$V1, keytype = "ENSEMBL", column="SYMBOL"))
    module_DOMINO_2$V1 = unname(symbols)
  }
  if ((file.size(paste("../../results/DOMINO/String/0/seed_",'COPD_new','/modules.out',sep = "")))==0){
    p_value = read.csv(paste("../../data/gene_wise_p/",'COPD-pvalue-new.txt',sep = ""),header = FALSE,sep = ' ')
    module_DOMINO_2 = p_value[p_value$V2<1e-4,]
  }
  
  # robust
  module_robust_1 = read.csv(paste("../../results/robust/",'COPD_old','_0_string.csv',sep = ""),header = T,sep = ',')
  module_robust_1 = data.frame(V1 = module_robust_1$vertex)
  module_robust_2 = read.csv(paste("../../results/robust/",'COPD_new','_0_string.csv',sep = ""),header = T,sep = ',')
  module_robust_2 = data.frame(V1 = module_robust_2$vertex)
  
  ovelapping = c(length(intersect(module_RFIM_1$V1,module_RFIM_2$V1))/nrow(module_RFIM_1),
                 length(intersect(module_DIAMOND_1$V1,module_DIAMOND_2$V1))/nrow(module_DIAMOND_1),
                 length(intersect(module_ModuleDiscoverer_1$V1,module_ModuleDiscoverer_2$V1))/nrow(module_ModuleDiscoverer_1),
                 length(intersect(module_DOMINO_1$V1,module_DOMINO_2$V1))/nrow(module_DOMINO_1),
                 length(intersect(module_robust_1$V1,module_robust_2$V1))/nrow(module_robust_1))
}

dat_p = data.frame(Overlap=ovelapping,method=method,phenotype=phenotype)
dat_p$method <- factor(dat_p$method, levels=c('DIAMOND','ModuleDicoverer','DOMINO','robust','RFIM'))


g7 = ggplot(dat_p, aes(x=method, y= Overlap, fill=method)) + geom_bar(stat="identity",color="white")+
  facet_wrap(~phenotype,scales = "free",nrow = 1)+theme_bw()+ylab('Overlap')+xlab("")+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#66a61e", "#e7298a"))+
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

#########
for (k in 1:1){
  if (k<8){p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ',')}
  if (k>=8){p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ' ')}
  # load the modue of each method
  # RFIM
  module_RFIM_1 = read.csv(paste("../../results/RFIM/network-",'COPD_old','_iRefIndex_f_0.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
  module_RFIM_2 = read.csv(paste("../../results/RFIM/network-",'COPD_new','_iRefIndex_f_0.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
  
  # DIAMOND
  module_DIAMOND_1 = read.csv(paste("../../results/DIAMOND/DIAMOND_",'COPD_old','_0_iRefIndex.csv',sep = ""),header = FALSE,sep = ',')
  module_DIAMOND_2 = read.csv(paste("../../results/DIAMOND/DIAMOND_",'COPD_new','_0_iRefIndex.csv',sep = ""),header = FALSE,sep = ',')
  
  # ModuleDicoverer
  module_ModuleDiscoverer_1 = read.csv(paste("../../results/ModuleDiscoverer/ModuleDiscoverer_",'COPD_old','_0_iRefIndex.csv',sep = ""),header = FALSE,sep = ',')
  module_ModuleDiscoverer_2 = read.csv(paste("../../results/ModuleDiscoverer/ModuleDiscoverer_",'COPD_new','_0_iRefIndex.csv',sep = ""),header = FALSE,sep = ',')
  
  # DOMINO
  if ((file.size(paste("../../results/DOMINO/iRefIndex/0/seed_",'COPD_old','/modules.out',sep = "")))>0){
    module_DOMINO_1 = read.csv(paste("../../results/DOMINO/iRefIndex/0/seed_",'COPD_old','/modules.out',sep = ""),header = FALSE,sep = ',',na.strings=c("","NA"))
    module_DOMINO_1 = melt(as.matrix(module_DOMINO_1))
    module_DOMINO_1 = module_DOMINO_1[complete.cases(module_DOMINO_1),]
    module_DOMINO_1 = data.frame(V1 = module_DOMINO_1$value)
    module_DOMINO_1$V1 = gsub("\\[","",module_DOMINO_1$V1)
    module_DOMINO_1$V1 = gsub("\\]","",module_DOMINO_1$V1)
    module_DOMINO_1$V1 = gsub(" ","",module_DOMINO_1$V1)
    symbols <- (mapIds(org.Hs.eg.db, keys = module_DOMINO_1$V1, keytype = "ENSEMBL", column="SYMBOL"))
    module_DOMINO_1$V1 = unname(symbols)
  }
  if ((file.size(paste("../../results/DOMINO/iRefIndex/0/seed_",'COPD_old','/modules.out',sep = "")))==0){
    p_value = read.csv(paste("../../data/gene_wise_p/",'COPD-pvalue-old.txt',sep = ""),header = FALSE,sep = ' ')
    module_DOMINO_1 = p_value[p_value$V2<1e-4,]
  }
  
  if ((file.size(paste("../../results/DOMINO/iRefIndex/0/seed_",'COPD_new','/modules.out',sep = "")))>0){
    module_DOMINO_2 = read.csv(paste("../../results/DOMINO/iRefIndex/0/seed_",'COPD_new','/modules.out',sep = ""),header = FALSE,sep = ',',na.strings=c("","NA"))
    module_DOMINO_2 = melt(as.matrix(module_DOMINO_2))
    module_DOMINO_2 = module_DOMINO_2[complete.cases(module_DOMINO_2),]
    module_DOMINO_2 = data.frame(V1 = module_DOMINO_2$value)
    module_DOMINO_2$V1 = gsub("\\[","",module_DOMINO_2$V1)
    module_DOMINO_2$V1 = gsub("\\]","",module_DOMINO_2$V1)
    module_DOMINO_2$V1 = gsub(" ","",module_DOMINO_2$V1)
    symbols <- (mapIds(org.Hs.eg.db, keys = module_DOMINO_2$V1, keytype = "ENSEMBL", column="SYMBOL"))
    module_DOMINO_2$V1 = unname(symbols)
  }
  if ((file.size(paste("../../results/DOMINO/iRefIndex/0/seed_",'COPD_new','/modules.out',sep = "")))==0){
    p_value = read.csv(paste("../../data/gene_wise_p/",'COPD-pvalue-new.txt',sep = ""),header = FALSE,sep = ' ')
    module_DOMINO_2 = p_value[p_value$V2<1e-4,]
  }
  
  # robust
  module_robust_1 = read.csv(paste("../../results/robust/",'COPD_old','_0_iRefIndex.csv',sep = ""),header = T,sep = ',')
  module_robust_1 = data.frame(V1 = module_robust_1$vertex)
  module_robust_2 = read.csv(paste("../../results/robust/",'COPD_new','_0_iRefIndex.csv',sep = ""),header = T,sep = ',')
  module_robust_2 = data.frame(V1 = module_robust_2$vertex)
  
  ovelapping = c(length(intersect(module_RFIM_1$V1,module_RFIM_2$V1))/nrow(module_RFIM_1),
                 length(intersect(module_DIAMOND_1$V1,module_DIAMOND_2$V1))/nrow(module_DIAMOND_1),
                 length(intersect(module_ModuleDiscoverer_1$V1,module_ModuleDiscoverer_2$V1))/nrow(module_ModuleDiscoverer_1),
                 length(intersect(module_DOMINO_1$V1,module_DOMINO_2$V1))/nrow(module_RFIM_1),
                 length(intersect(module_robust_1$V1,module_robust_2$V1))/nrow(module_robust_1))
}

dat_p = data.frame(Overlap=ovelapping,method=method,phenotype=phenotype)
dat_p$method <- factor(dat_p$method, levels=c('DIAMOND','ModuleDicoverer','DOMINO','robust','RFIM'))

g8 = ggplot(dat_p, aes(x=method, y= Overlap, fill=method)) + geom_bar(stat="identity",color="white")+
  facet_wrap(~phenotype,scales = "free",nrow = 1)+theme_bw()+ylab('Overlap')+xlab("")+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#66a61e", "#e7298a"))+
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


p1 = grid.arrange(g2,g3,g7,g5,g6,g8,nrow = 2)

ggsave(p1, file="../../figures/module_LCC_string_length.pdf", width=5, height=5, dpi = 300, units = "in")
#ggsave(p1, file="../../figures/module_LCC_iRefIndex_length.pdf", width=10, height=6.5, dpi = 300, units = "in")
