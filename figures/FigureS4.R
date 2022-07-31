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

phenotype = c('Asthma','Lung_cancer','Colorectal_cancer','Prostate_cancer','CVD_new','Diabetes_new')
random_type = c('RDPN','rewired','scale_free','shuffled','uniform')

query_term = c("../../data/Disgenet/C0004096_disease_gda_summary_Asthma.tsv",
               "../../data/Disgenet/C0278504_disease_gda_summary_lung_merged.tsv",
               "../../data/Disgenet/C0009402_disease_gda_summary_Colorectal_Carcinoma.tsv",
               "../../data/Disgenet/C0600139_disease_gda_summary_Prostate_carcinoma.tsv",
               "../../data/Disgenet/C0007222_disease_gda_summary_CVD.tsv",
               "../../data/Disgenet/C0011847_disease_gda_summary_Diabetes.tsv")

# metrics: (1) meainingfull enrichment (3) meainingfull overlap 

# String
String_comparison = data.frame(matrix(ncol = 5, nrow = 0))
for (k in 1:6){
  for (i in 1:5){
    # GSEA
    p_values = read.csv(paste("../../results/GSEA_meaningfulness/",phenotype[k],'_',random_type[i],'_String.txt',sep = ""),header = FALSE,sep = ',')
    p_ttest = t.test(p_values$V1[2:11], mu = p_values$V1[1])
    String_comparison = rbind(String_comparison,c(random_type[i],'KEGG enrichment',p_ttest$p.value,phenotype[k],p_values$V1[1]))

    # Overlap
    module_real = read.csv(paste("../../results/RFIM/network-",phenotype[k],'_String_f_0.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
    disgene = read.csv(query_term[k],header = TRUE, check.names = FALSE, sep = '\t')
    overalap_i = c()
    for (j in 0:9){
      module_rand = read.csv(paste("../../results/RFIM/network-",phenotype[k],'_String_',random_type[i],'_',j,'.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
      overalap_i = c(overalap_i,length(intersect(module_rand$V1,disgene$Gene))/min(nrow(module_rand),nrow(disgene)))
    }
    p_ttest = t.test(overalap_i, mu = length(intersect(module_real$V1,disgene$Gene)))
    String_comparison = rbind(String_comparison,c(random_type[i],'DisGeNT overlap',p_ttest$p.value,phenotype[k],length(intersect(module_real$V1,disgene$Gene))/min(nrow(module_real),nrow(disgene))))
  }
}
colnames(String_comparison) = c("Generator", "Metric", "Performance", "Phenotype","Original")
String_comparison$Performance = as.numeric(String_comparison$Performance)
String_comparison$Original = as.numeric(String_comparison$Original)
String_comparison$Performance = -log10((String_comparison$Performance))

g1 = ggplot(String_comparison[String_comparison$Metric=="KEGG enrichment",], aes(x=Generator, y=Performance,size=Original)) + geom_point(color="#984ea3")+
  scale_fill_manual(values = c("#984ea3","#377eb8"))+ylab('-log10(p-value)')+facet_grid(~Phenotype,scales = "free")+
  geom_abline(intercept = -log10(0.05), slope = 0, color="#999999", size=1,alpha=0.5)+theme_bw()+
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
    #legend.position = 'none',
    legend.box.background = element_rect(size = 0.2), # Box for legend
    legend.key.size = unit(4, unit = 'mm'),
    legend.text = element_text(size = 8),
    strip.background = element_blank()
  )
g2 = ggplot(String_comparison[String_comparison$Metric=="DisGeNT overlap",], aes(x=Generator, y=Performance,size=Original)) + geom_point(color="#984ea3")+
  scale_fill_manual(values = c("#984ea3","#377eb8"))+ylab('-log10(p-value)')+facet_grid(~Phenotype,scales = "free")+
  geom_abline(intercept = -log10(0.05), slope = 0, color="#999999", size=1,alpha=0.5)+theme_bw()+
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
    #legend.position = 'none',
    legend.box.background = element_rect(size = 0.2), # Box for legend
    legend.key.size = unit(4, unit = 'mm'),
    legend.text = element_text(size = 8),
    strip.background = element_blank()
  )

# iRefIndex
iRefIndex_comparison = data.frame(matrix(ncol = 5, nrow = 0))
for (k in 1:6){
  for (i in 1:5){
    # GSEA
    p_values = read.csv(paste("../../results/GSEA_meaningfulness/",phenotype[k],'_',random_type[i],'_iRefIndex.txt',sep = ""),header = FALSE,sep = ',')
    p_ttest = t.test(p_values$V1[2:11], mu = p_values$V1[1])
    iRefIndex_comparison = rbind(iRefIndex_comparison,c(random_type[i],'KEGG enrichment',p_ttest$p.value,phenotype[k],p_values$V1[1]))
    
    # Overlap
    module_real = read.csv(paste("../../results/RFIM/network-",phenotype[k],'_iRefIndex_f_0.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
    disgene = read.csv(query_term[k],header = TRUE, check.names = FALSE, sep = '\t')
    overalap_i = c()
    for (j in 0:9){
      module_rand = read.csv(paste("../../results/RFIM/network-",phenotype[k],'_iRefIndex_',random_type[i],'_',j,'.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
      overalap_i = c(overalap_i,length(intersect(module_rand$V1,disgene$Gene))/min(nrow(module_rand),nrow(disgene)))
    }
    p_ttest = t.test(overalap_i, mu = length(intersect(module_real$V1,disgene$Gene)))
    iRefIndex_comparison = rbind(iRefIndex_comparison,c(random_type[i],'DisGeNT overlap',p_ttest$p.value,phenotype[k],length(intersect(module_real$V1,disgene$Gene))/min(nrow(module_real),nrow(disgene))))
  }
}
colnames(iRefIndex_comparison) = c("Generator", "Metric", "Performance", "Phenotype","Original")
iRefIndex_comparison$Performance = as.numeric(iRefIndex_comparison$Performance)
iRefIndex_comparison$Original = as.numeric(iRefIndex_comparison$Original)
iRefIndex_comparison$Performance = -log10((iRefIndex_comparison$Performance))

g3 = ggplot(iRefIndex_comparison[iRefIndex_comparison$Metric=="KEGG enrichment",], aes(x=Generator, y=Performance,size=Original)) + geom_point(color="#984ea3")+
  scale_fill_manual(values = c("#984ea3","#377eb8"))+ylab('-log10(p-value)')+facet_grid(~Phenotype,scales = "free")+
  geom_abline(intercept = -log10(0.05), slope = 0, color="#999999", size=1,alpha=0.5)+theme_bw()+
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
    #legend.position = 'none',
    legend.box.background = element_rect(size = 0.2), # Box for legend
    legend.key.size = unit(4, unit = 'mm'),
    legend.text = element_text(size = 8),
    strip.background = element_blank()
  )
g4 = ggplot(iRefIndex_comparison[iRefIndex_comparison$Metric=="DisGeNT overlap",], aes(x=Generator, y=Performance,size=Original)) + geom_point(color="#984ea3")+
  scale_fill_manual(values = c("#984ea3","#377eb8"))+ylab('-log10(p-value)')+facet_grid(~Phenotype,scales = "free")+
  geom_abline(intercept = -log10(0.05), slope = 0, color="#999999", size=1,alpha=0.5)+theme_bw()+
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
    #legend.position = 'none',
    legend.box.background = element_rect(size = 0.2), # Box for legend
    legend.key.size = unit(4, unit = 'mm'),
    legend.text = element_text(size = 8),
    strip.background = element_blank()
  )


p1 = grid.arrange(g1,g2,g3,g4,nrow = 4)
ggsave(p1, file="../../figures/figureS4.pdf", width=10, height=10)
