library(ggplot2)
library(ggpubr)
library(gridExtra)

setwd("/Users/xu-wenwang/Dropbox/Projects/RFIM/version2.0/code/figures")

phenotype = c('Asthma','Lung_cancer','Colorectal_cancer','Prostate_cancer','CVD_new','Diabetes_new')
method = c('RFIM','DIAMOnD','ModuleDiscoverer','DOMINO','ROBUST','Hierarchical HotNet')
# metrics: (1) KEGG p-value

# String
String_comparison = data.frame(matrix(ncol = 4, nrow = 0))
for (k in 1:6){
  for (i in 1:6){
    if (i==1){
      p_values = read.csv(paste("../../results/GSEA_meaningfulness/",phenotype[k],'_RDPN_String.txt',sep = ""),header = F)
    } else {
      p_values = read.csv(paste("../../results/GSEA_meaningfulness/",method[i],'_',phenotype[k],'_String.txt',sep = ""),header = F)
    }
    String_comparison = rbind(String_comparison,c(method[i],'KEGG enrichment',p_values$V1[1],phenotype[k]))
  }
}
colnames(String_comparison) = c("Method", "Metric", "Performance", "Phenotype")
String_comparison$Method <- factor(String_comparison$Method,levels = c("DIAMOnD","DOMINO","Hierarchical HotNet",'ModuleDiscoverer',"ROBUST","RFIM"))
String_comparison$Performance = as.numeric(String_comparison$Performance)
g1 = ggplot(String_comparison, aes(x=Method, y=Performance,fill=Method)) + geom_bar(stat="identity")+
  facet_grid(Metric~Phenotype,scales = "free")+theme_bw()+ylab('-log10(p-value)')+
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
for (k in 1:6){
  for (i in 1:6){
    if (i==1){
      p_values = read.csv(paste("../../results/GSEA_meaningfulness/",phenotype[k],'_RDPN_iRefIndex.txt',sep = ""),header = F)
    } else {
      p_values = read.csv(paste("../../results/GSEA_meaningfulness/",method[i],'_',phenotype[k],'_iRefIndex.txt',sep = ""),header = F)
    }
    iRefIndex_comparison = rbind(iRefIndex_comparison,c(method[i],'KEGG enrichment',p_values$V1[1],phenotype[k]))
  }
}
colnames(iRefIndex_comparison) = c("Method", "Metric", "Performance", "Phenotype")
iRefIndex_comparison$Method <- factor(iRefIndex_comparison$Method,levels = c("DIAMOnD","DOMINO","Hierarchical HotNet",'ModuleDiscoverer',"ROBUST","RFIM"))
iRefIndex_comparison$Performance = as.numeric(iRefIndex_comparison$Performance)
g2 = ggplot(iRefIndex_comparison, aes(x=Method, y=Performance,fill=Method)) + geom_bar(stat="identity")+
  facet_grid(Metric~Phenotype,scales = "free")+theme_bw()+ylab('-log10(p-value)')+
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
ggsave(p1, file="../../figures/figureS2.pdf", width=7, height=5)
