library(ggplot2)
library(gridExtra)
library(reshape2)

setwd("/Users/xu-wenwang/Dropbox/Projects/K25/version2.0/code/figures")

phenotype = c('Asthma')

# String
String_comparison = matrix(0,10,10)
for (i in 1:10){
  for (j in 1:10){
    # Overlap
    module_i = read.csv(paste("../../results/RFIM/network-",phenotype,'_String_store_rea_',i,'.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
    module_j = read.csv(paste("../../results/RFIM/network-",phenotype,'_String_store_rea_',j,'.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
    String_comparison[i,j] = length(intersect(module_i$V1,module_j$V1))/min(nrow(module_i),nrow(module_j))
  }
}
diag(String_comparison) = NA
String_comparison = melt(String_comparison)
colnames(String_comparison) = c("realization_1", "realization_2", "Overlap")

g1 = ggplot(String_comparison, aes(x=factor(realization_1), y=factor(realization_2),fill=Overlap)) + 
  geom_tile(color='white')+scale_fill_gradientn(colours=c("#0571b0","#f7f7f7","#ca0020"))+
  theme_bw()+xlab("Realization")+ylab("Realization")+
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
    #legend.position = 'none',
    legend.box.background = element_rect(size = 0.2), # Box for legend
    legend.key.size = unit(4, unit = 'mm'),
    legend.text = element_text(size = 8),
    strip.background = element_blank()
  )

# iRefIndex
iRefIndex_comparison = matrix(0,10,10)
for (i in 1:10){
  for (j in 1:10){
    # Overlap
    module_i = read.csv(paste("../../results/RFIM/network-",phenotype,'_iRefIndex_store_rea_',i,'.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
    module_j = read.csv(paste("../../results/RFIM/network-",phenotype,'_iRefIndex_store_rea_',j,'.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
    iRefIndex_comparison[i,j] = length(intersect(module_i$V1,module_j$V1))/min(nrow(module_i),nrow(module_j))
  }
}
diag(iRefIndex_comparison) = NA
iRefIndex_comparison = melt(iRefIndex_comparison)
colnames(iRefIndex_comparison) = c("realization_1", "realization_2", "Overlap")

g2 = ggplot(iRefIndex_comparison, aes(x=factor(realization_1), y=factor(realization_2),fill=Overlap)) + 
  geom_tile(color='white')+scale_fill_gradientn(colours=c("#0571b0","#f7f7f7","#ca0020"))+
  theme_bw()+xlab("Realization")+ylab("Realization")+
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
    #legend.position = 'none',
    legend.box.background = element_rect(size = 0.2), # Box for legend
    legend.key.size = unit(4, unit = 'mm'),
    legend.text = element_text(size = 8),
    strip.background = element_blank()
  )

p1 = grid.arrange(g1,g2,nrow = 1)
ggsave(p1, file="../../figures/figure_order.pdf", width=8, height=3)
