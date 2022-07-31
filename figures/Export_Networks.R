setwd("/Users/xu-wenwang/Dropbox/Projects/K25/version2.0/code/figures")

gene_wise = c('Asthma_combined.csv','Breast_cancer_combined.csv','Lung_cancer_combined.csv',"Colorectal_cancer_combined.csv",
              "Gastric_cancer_combined.csv", "Ovarian_cancer_combined.csv", "Prostate_cancer_combined.csv",
              "COPD-pvalue-new.txt", "CVD-pvalue-new.txt","Diabetes-pvalue-new.txt")
phenotype = c('Asthma','Breast_cancer','Lung_cancer','Colorectal_cancer','Gastric_cancer','Ovarian_cancer',
              'Prostate_cancer','COPD_new','CVD_new','Diabetes_new')
# String
ppi_temp = read.csv("../../data/interactome/String_proceed.txt",header = TRUE,sep = '\t')
for (k in 1:10){
  if (k<8){p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ',')}
  if (k>=8){p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ' ')}
  module_RFIM = read.csv(paste("../../results/RFIM/network-",phenotype[k],'_String_f_0.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
  
  ppi = ppi_temp
  ppi = ppi[ppi$source%in%module_RFIM$V1&ppi$target%in%module_RFIM$V1,]
  ppi$edgeweight = 1
  unique_node = unique(c(ppi$source,ppi$target))
  node_info = data.frame(id=unique_node,weight=p_value$V2[match(unique_node,p_value$V1)])
  node_info$weight=node_info$weight<0.001
  
  write.table(ppi, file = paste("../../results/Module_network/ppi_string_",phenotype[k],'.csv',sep = ''), row.names = FALSE, col.names = TRUE, sep=",")
  write.table(node_info, file = paste("../../results/Module_network/node_string_",phenotype[k],'.csv',sep = ''), row.names = FALSE, col.names = TRUE, sep=",")
}

# iRefIndex
ppi_temp = read.csv("../../data/interactome/iRefIndex_proceed.txt",header = TRUE,sep = '\t')
for (k in 1:10){
  if (k<8){p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ',')}
  if (k>=8){p_value = read.csv(paste("../../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ' ')}
  module_RFIM = read.csv(paste("../../results/RFIM/network-",phenotype[k],'_iRefIndex_f_0.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
  
  ppi = ppi_temp
  ppi = ppi[ppi$source%in%module_RFIM$V1&ppi$target%in%module_RFIM$V1,]
  ppi$edgeweight = 1
  unique_node = unique(c(ppi$source,ppi$target))
  node_info = data.frame(id=unique_node,weight=p_value$V2[match(unique_node,p_value$V1)])
  node_info$weight=node_info$weight<0.001
  
  write.table(ppi, file = paste("../../results/Module_network/ppi_iRefIndex_",phenotype[k],'.csv',sep = ''), row.names = FALSE, col.names = TRUE, sep=",")
  write.table(node_info, file = paste("../../results/Module_network/node_iRefIndex_",phenotype[k],'.csv',sep = ''), row.names = FALSE, col.names = TRUE, sep=",")
}
