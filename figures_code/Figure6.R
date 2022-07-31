library(enrichplot)
library(DOSE)
library(gprofiler2)
library(gridExtra)

setwd("/Users/xu-wenwang/Dropbox/Projects/RFIM/version2.0/code/figures")

gene_wise = c('Asthma_combined.csv','Breast_cancer_combined.csv','Lung_cancer_combined.csv',"Colorectal_cancer_combined.csv",
              "Gastric_cancer_combined.csv", "Ovarian_cancer_combined.csv", "Prostate_cancer_combined.csv",
              "COPD-pvalue-new.txt", "CVD-pvalue-new.txt","Diabetes-pvalue-new.txt")
phenotype = c('Asthma','Breast_cancer','Lung_cancer','Colorectal_cancer','Gastric_cancer','Ovarian_cancer',
              'Prostate_cancer','COPD_new','CVD_new','Diabetes_new')
# String
for (k in 1:10){
  module = read.csv(paste("../../results/RFIM/network-",phenotype[k],'_String_f_0.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
  gostres <- gost(query = module$V1,organism = "hsapiens")
  p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
  pp <- publish_gostplot(p, width = NA, height = NA, filename = NULL )
  
  KEGG = which(gostres$result$source=="KEGG")
  pp2 = publish_gosttable(gostres, highlight_terms = gostres$result[KEGG[1:10],],
                          use_colors = TRUE, 
                          show_columns = c("source", "term_name", "term_size", "intersection_size"),
                          filename = NULL)
  ggsave(pp, file = paste("../../results/Module_network/gp_string_",phenotype[k],'.pdf',sep = ''), width=9, height=4, dpi = 100, units = "in")
  ggsave(pp2, file = paste("../../results/Module_network/gt_string_",phenotype[k],'.pdf',sep = ''), width=9, height=4, dpi = 100, units = "in")
}

# iRefIndex
for (k in 1:10){
  module = read.csv(paste("../../results/RFIM/network-",phenotype[k],'_iRefIndex_f_0.txt-spin1.txt',sep = ""),header = FALSE,sep = ',')
  gostres <- gost(query = module$V1,organism = "hsapiens")
  p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
  pp <- publish_gostplot(p, width = NA, height = NA, filename = NULL )
  
  KEGG = which(gostres$result$source=="KEGG")
  pp2 = publish_gosttable(gostres, highlight_terms = gostres$result[KEGG[1:10],],
                          use_colors = TRUE, 
                          show_columns = c("source", "term_name", "term_size", "intersection_size"),
                          filename = NULL)
  ggsave(pp, file = paste("../../results/Module_network/gp_iRefIndex_",phenotype[k],'.pdf',sep = ''), width=9, height=4, dpi = 100, units = "in")
  ggsave(pp2, file = paste("../../results/Module_network/gt_iRefIndex_",phenotype[k],'.pdf',sep = ''), width=9, height=4, dpi = 100, units = "in")
}
