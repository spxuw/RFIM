library(qapi)
setwd("/Users/xu-wenwang/Dropbox/Projects/K25/version2.0/code")

rm(list = ls())

# dataset
gene_wise = c('Asthma_combined.csv','Breast_cancer_combined.csv','Lung_cancer_combined.csv',"Colorectal_cancer_combined.csv",
              "Gastric_cancer_combined.csv", "Ovarian_cancer_combined.csv", "Prostate_cancer_combined.csv",
              "COPD-pvalue-old.txt", "COPD-pvalue-new.txt", "CVD-pvalue-old.txt","CVD-pvalue-new.txt","Diabetes-pvalue-old.txt", "Diabetes-pvalue-new.txt")
phenotype = c('Asthma','Breast_cancer','Lung_cancer','Colorectal_cancer','Gastric_cancer','Ovarian_cancer',
              'Prostate_cancer','COPD_old','COPD_new','CVD_old','CVD_new','Diabetes_old','Diabetes_new')
# string
ppi_temp = read.csv("../data/interactome/String_proceed.txt",header = TRUE,sep = '\t')
f_remove = c(0,0.1,0.2,0.3,0.4,0.5)
for (k in 1:13){
  for (thre in 1:6){
    if (k<8){p_value = read.csv(paste("../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ',')}
    if (k>=8){p_value = read.csv(paste("../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ' ')}
    p_value = p_value[order(p_value$V2),]
    # remove genes with na p-values and duplicated
    p_value = p_value[complete.cases(p_value),]
    p_value = p_value[!duplicated(p_value$V1),]
    rownames(p_value) = p_value$V1
    
    # choose the overlap proteins
    ppi = ppi_temp
    set.seed(thre)
    if (f_remove[thre]>0){ppi = ppi[-sample(1:nrow(ppi), as.integer(f_remove[thre]*nrow(ppi))),]}
    overlapped_1 = intersect(p_value$V1, unique(c(ppi$source, ppi$target)))
    p_value = p_value[p_value$V1%in%overlapped_1,]
    ppi = ppi[ppi$source%in%overlapped_1,]
    ppi = ppi[ppi$target%in%overlapped_1,]
    
    # remove isolate genes
    p_value = p_value[p_value$V1%in%unique(c(ppi$source,ppi$target)),]
    
    # copy the p-value then revise the local field/p-value
    p_value['p_revised'] = p_value$V2
    for (i in 1:nrow(p_value)){
      m1 = which(ppi$source==p_value$V1[i])
      m2 = which(ppi$target==p_value$V1[i])
      set1 = unique(c(ppi$target[m1],ppi$source[m2]))
      p_set1 = p_value$V2[(match(set1,p_value$V1))]
      if (length(p_set1)>0){p_value$p_revised[i] = min(p_value$V2[i],min(p_set1))}
    }
    p_value['h'] = norminv(1-p_value$V2, mu = 0, sigma = 1)
    p_value$h[p_value$h=='Inf'] = 2147483647
    p_value$h[1] = 2147483647 #set h of the lowest p-value gene as maximum anyway
    
    ppi$edgeweight = p_value[ppi$source,'h'] + p_value[ppi$target,'h']
    ppi$edgeweight[ppi$edgeweight<0] = 0
    node_info = data.frame(nodeindex=0:(nrow(p_value)-1),nodename=p_value$V1,nodeweight=p_value$p_revised)
    ppi$source <- as.character(factor(ppi$source, levels=node_info$nodename, labels=node_info$nodeindex))
    ppi$target <- as.character(factor(ppi$target, levels=node_info$nodename, labels=node_info$nodeindex))
    stats_expression = data.frame(N=nrow(node_info), E=nrow(ppi),pmin=min(node_info$nodeweight),
                                  p1min=min(node_info$nodeweight),pmax=max(node_info$nodeweight),
                                  wmin=min(ppi$edgeweight),wmax=max(ppi$edgeweight),wave=mean(ppi$edgeweight))
    NAMES = colnames(stats_expression)
    NAMES[1] <- paste0("#", NAMES[1])
    write.table(stats_expression, paste("../data_RFIM/",phenotype[k],"_String_f_",f_remove[thre],".txt",sep = ""), col.names=NAMES, sep=" ",quote = FALSE, row.names = FALSE)
    NAMES = colnames(node_info)
    NAMES[1] <- paste0("#", NAMES[1])
    write.table(node_info, paste("../data_RFIM/",phenotype[k],"_String_f_",f_remove[thre],".txt",sep = ""), col.names=NAMES, sep=" ",quote = FALSE, row.names = FALSE,append = TRUE)
    NAMES = colnames(ppi)
    NAMES[1] <- paste0("#", NAMES[1])
    write.table(ppi, paste("../data_RFIM/",phenotype[k],"_String_f_",f_remove[thre],".txt",sep = ""), col.names=NAMES, sep=" ",quote = FALSE, row.names = FALSE,append = TRUE)
  }
}

# ....................iRefIndex.....................................................
ppi_temp = read.csv("../data/Interactome/iRefIndex_proceed.txt",header = TRUE,sep = '\t')
f_remove = c(0,0.1,0.2,0.3,0.4,0.5)
for (k in 1:13){
  for (thre in 1:6){
    if (k<8){p_value = read.csv(paste("../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ',')}
    if (k>=8){p_value = read.csv(paste("../data/gene_wise_p/",gene_wise[k],sep = ""),header = FALSE,sep = ' ')}
    p_value = p_value[order(p_value$V2),]
    # remove genes with na p-values and duplicated
    p_value = p_value[complete.cases(p_value),]
    p_value = p_value[!duplicated(p_value$V1),]
    rownames(p_value) = p_value$V1
    
    # choose the overlap proteins
    ppi = ppi_temp
    set.seed(thre)
    if (f_remove[thre]>0){ppi = ppi[-sample(1:nrow(ppi), as.integer(f_remove[thre]*nrow(ppi))),]}
    overlapped_1 = intersect(p_value$V1, unique(c(ppi$source, ppi$target)))
    p_value = p_value[p_value$V1%in%overlapped_1,]
    ppi = ppi[ppi$source%in%overlapped_1,]
    ppi = ppi[ppi$target%in%overlapped_1,]
    
    # remove isolate genes
    p_value = p_value[p_value$V1%in%unique(c(ppi$source,ppi$target)),]
    
    # copy the p-value then revise the local field/p-value
    p_value['p_revised'] = p_value$V2
    for (i in 1:nrow(p_value)){
      m1 = which(ppi$source==p_value$V1[i])
      m2 = which(ppi$target==p_value$V1[i])
      set1 = unique(c(ppi$target[m1],ppi$source[m2]))
      p_set1 = p_value$V2[(match(set1,p_value$V1))]
      if (length(p_set1)>0){p_value$p_revised[i] = min(p_value$V2[i],min(p_set1))}
    }
    p_value['h'] = norminv(1-p_value$V2, mu = 0, sigma = 1)
    p_value$h[p_value$h=='Inf'] = 2147483647
    p_value$h[1] = 2147483647 #set h of the lowest p-value gene as maximum anyway
    
    ppi$edgeweight = p_value[ppi$source,'h'] + p_value[ppi$target,'h']
    ppi$edgeweight[ppi$edgeweight<0] = 0
    node_info = data.frame(nodeindex=0:(nrow(p_value)-1),nodename=p_value$V1,nodeweight=p_value$p_revised)
    ppi$source <- as.character(factor(ppi$source, levels=node_info$nodename, labels=node_info$nodeindex))
    ppi$target <- as.character(factor(ppi$target, levels=node_info$nodename, labels=node_info$nodeindex))
    stats_expression = data.frame(N=nrow(node_info), E=nrow(ppi),pmin=min(node_info$nodeweight),
                                  p1min=min(node_info$nodeweight),pmax=max(node_info$nodeweight),
                                  wmin=min(ppi$edgeweight),wmax=max(ppi$edgeweight),wave=mean(ppi$edgeweight))
    NAMES = colnames(stats_expression)
    NAMES[1] <- paste0("#", NAMES[1])
    write.table(stats_expression, paste("../data_RFIM/",phenotype[k],"_iRefIndex_f_",f_remove[thre],".txt",sep = ""), col.names=NAMES, sep=" ",quote = FALSE, row.names = FALSE)
    NAMES = colnames(node_info)
    NAMES[1] <- paste0("#", NAMES[1])
    write.table(node_info, paste("../data_RFIM/",phenotype[k],"_iRefIndex_f_",f_remove[thre],".txt",sep = ""), col.names=NAMES, sep=" ",quote = FALSE, row.names = FALSE,append = TRUE)
    NAMES = colnames(ppi)
    NAMES[1] <- paste0("#", NAMES[1])
    write.table(ppi, paste("../data_RFIM/",phenotype[k],"_iRefIndex_f_",f_remove[thre],".txt",sep = ""), col.names=NAMES, sep=" ",quote = FALSE, row.names = FALSE,append = TRUE)
  }
}