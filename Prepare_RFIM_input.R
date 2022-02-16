library(qapi)
setwd("/path of code")

rm(list = ls())

ppi_temp = read.csv("../data/Interactome/9606.protein.links.v11.5.txt",header = TRUE,sep = ' ')
mapping_ppi = read.csv("../data/Interactome/9606.protein.info.v11.5.txt",header = TRUE, sep = '\t')
ppi_temp$combined_score = ppi_temp$combined_score/max(ppi_temp$combined_score)
ppi_temp = ppi_temp[ppi_temp$combined_score>0.7,]
colnames(ppi_temp) = c('source','target','edgeweight')
ppi_temp$source <- as.character(factor(ppi_temp$source, levels=mapping_ppi$string_protein_id, labels=mapping_ppi$preferred_name))
ppi_temp$target <- as.character(factor(ppi_temp$target, levels=mapping_ppi$string_protein_id, labels=mapping_ppi$preferred_name))
ppi_temp = ppi_temp[!duplicated(data.frame(t(apply(ppi_temp[,1:2],1,sort)))),]


phenotype = 'Asthma'

ff = c(0,0.1,0.2,0.3,0.4,0.5) # incompleteness of interactome
# p-values association
for (thre in 1:6){
  p_old = read.csv(paste("../data/",phenotype,"/.Pvalue.csv",sep = ""),header = FALSE,sep = ',')
  p_old = p_old[complete.cases(p_old), ]
  p_old = p_old[!duplicated(p_old$V1),]
  
  ppi = ppi_temp
  set.seed(thre)
  if (ff[thre]>0){ppi = ppi[-sample(1:nrow(ppi), as.integer(ff[thre]*nrow(ppi))),]}
  overlapped_1 = intersect(p_old$V1, unique(c(ppi$source, ppi$target)))
  p_old = p_old[p_old$V1%in%overlapped_1,]
  ppi = ppi[ppi$source%in%overlapped_1,]
  ppi = ppi[ppi$target%in%overlapped_1,]
  p_old = p_old[p_old$V1%in%unique(c(ppi$source,ppi$target)),]
  
  p_old['V4'] = p_old$V2
  for (i in 1:nrow(p_old)){
    m1 = which(ppi$source==p_old$V1[i])
    m2 = which(ppi$target==p_old$V1[i])
    set1 = unique(c(ppi$target[m1],ppi$source[m2]))
    p_set1 = p_old$V2[(match(set1,p_old$V1))]
    if (length(p_set1)>1){p_old$V4[i] = min(p_old$V2[i],min(p_set1))}
  }
  p_old['V3'] = norminv(1-p_old$V2, mu = 0, sigma = 1)
  p_old$V3[p_old$V3=='Inf'] = 2147483647
  p_old$V3[1] = 2147483647
  
  for (i in 1:nrow(ppi)){
    m1 = which(p_old$V1==ppi$source[i])
    m2 = which(p_old$V1==ppi$target[i])
    ppi$edgeweight[i] = max(0, p_old$V3[m1]+p_old$V3[m2])
  }
  
  node_info = data.frame(nodeindex=0:(nrow(p_old)-1),nodename=p_old$V1,nodeweight=p_old$V4)
  ppi$source <- as.character(factor(ppi$source, levels=node_info$nodename, labels=node_info$nodeindex))
  ppi$target <- as.character(factor(ppi$target, levels=node_info$nodename, labels=node_info$nodeindex))
  stats_expression = data.frame(N=nrow(node_info), E=nrow(ppi),pmin=min(node_info$nodeweight),
                                p1min=min(node_info$nodeweight),pmax=max(node_info$nodeweight),
                                wmin=min(ppi$edgeweight),wmax=max(ppi$edgeweight),wave=mean(ppi$edgeweight))
  NAMES = colnames(stats_expression)
  NAMES[1] <- paste0("#", NAMES[1])
  write.table(stats_expression, paste("../data/",phenotype,"_string_f_",ff[thre],".txt",sep = ""), col.names=NAMES, sep=" ",quote = FALSE, row.names = FALSE)
  NAMES = colnames(node_info)
  NAMES[1] <- paste0("#", NAMES[1])
  write.table(node_info, paste("../data/",phenotype,"_string_f_",ff[thre],".txt",sep = ""), col.names=NAMES, sep=" ",quote = FALSE, row.names = FALSE,append = TRUE)
  NAMES = colnames(ppi)
  NAMES[1] <- paste0("#", NAMES[1])
  write.table(ppi, paste("../data/",phenotype,"_string_f_",ff[thre],".txt",sep = ""), col.names=NAMES, sep=" ",quote = FALSE, row.names = FALSE,append = TRUE)
}

