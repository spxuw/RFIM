library(ggplot2)
library(gridExtra)
library(scales)

setwd("/Users/xu-wenwang/Dropbox/Projects/RFIM/version2.0/code/figures")

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# DIAMOND
Time_diamond_string = read.csv("../../results/DIAMOND/DIAMOND_time_string.csv",header = FALSE)
Time_diamond_iRefIndex = read.csv("../../results/DIAMOND/DIAMOND_time_iRefIndex.csv",header = FALSE)
Time_diamond_string = Time_diamond_string[seq(1,78,6),]
Time_diamond_iRefIndex = Time_diamond_iRefIndex[seq(1,78,6),]

# ModuleDiscoveer
Time_modulediscover_string = read.csv("../../results/ModuleDiscoverer/ModuleDiscoverer_time_string_v1.csv",header = FALSE)
Time_modulediscover_iRefIndex = read.csv("../../results/ModuleDiscoverer/ModuleDiscoverer_time_iRefIndex_v1.csv",header = FALSE)
Time_modulediscover_string = Time_modulediscover_string[seq(1,78,6),]
Time_modulediscover_iRefIndex = Time_modulediscover_iRefIndex[seq(1,78,6),]

# RFIM
Time_RFIM_string = read.csv("../../results/RFIM/String_time_original.txt",header = FALSE)
Time_RFIM_iRefIndex = Time_RFIM_string[seq(79,156,6),]
Time_RFIM_string = Time_RFIM_string[seq(1,78,6),]

# DOMINO
Time_DOMINO_string = read.csv("../DOMINO/time_string.txt",header = FALSE)
Time_DOMINO_iRefIndex = read.csv("../DOMINO/time_iRefIndex.txt",header = FALSE)
Time_DOMINO_string = Time_DOMINO_string[seq(1,78,6),]
Time_DOMINO_iRefIndex = Time_DOMINO_iRefIndex[seq(1,78,6),]

# robust
Time_robust_string = read.csv("../../results/robust/robust_time_string.csv",header = FALSE)
Time_robust_iRefIndex = Time_robust_string[seq(79,156,6),]
Time_robust_string = Time_robust_string[seq(1,78,6),]

# hierarchical-hotnet
Time_hotnet_string_1 = read.csv("../hierarchical-hotnet/time_string_1.txt",header = FALSE)
Time_hotnet_string_2 = read.csv("../hierarchical-hotnet/time_string_2.txt",header = FALSE)
Time_hotnet_string_3 = read.csv("../hierarchical-hotnet/time_string_3.txt",header = FALSE)
Time_hotnet_string_4 = read.csv("../hierarchical-hotnet/time_string_4.txt",header = FALSE)
Time_hotnet_string_5 = read.csv("../hierarchical-hotnet/time_string_5.txt",header = FALSE)
Time_hotnet_iRefIndex_1 = read.csv("../hierarchical-hotnet/time_iRefIndex_1.txt",header = FALSE)
Time_hotnet_iRefIndex_2 = read.csv("../hierarchical-hotnet/time_iRefIndex_2.txt",header = FALSE)
Time_hotnet_iRefIndex_3 = read.csv("../hierarchical-hotnet/time_iRefIndex_3.txt",header = FALSE)
Time_hotnet_iRefIndex_4 = read.csv("../hierarchical-hotnet/time_iRefIndex_4.txt",header = FALSE)
Time_hotnet_iRefIndex_5 = read.csv("../hierarchical-hotnet/time_iRefIndex_5.txt",header = FALSE)
Time_hotnet_string = Time_hotnet_string_1$V1+Time_hotnet_string_2$V1+Time_hotnet_string_3$V1+Time_hotnet_string_4$V1+Time_hotnet_string_5$V1
Time_hotnet_iRefIndex = Time_hotnet_iRefIndex_1$V1+Time_hotnet_iRefIndex_2$V1+Time_hotnet_iRefIndex_3$V1+Time_hotnet_iRefIndex_4$V1+Time_hotnet_iRefIndex_5$V1


data_time = data.frame(time_cost = c(Time_diamond_string,Time_diamond_iRefIndex,
                      Time_modulediscover_string,Time_modulediscover_iRefIndex,Time_RFIM_string,Time_RFIM_iRefIndex,
                      Time_DOMINO_string,Time_DOMINO_iRefIndex,Time_robust_string,Time_robust_iRefIndex,
                      Time_hotnet_string,Time_hotnet_iRefIndex),
             method=c(rep('DIAMOnD',26),rep('ModuleDiscoverer',26),rep('RFIM',26),
                      rep('DOMINO',26),rep('robust',26),rep('hierarchical-hotnet',26)),
             interactome = rep(c(rep('STRING',13),rep('iRefIndex',13)),6))

data_time_2 <- data_summary(data_time, varname="time_cost", groupnames=c("method",'interactome'))
data_time_2$method <- factor(data_time_2$method,levels = c("DIAMOnD","DOMINO","hierarchical-hotnet",'ModuleDiscoverer',"robust","RFIM"))

g1 = ggplot(data_time_2, aes(x = method, time_cost,group=interactome,fill=interactome)) + geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=time_cost-sd, ymax=time_cost+sd), width=.2,position=position_dodge(.9))+
  scale_fill_manual(values = c("#1b9e77","#d95f02"))+
  ylab('Computational time (in sec.)')+xlab('')+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+ 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5),
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        axis.ticks = element_line(colour = "black", size = 0.5))

ggsave(g1, file="../../figures/figure2.pdf", width=6, height=5)

