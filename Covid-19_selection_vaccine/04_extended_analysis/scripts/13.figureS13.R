rm(list = ls())
root_path <- "C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine"
setwd(root_path)
source("load_package.R")
source("plot_theme.R")

outpath0 <- "03_output/figure_2_3"
outpath <-  paste0(outpath0,"/","/muts_nonsyn_vaccine_202003_202110")
outpath1 <-  paste0(outpath,"/","vero_result",sep="")
outpath2 <-  paste0(outpath0,"/","muts_nonsyn_NI_202003_202110")
outpath3 <-   paste0(outpath2,"/","NI_result",sep="")

outpath4 <-  "04_extended_analysis/plot_output/Figure_S13/"
dir.create(outpath4)

#Plot
setwd(paste0(root_path,"/",outpath1,sep=""))
#rm(list = ls())
datapath <- paste0(root_path,"/",outpath1,sep="")
datafile <- list.files(datapath, pattern = "*.csv$", full.names = TRUE)
files=list()
read.txt <- function(x) {
  des <- readLines(x)                   
  return(paste(des, collapse = ""))     
}
review <- lapply(datafile, read.csv)
docname <- list.files(datapath, pattern = "*.csv$")
myfiles = do.call(rbind, lapply(docname, function(x) read.csv(x, stringsAsFactors = FALSE)))
D1_vero <- myfiles




D1_vero$gene <-factor(D1_vero$gene,levels= c('S','S1','S2','E','M','N','ORF1a','ORF1b','ORF3a','ORF6','ORF7a','ORF7b','ORF8','ORF9b',
                                         ordered=TRUE))

D1_vero$coverage <-factor(D1_vero$coverage,levels= c('Low coverage (25)','Medium coverage (50)','High coverage (75)',
                                                 ordered=TRUE))
D1_vero <- na.omit(D1_vero)


p1  <- ggplot(D1_vero,aes(x=gene, y=vero_fit,color=coverage,group=coverage))+ 
  geom_point(aes(x=gene,y=vero_fit,fill=coverage),position=position_dodge(0.75),size=7,alpha=1,shape=21,color="black")+
  geom_errorbar(aes(x=gene, y=vero_fit,ymin=vero_low,ymax=vero_high),
                width=0.2,stat="identity",position = position_dodge(0.75),lwd=0.6,color="black")+
  scale_y_continuous(limits = c(-24,24),breaks = seq(-24,24,12))+
  labs(x="", y= "Effect on number of nonsynonymous\n mutations" )+  
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",linewidth=0.5)+
  scale_color_manual(name = 'Adjusted vaccine coverage',values = c("#92C5DE","#4393C3","#2166AC"))+
  scale_fill_manual(name = 'Adjusted vaccine coverage',values = c("#92C5DE","#4393C3","#2166AC"))+ 
  mytheme+
  ggtitle("Adjusted vaccine coverage")+
  theme(strip.text = element_text(face="bold",size=18))+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  #theme(legend.position = "none")+
  guides(color=guide_legend(nrow=1))
#ggtitle("Correlation between Ka/Ks and international travel")

p1


setwd(paste0(root_path,"/",outpath3,sep=""))
datapath <- paste0(root_path,"/",outpath3,sep="")
datafile <- list.files(datapath, pattern = "*.csv$", full.names = TRUE)
files=list()
read.txt <- function(x) {
  des <- readLines(x)                   
  return(paste(des, collapse = ""))     
}
review <- lapply(datafile, read.csv)
docname <- list.files(datapath, pattern = "*.csv$")
myfiles = do.call(rbind, lapply(docname, function(x) read.csv(x, stringsAsFactors = FALSE)))
D1_Ni <- myfiles

D1_Ni$gene <-factor(D1_Ni$gene,levels= c('S','S1','S2','E','M','N','ORF1a','ORF1b','ORF3a','ORF6','ORF7a','ORF7b','ORF8','ORF9b',
                                         ordered=TRUE))

D1_Ni$coverage <-factor(D1_Ni$coverage,levels= c('Low coverage (25)','Medium coverage (50)','High coverage (75)',
                                                 ordered=TRUE))

D1_Ni <- na.omit(D1_Ni)


p2  <- ggplot(D1_Ni,aes(x=gene, y=NI_fit,color=coverage,group=coverage))+ 
  geom_point(aes(x=gene,y=NI_fit,fill=coverage),position=position_dodge(0.5),size=7,shape=21,color="black")+
  geom_errorbar(aes(x=gene, y=NI_fit,ymin=NI_low,ymax=NI_high),
                width=0.2,stat="identity",position = position_dodge(0.5),lwd=0.6,alpha=1,color="black")+
  scale_y_continuous(limits = c(-8,24),breaks = seq(-8,24,8))+
  labs(x="", y= "Effect on number of nonsynonymous\n mutations" )+  
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",linewidth=0.5)+
  scale_color_manual(name = 'Natural immunity coverage',values = c("#F4A582","#D6604D","#B2182B"))+ 
  scale_fill_manual(name = 'Natural immunity coverage',values = c("#F4A582","#D6604D","#B2182B"))+ 
  
  mytheme+
  ggtitle("Natural immunity")+
  theme(strip.text = element_text(face="bold",size=18))+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  guides(color=guide_legend(nrow=1))

#ggtitle("Correlation between Ka/Ks and international travel")

p2



p<-ggarrange(p1, p2,ncol =1, nrow = 2,widths = c(1,1),heights = c(1,1),
             labels = c("A","B"),  vjust=1.1,align = "hv",
             #legend = "bottom",
             font.label = list(size = 30, face = "bold"))
p

ggsave(p,filename = file.path(root_path,"/",outpath4,"/","figure_S13.pdf"),
       width = 15.5,height =13)


##Perform a significance test 
##Determining if there is a significant difference in associations between adjusted vaccine coverage and selection pressure on S and non-S proteins.
setwd(paste0(root_path,"/",outpath,sep=""))
datapath <- paste0(root_path,"/",outpath,sep="")
datafile <- list.files(datapath, pattern = "*.csv$", full.names = TRUE)
files=list()
read.txt <- function(x) {
  des <- readLines(x)                   
  return(paste(des, collapse = ""))     
}
review <- lapply(datafile, read.csv)
docname <- list.files(datapath, pattern = "*.csv$")
myfiles = do.call(rbind, lapply(docname, function(x) read.csv(x, stringsAsFactors = FALSE)))
D1_vero_all <- myfiles

D1_vero_all$gene <-factor(D1_vero_all$gene,levels= c('S','S1','S2','E','M','N','ORF1a','ORF1b','ORF3a',
                                               'ORF6','ORF7a','ORF7b','ORF8','ORF9b',ordered=TRUE))
D1_vero_all <- na.omit(D1_vero_all)

S_ind1 <- which(D1_vero_all$gene=="S"|D1_vero_all$gene=="S1"|D1_vero_all$gene=="S2")

S_vero <- D1_vero_all[S_ind1,]
non_S_vero <- D1_vero_all[-S_ind1,]

t.test(S_vero$dat_vero_fit , non_S_vero$dat_vero_fit, var.equal = T)


##Perform a significance test 
##Determining if there is a significant difference in associations between natural immunity and selection pressure on S and non-S proteins.
setwd(paste0(root_path,"/",outpath2,sep=""))
datapath <- paste0(root_path,"/",outpath2,sep="")
datafile <- list.files(datapath, pattern = "*.csv$", full.names = TRUE)
files=list()
read.txt <- function(x) {
  des <- readLines(x)                   
  return(paste(des, collapse = ""))     
}
review <- lapply(datafile, read.csv)
docname <- list.files(datapath, pattern = "*.csv$")
myfiles = do.call(rbind, lapply(docname, function(x) read.csv(x, stringsAsFactors = FALSE)))
D1_Ni_all <- myfiles
D1_Ni_all$gene <-factor(D1_Ni_all$gene,levels= c('S','S1','S2','E','M','N','ORF1a','ORF1b','ORF3a',
                                                     'ORF6','ORF7a','ORF7b','ORF8','ORF9b',ordered=TRUE))
D1_Ni_all <- na.omit(D1_Ni_all)

S_ind2 <- which(D1_Ni_all$gene=="S"|D1_Ni_all$gene=="S1"|D1_Ni_all$gene=="S2")
S_Ni <- D1_Ni_all[S_ind2,]
non_S_Ni <- D1_Ni_all[-S_ind2,]

t.test(S_Ni$dat_NI_fit, non_S_Ni$dat_NI_fit, var.equal = T)


