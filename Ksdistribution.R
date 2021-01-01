rm(list=ls())
library(getopt)
command <- matrix(
   c("input","i",1,"character","Need Ks distribution file",
     "output","o",2,"character", "Option: outfile name, default name is Ks.distribution.pdf",
     "xlim","x",1,"integer","xlim upper limit value",
     "help","h",0,"logical","help information"),
     byrow=T,ncol=5)
args <- getopt(spec = command)
if (is.null(args$input)) {
  cat(paste(getopt(spec=command, usage = T), "\n"))
  quit()
}
if (is.null(args$xlim)){
   args$xlim = 3
}
#if(is.null(args$help)){
#  print -> print.default
#}else{
#  cat(paste(getopt(spec=command, usage = T), "\n"))
#  quit()
#}
library(ggplot2)
data<-read.table(args$input,header=TRUE)
p <- ggplot(data,aes(Ks,fill=Species,color=Species,alpha=0.8))  +
  geom_density()  + xlim(0,args$xlim)  + theme_classic()  + 
  guides(alpha=FALSE)  +
  theme(axis.title = element_text(size=16),axis.text=element_text(size=16),legend.position =  "top")
if(is.null(args$output)){
  ggsave(p,filename = "Ks.distribution.pdf",width = 12,height = 9)
}else{
  ggsave(p,filename = args$output, width = 12,height = 9)
}
