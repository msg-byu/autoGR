aflow <- read.delim("../data/R_workspace/MP_irkpd.csv",header=F,sep="\t")
gr_server <- read.delim("../data/R_workspace/GR_server_irkpd.csv",header=F,sep="\t")
gr_auto <- read.delim("../data/R_workspace/GR_auto_irkpd.csv",header=F,sep="\t")
library("scales")
mycol=c("green","red","blue")
opa.a <- 0.2
span.sm <- 0.5

x<-aflow
x<-log10(x)
ff1<-loess.smooth(x[,1],x[,2],span=span.sm,evaluation=200)

x<-gr_server
x<-log10(x)
ff3<-loess.smooth(x[,1],x[,2],span=span.sm,evaluation=200)

x<-gr_auto
x<-log10(x)
ff4<-loess.smooth(x[,1],x[,2],span=span.sm,evaluation=200)

write.csv(ff1,file="../data/R_workspace/MP_irkpd_loess.csv",row.names=FALSE)
write.csv(ff3,file="../data/R_workspace/GR_server_irkpd_loess.csv",row.names=FALSE)
write.csv(ff4,file="../data/R_workspace/GR_auto_irkpd_loess.csv",row.names=FALSE)