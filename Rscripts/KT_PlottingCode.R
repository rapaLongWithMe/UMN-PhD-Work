library(dplyr)
library(ggplot2)
library(Rmisc)
library(reshape2)

#Set working directory and read in timecourse normalized log2 transformed fpkm data
setwd("/Users/greenham/Documents/Brassica_Genomics/JGI_Morphotypes/RNAseq/")

TCfilter=read.csv("R500_filtered.csv")

#Get dimensions to fill in aggregate length
dim(TCfilter)
TCmean=aggregate(TCfilter[,3:20861],list(TCfilter$Time),mean)
TCsd=aggregate(TCfilter[,3:20861],list(TCfilter$Time),sd)

colnames(TCmean)[1]<-"Time"
colnames(TCsd)[1]<-"Time"

#Transpose to make it easier to pull out genes.
TCmeanT<-as.data.frame(t(TCmean))
TCsdT<-as.data.frame(t(TCsd))

#Make the 'Time' row the column names
colnames(TCmeanT)<-unlist(TCmeanT["Time",])
TCmeanT<-TCmeanT[!row.names(TCmeanT)=='Time',]

colnames(TCsdT)<-unlist(TCsdT["Time",])
TCsdT<-TCsdT[!row.names(TCsdT)=='Time',]

#Remove time point 8 for plotting
drop<-c("CTP8","WTP8")
TCmeanT<-TCmeanT[,!names(TCmeanT)%in%drop]
TCsdT<-TCsdT[,!names(TCsdT)%in%drop]

#Convert rownames to GeneID column
GeneID<-rownames(TCmeanT)
rownames(TCmeanT)<-NULL
TCmeanT<-cbind(GeneID,TCmeanT)

GeneID<-rownames(TCsdT)
rownames(TCsdT)<-NULL
TCsdT<-cbind(GeneID,TCsdT)

#Subset the gene to plot.
GeneM<-TCmeanT[TCmeanT$GeneID=='BraA03g51570R',]
GeneSD<-TCsdT[TCsdT$GeneID=='BraA03g51570R',]

#Convert to long format
GeneMlong<-melt(GeneM,id.vars=c("GeneID"),variable.name="Time",value.name="Mean")
GeneSDlong<-melt(GeneSD,id.vars=c("GeneID"),variable.name="Time",value.name="SD")
GeneMSD<-join(GeneMlong,GeneSDlong)

#Adding new Treatment and Time columns
Treatment<-c("Cold","Cold","Cold","Cold","Cold","Cold","Cold","Control","Control","Control","Control","Control","Control","Control")
Time<-c("1","2","3","4","5","6","7","1","2","3","4","5","6","7")
GeneMSD<-select(GeneMSD,-c(Time))
GeneMSDall<-cbind(GeneMSD,Time,Treatment)
GeneMSDall$Mean<-as.numeric(as.character(GeneMSDall$Mean))
GeneMSDall$SD<-as.numeric(as.character(GeneMSDall$SD))

#Insert desired title.
plot<-ggplot(GeneMSDall,aes(x=Time,y=Mean,group=Treatment,colour=Treatment)) + geom_line(size=1) + ggtitle("R500_CBF2") + xlab("Time Point") + ylab("log2(FPKM)") + geom_ribbon(aes(ymax=Mean+SD,ymin=Mean-SD,fill=Treatment),linetype=0,alpha=1/5) + scale_fill_manual(values=c("#3399FF","#CC0066")) + scale_color_manual(values=c("#3399FF","#CC0066")) + theme_bw() + theme(legend.title=element_text(size=14),legend.text=element_text(size=12),panel.grid.major=element_blank(), panel.grid.minor=element_blank(),panel.background=element_blank(),panel.border = element_rect(colour="black",size=2), axis.text.y=element_text(size=12), axis.title.y=element_text(size=14),axis.text.x=element_text(size=12,angle=45,hjust=1),axis.title.x=element_text(size=14))
pdf(file="R500_ColdTC_CBF2.pdf",width=8,height=6)
plot
dev.off()

#Plot control only
control<-subset(GeneMSDall,Treatment=="Control")
plot<-ggplot(control,aes(x=Time,y=Mean,group=Treatment,colour=Treatment)) + 
  geom_line(size=1) + ggtitle("HN53_PHYB") + 
  xlab("Time Point") + ylab("log2(FPKM)") + 
  geom_ribbon(aes(ymax=Mean+SD,ymin=Mean-SD,fill=Treatment),linetype=0,alpha=1/5) + 
  scale_fill_manual(values=c("#3399FF")) + 
  scale_color_manual(values=c("#3399FF")) + 
  theme_bw() + 
  theme(legend.title=element_text(size=14),legend.text=element_text(size=12),panel.grid.major=element_blank(), panel.grid.minor=element_blank(),panel.background=element_blank(),panel.border = element_rect(colour="black",size=2), axis.text.y=element_text(size=12), axis.title.y=element_text(size=14),axis.text.x=element_text(size=12,angle=45,hjust=1),axis.title.x=element_text(size=14))


pdf(file="HN53_ControlTC_PHYB.pdf",width=8,height=6)
plot
dev.off()

