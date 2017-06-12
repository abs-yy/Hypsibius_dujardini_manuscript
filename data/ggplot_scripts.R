library(ggsci)
library(ggplot2)
library(reshape)
library(scales)

## Figure 2A
data<-read.delim("https://raw.githubusercontent.com/abs-yy/Hypsibius_dujardini_manuscript/master/data/Fig2A_HGT-content-in-metazoa.txt", header=T)
a<-melt(data, id.var=c("Category", "Organism"))
a<-a[-c(29,58,87),]
a$variable<-ifelse(a$variable=="Ab.initio", "Ab initio prediction", ifelse(a$variable=="This.work", "This work", "ENSEMBL gene set"))
colnames(a)<-c("Category", "Species", "Origin", "HGT_percentage")
a$Origin<-factor(a$Origin, levels=c("ENSEMBL gene set", "Ab initio prediction", "This work"))
ggplot(a, aes(x=Species, y=HGT_percentage, fill=Origin)) + geom_bar(position=position_dodge(), stat="identity", color="black") + facet_grid(~Category, scale="free", space="free") + theme_bw() + scale_y_continuous(breaks=seq(0,15,1)) + ylab("HGT percentage (%)") + xlab("Organism")  + theme(strip.background=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=12, face="italic"), axis.text.y=element_text(size=12), strip.text=element_text(size=12), legend.position=c(.5,.8), legend.title=element_text(size=16, face="bold"), legend.text=element_text(size=15)) + scale_fill_npg()
ggsave("Fig_2A.tiff", dpi=300)
#Saving 12.5 x 7.38 in image


## Figure 2B
data<-read.csv("https://raw.githubusercontent.com/abs-yy/Hypsibius_dujardini_manuscript/master/data/Fig2B_HGT_content_in_Hdujardini.csv", header=T)
a<-melt(data, id.var=c("Category", "Likelihood"))
a$variable<-ifelse(a$variable=="Expression_and_Conservation", "Expression and conservation", ifelse( a$variable=="Conservation_only", "Conservation only", ifelse(a$variable=="Expression_only", "Expression only", "Neither expression nor conservation" )))
colnames(a)<-c("assignment", "Likelihood", "Category", "Number")
a$assignment<-factor(a$assignment, levels=c("Prokaryote", "Non metazoan", "Virus", "Complex HGTs", "Metazoan", "Complex"))
ggplot(a, aes(x=assignment, y=Number, fill=Category)) + geom_bar(stat="identity") + xlab("Phylogenetic assignment of HGT candidates") + ylab("Number of HGT candidates") + facet_grid(~Likelihood, space="free", scale="free") + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10), strip.background=element_blank(), legend.position=c(.7,.8), legend.title=element_text(size=16, face="bold"), legend.text=element_text(size=15), strip.text=element_text(size=12)) + scale_fill_npg()
ggsave("Fig2B.tiff", dpi=300)
#Saving 9.27 x 5.72 in image

## Figure 3B
data<-read.delim("https://raw.githubusercontent.com/abs-yy/Hypsibius_dujardini_manuscript/master/data/Fig3B_gene-expression-DEGs-FCandTUN-annotated.txt", header=T)
data$Sample<-ifelse(data$Sample=="Hdujardini", "H. dujardini", ifelse(data$Sample=="Rvarieornatus.fast", "R. varieornatus (Fast)", "R. varieornatus (Slow)"))
data$Sample<-factor(data$Sample, levels=c("R. varieornatus (Fast)", "R. varieornatus (Slow)", "H. dujardini"))
data$Annotation<-factor(data$Annotation, levels=c("CAHS", "SAHS", "MAHS", "RvLEAM", "GST", "SOD", "CATE", "HSP", "Other"))
ggplot(data) + geom_point(data=data[data$Annotation=="Other",], aes(x=Tun, y=FC, shape=Conservation), color="gray")   +geom_point(data=data[data$Annotation!="Other",], aes(x=Tun, y=FC, color=Annotation, shape=Conservation )) + scale_x_continuous(trans="log1p", breaks=c(1,5,10,50,100,500,1000,5000,10000,50000), lim=c(1,NA), labels=comma) + scale_y_continuous(labels=comma, trans="log1p", breaks=c(1,5,10,50,100,500,1000,5000), lim=c(1, NA)) + facet_grid(~Sample) + theme_bw() + theme(strip.background=element_blank(), strip.text=element_text(size=10, face="italic"), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10), axis.text.y=element_text(size=10)) + xlab("Expression in Tun state") + ylab("Fold change ( Tun/Active+0.1)")
ggsave("Fig_3B.tiff", dpi=300)
#Saving 9.57 x 4.72 in image


