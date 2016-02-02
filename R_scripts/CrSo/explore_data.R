setwd("~/Documents/Scripts/R/CrSo2015")
library(plotrix)
library(ggplot2)
Nmax=read.table("Nmax.txt",header=F,sep="\t")
rNmax=read.table("rNmax.txt",header=F,sep="\t")
averNmax=read.table("ave_rNmax.txt",header=F,sep="\t")

## histograms
#hist(rNmax$V2,col=2,breaks=10,density=10)
#hist(averNmax$V2,col=2,breaks=3,density=10)
#hist(Nmax$V2,col=2,breaks=10)
#x=list(averNmax$V2,Nmax$V2)
#multhist(x,breaks=10)

names=c(rep('real Genes',times=length(Nmax$V2)),rep('mean simulated Genes',times=length(averNmax$V2)))
values = c(Nmax$V2,averNmax$V2)
df = data.frame(names,values)
names(df) <- c("Sample","Value")
attach(df)
head (df)
ggplot(df, aes(x=Value,group=Sample,fill=Sample))+
  ylim(0,0.01)+xlim(0,300)+ 
  geom_density(position="identity", alpha=0.5,right=TRUE) +
  theme_bw() +
  geom_vline(df,aes(xintercept=mean(Value, na.rm=T)),linetype="dashed",size=1)+
  xlab("Nmax values")


avevalues = c(Nmax$V3,averNmax$V3,rNmax$V3)
avedf = data.frame(names,avevalues)
names(avedf) <- c("Sample","aveValue")
attach(avedf)
ggplot(avedf, aes(x=aveValue,group=Sample,fill=Sample))+
  ylim(0,0.15)+xlim(0,20)+ 
  geom_density(position="identity", alpha=0.5,right=TRUE) +
  theme_bw() +
  geom_vline(linetype="dashed",size=1)+
  xlab("average Nmax values")

summary(Nmax$V2)
summary(Nmax$V3)
summary(rNmax$V2)

## normality tests
shapiro.test(averNmax$V2)
ks.test(rNmax$V2,"pnorm",mean(rNmax$V2),sqrt(var(rNmax$V2)))
qqnorm(averNmax$V2)
qqline(averNmax$V2)


##@@ SS data @@##

ss=read.table("test_ss.txt",header=F,sep="\t")
rss=read.table("test_rss.txt",header=F,sep="\t")

# ss=read.table("ss.txt",header=F,sep="\t")
# rss=read.table("./random/ave_ss.txt",header=F,sep="\t")

ssnames=c(rep('real Genes',times=length(ss$V2)),rep('simulated Genes',times=length(rss$V2)))
ssvalues = c(ss$V2,rss$V2)
dfss = data.frame(ssnames,ssvalues)
names(dfss) <- c("Sample","Value")
attach(dfss)
ggplot(dfss, aes(x=Value,group=Sample,fill=Sample))+
  ylim(0,2)+xlim(0,15)+ 
  geom_density(position="identity", alpha=0.5,right=TRUE) +
  theme_bw() +
 # geom_vline(linetype="dashed",size=1)+
  xlab("SS values")
  


