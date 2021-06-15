
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#009E73")
cbPalette1 <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
cbPalette2 <- c("#a2b9bc","#b2ad7f","#878f99","#6b5b95","#feb236","#d64161","#ff7b25","#d6cbd3","#eca1a6","#bdcebe","#ada397")

library(reshape2)
library(ggplot2)
library(scales)
library(ggbeeswarm)


## heterozygosity comparison

ht = read.table("heterCombined.txt1")
head(ht)
unique(ht$V1)

plot = ggplot(data = ht, aes(x=V1, y=round(V4,4), fill = V6, shape=V5))+ geom_beeswarm(cex=1.5,size=2.5,priority = "density") + scale_shape_manual(values = c(21,22))
plot = plot + scale_y_continuous(limits = c(0,0.01),name = "heterozygosity", expand = c(0.0001,0.0001), labels = percent) + labs(fill="",shape="")
plot = plot + scale_fill_manual(values = c(cbPalette[2], cbPalette[1]))+ scale_x_discrete(name = "", limits=c("animals", "mammals", "ruminants","rhinoceros"))
plot = plot + theme(panel.background = element_blank(),
             panel.grid = element_blank(),
             axis.line.x = element_line(),
             axis.line.y = element_line(),
             legend.position = "top",
             axis.text = element_text(size = 12))



## psmc

minBreak= function(x,y){
  t=c()
  for (i in x:y){
    for (m in 1:9){
      n=m*10^i
      t=c(t,n)
    }
  }
  return(t)
}

minLabels = function(x,y){
  t=c()
  for (i in x:y){
    t=c(t,paste0("1e+",i))
    for (m in 2:9){
      t=c(t,"")
    }
  }
  return(t)
}


scaleFUN <- function(x) sprintf("%.1f", x)

data = read.table("psmc.difParaAdjusted")
head(data)

data_a = subset(data, V2=="one")

plota = ggplot(data=data_a, aes(V4,V5))+geom_step(aes(color=V1),size=1.5) +scale_color_manual(values = cbPalette) + labs(color = "species")
plota = plota + scale_x_log10(name = "g=12, u=1.2e-8",expand = c(0,0), limits=c(10000,10000000),breaks=minBreak(4,6), labels = minLabels(4,6), sec.axis=dup_axis())
plota = plota +scale_y_continuous(name = "Effective population size (10^4)", limits=c(0,9.5), expand = c(0,0), sec.axis = dup_axis())
plota = plota  + theme(panel.background = element_blank(),
             panel.grid = element_blank(), 
             axis.line = element_line(), 
             axis.text.x = element_blank(), 
             axis.text.y.right = element_blank(), 
             axis.title.x.top = element_blank(),
             axis.title.y.right = element_blank(),
             legend.position = "none"
             )

plota

data_b = subset(data, V2=="cite")

plotb = ggplot(data=data_b, aes(V4,V5))+geom_step(aes(color=V1),size=1.5) +scale_color_manual(values = cbPalette) + labs(color = "species")
plotb = plotb + scale_x_log10(name = "g=12, u=2.34e-8",expand = c(0,0), limits=c(10000,10000000),breaks=minBreak(4,6), labels = minLabels(4,6), sec.axis=dup_axis())
plotb = plotb +scale_y_continuous(name = "Effective population size (10^4)", limits=c(0,7), labels = scaleFUN, expand = c(0,0), sec.axis = dup_axis())
plotb = plotb  + theme(panel.background = element_blank(),
                       panel.grid = element_blank(), 
                       axis.line = element_line(), 
                       axis.text.x = element_blank(), 
                       axis.text.y.right = element_blank(), 
                       axis.title.x.top = element_blank(),
                       axis.title.y.right = element_blank(),
                       legend.position = "none"
)

plotb

data_c = subset(data, V2=="thr")

plotc = ggplot(data=data_c, aes(V4,V5))+geom_step(aes(color=V1),size=1.5) +scale_color_manual(values = cbPalette) + labs(color = "species")
plotc = plotc + scale_x_log10(name = "g=12, u=3.3e-8",expand = c(0,0), limits=c(10000,10000000),breaks=minBreak(4,6), labels = minLabels(4,6), sec.axis=dup_axis())
plotc = plotc +scale_y_continuous(name = "Effective population size (10^4)", limits=c(0,4), labels = scaleFUN, expand = c(0,0), sec.axis = dup_axis())
plotc = plotc  + theme(panel.background = element_blank(),
                       panel.grid = element_blank(), 
                       axis.line = element_line(), 
                       axis.text.x.top  = element_blank(), 
                       axis.text.y.right = element_blank(), 
                       axis.title.x.top = element_blank(),
                       axis.title.y.right = element_blank(),
                       legend.position = "none"
)

plotc


psmc_e = read.table("psmc.palaearcticAdjusted")
psmc_e = subset(psmc_e, V2==0)
head(psmc_e)

pe_psmc=ggplot(data=psmc_e, aes(V3,V4))+geom_step(aes(color=V1),size=1.5)+scale_color_manual(values = cbPalette)
pe_psmc=pe_psmc+scale_x_log10(name = "g=10, u=1.95e-8", expand = c(0,0), limits=c(10000,10000000),breaks=minBreak(4,6), labels = minLabels(4,6),sec.axis=dup_axis())+scale_y_continuous(name = expression(Effective~population~size~(10^{4})),labels = scaleFUN, expand = c(0,0), sec.axis = dup_axis(),limits = c(0,7))
pe_psmc=pe_psmc+theme(
  panel.background = element_blank(), 
  panel.grid = element_blank(), 
  axis.line = element_line(),
  axis.text.x = element_blank(),
  axis.text.y.right = element_blank(),
  axis.title.x.top = element_blank(),
  axis.title.y.right = element_blank(),
  legend.position = "none"
)

pe_psmc

gridExtra::grid.arrange(plota,plotb,pe_psmc, plotc, nrow=4, heights = c(24,24,24,28))


## snpEff

rd = read.csv("mammalMean.txt5", header = F)

plot = ggplot(data = rd, aes(x=V5, y=V6)) +geom_point(aes(color=V1, shape=V1),size=4)+geom_text_repel(aes(label=V2),size=4)
plot = plot + scale_x_continuous(name = "mean rate of Missense / Slient") + scale_y_continuous(name = "mean value of LoF mutation rate")
plot = plot +  geom_smooth(method = "lm", formula = y~x, color="black",size=1, se=F) + scale_shape_manual(values = c(19,15))+scale_color_manual(values = c("#D55E00","#999999"))
plot = plot  + annotate("text", x=0.6, y=0.015, label = "atop(italic(R) ^ 2 == 0.61, 'P value = 9.43e-5')", parse=T, size=6)
plot +  theme(panel.background = element_blank(),
             panel.grid = element_blank(),
             axis.line  = element_line(),
             axis.text = element_text(size = 12),
             axis.title = element_text(size = 12),
             legend.position = "none"
)

## roh

roh=read.table("roh.sum",header = F)
roh$V8=rep("MB",48)
roh$V8 = paste0(roh$V4/1000,roh$V8)

rohTmp=roh[with(roh,V2==500),]
rohTmp=rohTmp[with(rohTmp,rev(order(V3))),]$V1

roh$V1=factor(roh$V1,levels = rev(rohTmp))
roh$V1=factor(roh$V1,levels = c("sumatranRhino", "indianRhino", "blackRhino","whiteRhino"))

rohTmp=roh[with(roh,V2==500),]$V1

rohN=subset(roh,V2=="nondenovo")
rohpN=ggplot(data = rohN, aes(V3, V7/1000))+geom_col(aes(fill=V8), width=0.7, position = position_stack(reverse = T))+coord_flip(expand = T)
rohpN=rohpN+scale_x_discrete(name="")+scale_y_continuous(name = "sum of ROH length (MB)",sec.axis = dup_axis(),limits = c(0,800))
rohpN=rohpN+theme(axis.line.x = element_line(),panel.background = element_blank())+labs(fill="ROH categories")+scale_fill_manual(values = cbPalette)

rohD=subset(roh,V2=="denovo" & V1=="all_sites")
rohpD=ggplot(data = rohD, aes(V1, V2/1000))+geom_col(aes(fill=V8), width=0.7, position = position_stack(reverse = T))+coord_flip(expand = T)
rohpD=rohpD+scale_x_discrete(name="")+scale_y_continuous(name = "sum of ROH length (MB)",sec.axis = dup_axis(),limits = c(0,800))
rohpD=rohpD+theme(axis.line.x = element_line(),panel.background = element_blank())+labs(fill="ROH categories")+scale_fill_manual(values = cbPalette)

gridExtra::grid.arrange(rohpD,rohpN,nrow=2)



