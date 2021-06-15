
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#009E73")
cbPalette1 <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

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



