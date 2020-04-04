# =================== Graphing Script
# Plot results across forecasting techniques and parametrised/trained from various samples, procedure applied on s2 (delinquents-only)


# ============= Library setup

#general
library(Hmisc)
library(ETLUtils)
library(ffbase)
library(ff)
library(tidyr)
library(dplyr)
library(data.table)
options(scipen=999)

#for plots
library(ggplot2)
library(scales)
library(ggthemes)
library(extrafont)
library(RColorBrewer)


# ============= Initialization
exp_version <- "-excl_closed" #toggle between specific experiment set name or null string for the default cases
# load datasets and stack

#### RANDOM DEFAULTS : Exp-truncation
it.name <- paste0("v2_5j(ii)", exp_version)
unpack.ffdf(paste0("EL",it.name))
dat.main <- data.table(Setup="d_s12-Exp", dat.EL.par)
it.name <- paste0("v2_5k(ii)", exp_version)
unpack.ffdf(paste0("EL",it.name))
dat.main <- rbind(dat.main,data.table(Setup="e_s22-Exp", dat.EL.par))
it.name <- paste0("v2_5l(ii)", exp_version)
unpack.ffdf(paste0("EL",it.name))
dat.main <- rbind(dat.main,data.table(Setup="f_s32-Exp", dat.EL.par))

#### MARKOVIAN DEFAULTS
it.name <- paste0("v2_6a(ii)", exp_version)
unpack.ffdf(paste0("EL",it.name))
dat.main <- rbind(dat.main,data.table(Setup="g_s12-Markov", dat.EL.par))
it.name <- paste0("v2_6b(ii)", exp_version)
unpack.ffdf(paste0("EL",it.name))
dat.main <- rbind(dat.main,data.table(Setup="h_s22-Markov", dat.EL.par))
it.name <- paste0("v2_6c(ii)", exp_version)
unpack.ffdf(paste0("EL",it.name))
dat.main <- rbind(dat.main,data.table(Setup="i_s32-Markov", dat.EL.par))

# - feature engineering
# since losses are discounted at various points (based on (g,d)-configuration), simply summing performing and defaulting won't do.
# However, if we take (g1,0)-configuration from any random setup, then it refers to defaulting all 
# loans at the very start - effectively the defaulting balance is then the same as the sum of all loan principals.
TotalBalance <- subset(dat.main, Measure == "CD" & Setup == "h_s22-Markov")[order(Threshold),][1, Bal_Def]
dat.main[, Balance := TotalBalance]
dat.main[, LossRate := Loss / Balance]




# ============= Loss rate graph

# ----- Manicured plot : CD-only
# -- data prep
toplot <- subset(dat.main, Measure == "CD")
setDT(toplot, key=c("Setup", "Threshold"))
plot.data2 <- toplot #for plotting points
Setups <- sort(unique(toplot$Setup))

#decrease overall number of points
treat.data <- subset(plot.data2, Threshold > 0) 
plot.data2 <- subset(plot.data2, !(Threshold > 0))
sel.rows <- seq(1,nrow(treat.data),by=4)
plot.data2 <- rbind(plot.data2, treat.data[sel.rows,])

# find global minima
plot.data3 <- toplot[order(Setup),list(LossRate=min(LossRate), ThresPosition = match(min(Loss),Loss) ),
                     by=list(Setup)]
plot.data3[, Threshold := sapply(1:.N, function(i,j) { 
  return( subset(toplot, Setup==Setups[i])[j[i], Threshold]) 
}, j=ThresPosition)]

#select only certain minima for main plot (others to be shown in zoomed plot)
plot.data4 <- subset(plot.data3, Setup %in% c("g_s12-Markov","h_s22-Markov"))
chosenFont <- "Times New Roman"

# create Legend text, with found minima in brackets
label.vec <- list(); i<- 1;
label.vec[[i]] <- parse(text=paste0("Random~(italic(k)%~%Exp)-italic(s)[12]~~(italic(d)^ '*' ==", plot.data3[i,Threshold], ")")); i <- i + 1
label.vec[[i]] <- parse(text=paste0("Random~(italic(k)%~%Exp)-italic(s)[22]~~(italic(d)^ '*' ==", plot.data3[i,Threshold], ")")); i <- i + 1
label.vec[[i]] <- parse(text=paste0("Random~(italic(k)%~%Exp)-italic(s)[32]~~(italic(d)^ '*' ==", plot.data3[i,Threshold], ")")); i <- i + 1
label.vec[[i]] <- parse(text=paste0("Markov-italic(s)[12]~~(italic(d)^ '*' ==", plot.data3[i,Threshold], ")")); i <- i + 1
label.vec[[i]] <- parse(text=paste0("Markov-italic(s)[22]~~(italic(d)^ '*' ==", plot.data3[i,Threshold], ")")); i <- i + 1
label.vec[[i]] <- parse(text=paste0("Markov-italic(s)[32]~~(italic(d)^ '*' ==", plot.data3[i,Threshold], ")")); i <- i + 1
size.v <- c(0.6, 0.9, 0.6, 0.6, 0.9, 0.6)
ltype.v <- c("dotted", "solid", "dotted", "dotted", "solid", "dotted")
shape.v <- c(1,0,2, 16,17,18) # 3,4,8,
col.v <- (brewer.pal(n=8, name="Set1"))[c(1,3,4,2,5,7)]

# -- main plot
min.y <- 0.15
max.y <- 0.48
plot.full <- ggplot(toplot) + theme_minimal() + 
  theme(text=element_text(size=12, family=chosenFont), legend.position="bottom") + 
  labs(y="Loss rate (%)", x = bquote({Thresholds~italic(d)~on~italic(g)[1]})) + 
  geom_line(aes(x=Threshold, y=LossRate, colour=Setup, size=Setup, linetype=Setup)) + 
  geom_point(aes(x=Threshold, y=LossRate, colour=Setup, shape=Setup), size=2.5, alpha=0.8, data=plot.data2) + 
  geom_segment(aes(x=0, xend=Threshold, y=LossRate, yend=LossRate, color=Setup), linetype="dashed", 
               data=plot.data4) + 
  geom_segment(aes(x=Threshold,xend=Threshold, y=min.y, yend=LossRate, color=Setup), linetype="dashed", 
               data=plot.data4) +  
  geom_point(aes(x=Threshold,y=LossRate), data=plot.data4, size=8, colour="gray15", shape=1) + 
  geom_text(aes(x=Threshold,y=LossRate,label=Threshold), data=plot.data4, size=2.5, colour="gray15") + 
  # scale options
  scale_color_manual(name = "Setup", labels=label.vec, values=col.v) + 
  scale_size_manual(name="Setup", labels=label.vec, values=size.v) + 
  scale_linetype_manual(name="Setup", labels=label.vec, values=ltype.v) + 
  scale_shape_manual(values=shape.v, name="Setup", labels=label.vec) + 
  guides(col=guide_legend(nrow=3)) + coord_cartesian(xlim=c(0,168),ylim=c(min.y,max.y)) +
  scale_y_continuous(limits=c(min.y,max.y), breaks=pretty_breaks(), 
                     labels=percent) + 
  scale_x_continuous(breaks=pretty_breaks())
plot.full


# zoom bounding box
xlim <- c(6,16); ylim <- c(0.235, 0.25)

# zoomed plot
plot.zoom <- 
  plot.full + coord_cartesian(xlim=xlim, ylim=ylim) + 
  geom_point(aes(x=Threshold,y=LossRate, color=Setup, shape=Setup), data=toplot, size=2.5) +
  geom_point(aes(x=Threshold,y=LossRate), data=plot.data3, size=8, colour="gray15", shape=1) + 
  geom_segment(aes(x=0,xend=Threshold, y=LossRate, yend=LossRate, color=Setup), linetype="dashed", 
               data=plot.data3) + 
  geom_segment(aes(x=Threshold,xend=Threshold, y=min.y, yend=LossRate, color=Setup), linetype="dashed", 
               data=plot.data3) +
  theme(legend.position="none", axis.text.x = element_text(margin=unit(c(0,0,0,0),"mm"),size=9),
        axis.text.y=element_text(margin=unit(c(0,0,0,0),"mm"),size=9),axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(color='black', fill="white"), 
        plot.background=element_rect(color="white"),
        plot.margin = unit(c(0,0,0,0),"mm"))
plot.zoom

# merge plots
use.ymin <-  0.3
use.ymax <- max.y

plot.merged <- plot.full + annotation_custom(grob = ggplotGrob(plot.zoom), xmin = 10, xmax=100,
                                             ymin = use.ymin, ymax = use.ymax)
plot.merged

dpi <- 200
# note the 9-version graph was previously set to LossThreshv2
ggsave(plot.merged, file=paste0("LossThreshv3_CD-On_s2_rate",exp_version,".png"),width=1200/dpi, height=1100/dpi,dpi=dpi)



