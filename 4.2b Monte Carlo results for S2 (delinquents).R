# =================== Graphing Script
# Plot Monte Carlo results across forecasting techniques
# - Forecasts trained from: S2
# - Optimisation applied on: S2


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
# number of MC trials: for naming purposes
MC.size <- 500
# load datasets and stack
it.name <- paste0("MonteCarlo_EL", MC.size, "-v1_2a")
unpack.ffdf(paste0(it.name))
dat.main <- data.table(Setup="a_s22-Exp", dat.EL.MC)
it.name <- paste0("MonteCarlo_EL", MC.size,"-v1_2b")
unpack.ffdf(paste0(it.name))
dat.main <- rbind(dat.main,data.table(Setup="b_s22-Markov", dat.EL.MC))


# - feature engineering
# since losses are discounted at various points (based on (g,d)-configuration), simply summing performing and defaulting won't do.
# However, if we take (g1,0)-configuration from any random setup, then it refers to defaulting all 
# loans at the very start - effectively the defaulting balance is then the same as the sum of all loan principals.
TotalBalance <- subset(dat.main, Measure == "CD" & Setup == "a_s22-Exp" & Iteration == 1 & Threshold == 0)[1, Bal_Def]
dat.main[, Balance := TotalBalance]
dat.main[, LossRate := Loss / Balance]




# ============= Average Loss rate graph

# ----- Manicured plot : CD-only

# -- data prep
toplot <- subset(dat.main, Measure == "CD")
setDT(toplot, key=c("Setup", "Iteration", "Threshold"))

# sampling distribution of the loss rate, respective to a given threshold
# from https://statweb.stanford.edu/~owen/mc/Ch-intro.pdf
toplot[, Mean_LossRate := mean(LossRate, na.rm=T), by=list(Setup, Threshold)]
toplot[, SampleVar_LossRate := var(LossRate, na.rm=T), by=list(Setup, Threshold)] # s^2
toplot[, SampleSize_LossRate := max(Iteration, na.rm=T), by=list(Setup, Threshold)] # n
toplot[, Var_Mean_LossRate := SampleVar_LossRate / SampleSize_LossRate] # s^2 / n
toplot[, std_Mean_LossRate := sqrt(Var_Mean_LossRate)]
# 95% confidence intervals for the means of the sample distributions of the loss rate
toplot[, Mean_LossRate_Upper := Mean_LossRate + 2.58*std_Mean_LossRate] # 2.58 for 99% CI; 1.96 for 95% CI
toplot[, Mean_LossRate_Lower := Mean_LossRate - 2.58*std_Mean_LossRate]

toplot.aggr <- unique(toplot[,list(Setup, Threshold, Mean_LossRate, std_Mean_LossRate, Mean_LossRate_Upper, Mean_LossRate_Lower)])
toplot.aggr[, LossRate := Mean_LossRate] # just to satisfy ggplot2's moaning

# -- find global minima
Setups <- sort(unique(toplot.aggr$Setup))
plot.data3 <- toplot.aggr[order(Setup),list(LossRate=min(Mean_LossRate), ThresPosition = match(min(Mean_LossRate),Mean_LossRate) ),
                          by=list(Setup)]
plot.data3[, Threshold := sapply(1:.N, function(i,j) { 
  return( subset(toplot.aggr, Setup==Setups[i])[j[i], Threshold]) 
}, j=ThresPosition)]

# -- create Legend text, with found minima in brackets
label.vec <- list(); i <- 1
label.vec[[i]] <- parse(text=paste0(MC.size,"~Monte~Carlo~trials:~Random~(italic(k)%~%Exp)-italic(s)[22]~~(italic(d)^ '*' ==", plot.data3[i,Threshold], ")")); i <- i + 1
label.vec[[i]] <- parse(text=paste0(MC.size,"~Monte~Carlo~trials:~Markov-italic(s)[22]~~(italic(d)^ '*' ==", plot.data3[i,Threshold], ")")); i <- i + 1
ltype.v <- c("solid", "solid")
shape.v <- c(1,16) # 3,4,8,
col.v <- (brewer.pal(n=8, name="Set1"))[c(3,5)]


# -- main plot
min.y <- 0.165
max.y <- 0.41
chosenFont <- "Times New Roman"
plot.full <- ggplot(toplot, aes(x=Threshold, y=LossRate, group=Setup)) + theme_minimal() + 
  theme(text=element_text(size=12, family=chosenFont), legend.position="bottom") + 
  labs(y="Mean loss rate (%)", x = bquote({Thresholds~italic(d)~on~italic(g)[1]})) + 
  geom_line(aes(x=Threshold, y=Mean_LossRate, colour=Setup, linetype=Setup), size=0.25, data=toplot.aggr) + 
  #geom_point(aes(x=Threshold, y=Mean_LossRate, colour=Setup, shape=Setup), size=2.5, alpha=0.8) + 
  geom_ribbon( data=toplot.aggr, aes(ymin=Mean_LossRate_Lower, ymax=Mean_LossRate_Upper, x=Threshold, fill=Setup), alpha=0.2, size=0) + 
  geom_segment(aes(x=0, xend=Threshold, y=LossRate, yend=LossRate, color=Setup), linetype="dashed", 
               data=plot.data3) + 
  geom_segment(aes(x=Threshold,xend=Threshold, y=min.y, yend=LossRate, color=Setup), linetype="dashed", 
               data=plot.data3) +  
  geom_point(aes(x=Threshold,y=LossRate, colour=Setup), data=plot.data3, size=3, shape=16) + 
  #geom_text(aes(x=Threshold,y=LossRate,label=Threshold), data=plot.data3, size=2.5, colour="gray15") + 
  # scale options
  scale_color_manual(name = "Setup", labels=label.vec, values=col.v) + 
  scale_fill_manual(name = "Setup", labels=label.vec, values=col.v) + 
  #scale_size_manual(name="Setup", labels=label.vec, values=size.v) + 
  scale_linetype_manual(name="Setup", labels=label.vec, values=ltype.v) + 
  scale_shape_manual(values=shape.v, name="Setup", labels=label.vec) + 
  guides(col=guide_legend(nrow=2)) + coord_cartesian(xlim=c(0,75),ylim=c(min.y,max.y)) +
  scale_y_continuous(limits=c(min.y,max.y), breaks=pretty_breaks(), 
                     labels=percent) + 
  scale_x_continuous(breaks=pretty_breaks())
plot.full


# -- Zoomed plot 1
# zoom bounding box
xlim <- c(8.75,12); ylim <- c(0.24265, 0.24325)
# zoomed plot definition
plot.zoom <- 
  plot.full + coord_cartesian(xlim=xlim, ylim=ylim) + 
  geom_point(aes(x=Threshold-0.02,y=LossRate, color=Setup), size=0.1, position=position_jitterdodge(0.1), alpha=0.3) +
  geom_boxplot(aes(group=Threshold, x=Threshold-0.21), outlier.shape = NA, alpha=0.2, width=0.12, colour="grey30",
               data=subset(toplot, Setup == "a_s22-Exp" & Threshold %in% 8:12)) + 
  geom_segment(aes(x=0,xend=Threshold, y=LossRate, yend=LossRate, color=Setup), linetype="dashed", data=plot.data3) + 
  geom_segment(aes(x=Threshold,xend=Threshold, y=min.y, yend=LossRate, color=Setup), linetype="dashed", data=plot.data3) +
  geom_errorbar(aes(ymin=Mean_LossRate_Lower, ymax=Mean_LossRate_Upper, x=Threshold), width=0.15, data=toplot.aggr, colour="grey30") + 
  theme(legend.position="none", axis.text.x = element_text(margin=unit(c(0,0,0,0),"mm"),size=9),
        axis.text.y=element_text(margin=unit(c(0,0,0,0),"mm"),size=9),axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(color='black', fill="white"), 
        plot.background=element_rect(color="white"),
        plot.margin = unit(c(0,0,0,0),"mm"))
plot.zoom


# -- Zoomed plot 2
# zoom bounding box
xlim <- c(20.6, 25); ylim <- c(0.1737, 0.183)
# zoomed plot definition
plot.zoom2 <- 
  plot.full + coord_cartesian(xlim=xlim, ylim=ylim, default=F) + 
  geom_point(aes(x=Threshold-0.625,y=LossRate, color=Setup), size=0.1, position=position_jitterdodge(0.18), alpha=0.5) +
  geom_boxplot(aes(group=Threshold, x=Threshold-0.45), outlier.shape = NA, alpha=0.2, width=0.2, colour="grey30",
               data=subset(toplot, Setup == "b_s22-Markov" & Threshold %in% 20:30)) + 
  geom_segment(aes(x=0,xend=Threshold, y=LossRate, yend=LossRate, color=Setup), linetype="dashed", data=plot.data3) + 
  geom_segment(aes(x=Threshold,xend=Threshold, y=min.y, yend=LossRate, color=Setup), linetype="dashed", data=plot.data3) +
  geom_errorbar(aes(ymin=Mean_LossRate_Lower, ymax=Mean_LossRate_Upper, x=Threshold), width=0.15, data=toplot.aggr, colour="grey30") + 
  theme(legend.position="none", axis.text.x = element_text(margin=unit(c(0,0,0,0),"mm"),size=9),
        axis.text.y=element_text(margin=unit(c(0,0,0,0),"mm"),size=9),axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(color='black', fill="white"), 
        plot.background=element_rect(color="white"),
        plot.margin = unit(c(0,0,0,0),"mm"))
plot.zoom2



# -- merge plots
# zoomed plot 1
use.ymin <-  0.3
use.ymax <- max.y +0.01

plot.merged <- plot.full + annotation_custom(grob = ggplotGrob(plot.zoom), xmin = 7, xmax=47,
                                             ymin = use.ymin, ymax = use.ymax)
plot.merged

# zoomed plot 2
use.ymin <-  0.19
use.ymax <- 0.30

plot.merged2 <- plot.merged + annotation_custom(grob = ggplotGrob(plot.zoom2), xmin = 20, xmax=60,
                                                ymin = use.ymin, ymax = use.ymax)
plot.merged2

dpi <- 200
# note the 9-version graph was previously set to LossThreshv2
ggsave(plot.merged2, file=paste0("LossThreshv3_CD-On_s2_rate_MC",MC.size,".png"),width=1200/dpi, height=1100/dpi,dpi=dpi)


