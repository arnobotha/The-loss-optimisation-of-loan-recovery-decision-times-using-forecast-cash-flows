# =================== Random Defaults Technique: Truncation Parameter k Estimation
# Random Defaults with Randomized (k,g1)-truncation with k drawn from a fitted distribution
# This script fits theoretical distributions to the empirical distribution of maximum delinquency observed per account

# Dependencies: 2.1b, 2.1c, 2.2a

# --- Partition data into delinquents-only (S2) and write-offs only (S3) samples, in accordance with Paper2/Thesis
samp2 <- subset(dat.use.sample, LoanAge == 1 & LoanID %in% which(vec.Del.Ind == 1) & Max_CDAlt > 0)[,list(LoanID, Datex, Max_CDAlt = as.double(Max_CDAlt), Write.Off.Ind, Sample="Delinquents")]
samp3 <- subset(dat.use.sample, LoanAge == 1 & LoanID %in% which(vec.Woff == 1)  & Max_CDAlt > 0)[,list(LoanID, Datex, Max_CDAlt = as.double(Max_CDAlt), Write.Off.Ind, Sample="Write-offs")]



# ========== Distribution Analysis + Proper Graph

# Distribution of maximum missed payments of various segments of the full sample
# Only consider delinquents and write-offs
toplot1 <- subset(dat.use.sample, LoanAge == 1 & LoanID %in% which(vec.Del.Ind == 1))[,list(LoanID, Datex, Max_CDAlt, Write.Off.Ind, Sample="Delinquents")]
toplot2 <- subset(dat.use.sample, LoanAge == 1 & LoanID %in% which(vec.Woff == 1))[,list(LoanID, Datex, Max_CDAlt, Write.Off.Ind, Sample="Write-offs")]
toplot <- rbind(toplot1,toplot2); # only bind together again for plotting purposes

chosenFont <- "Times New Roman"

# interesting differences in shape - both right-skewed in tail.
g.main1 <- ggplot(toplot, aes(x=Max_CDAlt)) + theme_minimal() + 
  geom_histogram(aes(fill=Sample,colour=Sample, y= ..density..), alpha=0.6, bins=50) + 
  scale_colour_manual(values=c("#a6bddb", "#1c9099")) + 
  scale_fill_manual(values=c("#a6bddb","#1c9099")) + 
  theme(text=element_text(size=12,family=chosenFont), legend.position="bottom") +
  labs(y="Density", x=bquote({Maximum~delinquency~per~account~-~max~g[1](italic(t))})) + 
  scale_x_continuous(breaks=pretty_breaks(), labels=comma) + 
  scale_y_continuous(breaks=pretty_breaks(), labels=comma)

g.sub1 <- ggplot(toplot, aes(x=Max_CDAlt)) + theme_minimal() + 
  geom_density(aes(colour=Sample,fill=Sample), alpha=0.5) + 
  scale_colour_manual(values=c("#a6bddb", "#1c9099")) + 
  scale_fill_manual(values=c("#a6bddb","#1c9099")) + 
  labs(y="", x="", title="Empirical Density Curves") + 
  scale_x_continuous(breaks=pretty_breaks(), labels=comma, limits=c(0, 80)) + 
  scale_y_continuous(breaks=pretty_breaks(), labels=comma) + 
  theme(text=element_text(size=12,family=chosenFont),
        legend.position="none", axis.text.x = element_text(margin=unit(c(0,0,0,0),"mm"),size=9),
        axis.text.y=element_text(margin=unit(c(0,0,0,0),"mm"),size=9),axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(color='black', fill="white"), 
        plot.background=element_rect(color="white"),
        plot.margin = unit(c(0,0,0,0),"mm"),
        plot.title = element_text(hjust=0.55, vjust=-10, margin=margin(t=-12))) 

g.full1 <- g.main1 + annotation_custom(grob = ggplotGrob(g.sub1), xmin = 40, xmax=130, ymin = 0.05, ymax = 0.15)
g.full1

dpi <- 200
ggsave(g.full1, file="Histograms and EPDFs of Max_CDAlt by sample.png",width=1200/dpi, height=1000/dpi,dpi=dpi)





# ========== Distribution Fitting

# ------- custom pdfs, cdfs, quantile functions
# define Gumbel
dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b)) # PDF
pgumbel <- function(q, a, b) exp(-exp((a-q)/b)) # CDF
qgumbel <- function(p, a, b) a-b*log(-log(p)) # inverse CDF (qunatile function)

# define dagum
ddagum <- function(x, p1, a, b) a*p1/x*( (x/b)^(a*p1) / (( (x/b)^a + 1 )^(p1+1) )) # PDF
pdagum <- function(q, p1, a, b) (1 + (q/b)^(-a))^(-p1) # CDF
qdagum <- function(p, p1, a, b) b*(p^(-1*p1) - 1)^(-1/a) # inverse CDF (quantile function)
rdagum <- function(n, p1, a, b) b*(runif(n)^(-1/p1)-1)^(-1/a)

# ---- custom functions
ks.test.vectorised <- function(x, fitted, pnames, sig.level=0.05) {
  
  test.stat <- unname(sapply(1: length(fitted), function(i) {
    est <- fitted[[i]]$estimate
    suppressWarnings(
      h.test <- do.call(ks.test, c(list(x=x, y=pnames[i]), split(unname(est), names(est))))
    )
    return (h.test$statistic)
  }))
  
  pvalues <- sapply(1: length(fitted), function(i) {
    est <- fitted[[i]]$estimate
    suppressWarnings(
      h.test <- do.call(ks.test, c(list(x=x, y=pnames[i]), split(unname(est), names(est))))
    )
    return (h.test$p.value)
  })
  
  h0.decisions <- ifelse(pvalues < sig.level, "rejected", "not rejected")
  
  return(list(PValues=pvalues, TestStats=test.stat, Decisions=h0.decisions))
}

ad.test.vectorised <- function(x, fitted, pnames, nullnames, sig.level=0.05) {
  
  test.stat <- unname(sapply(1:length(fitted), function(i) {
    est <- fitted[[i]]$estimate
    res <- do.call(ad.test, c( list(x=x, null=pnames[i], nullname=nullnames[i]), split(unname(est), names(est)) ) ) 
    return (res$statistic)
  }))
  
  pvalues <- sapply(1:length(fitted), function(i) {
    est <- fitted[[i]]$estimate
    res <- do.call(ad.test, c( list(x=x, null=pnames[i], nullname=nullnames[i]), split(unname(est), names(est)) ) ) 
    return (res$p.value)
  })
  
  h0.decisions <- ifelse(pvalues < sig.level, "rejected", "not rejected")
  
  return(list(PValues=pvalues, TestStats=test.stat, Decisions=h0.decisions))
}




# ==== S2 (delinquents only)

# - initial
hist(samp2$Max_CDAlt, breaks=50); describe(samp2$Max_CDAlt)

# - Fitting theoretical distributions to the maximum g_1 observations in S2
s1fit_g <- fitdist(samp2$Max_CDAlt, dist="gamma")
summary(s1fit_g)
s1fit_pois <- fitdist(samp2$Max_CDAlt, dist="pois")
summary(s1fit_pois)
s1fit_logn <- fitdist(samp2$Max_CDAlt, dist="lnorm")
summary(s1fit_logn)
s1fit_exp <- fitdist(samp2$Max_CDAlt, dist="exp")
summary(s1fit_exp)
s1fit_pareto <- fitdist(samp2$Max_CDAlt, dist="pareto", start=list(shape=100, scale=1) )
summary(s1fit_pareto)
s1fit_weibull <- fitdist(samp2$Max_CDAlt, dist="weibull")
summary(s1fit_weibull)
s1fit_logis <- fitdist(samp2$Max_CDAlt, dist="logis")
summary(s1fit_logis)
s1fit_cauchy <- fitdist(samp2$Max_CDAlt, dist="cauchy")
summary(s1fit_cauchy)
s1fit_dagum <- fitdist(samp2$Max_CDAlt, dist="dagum",  start=list(p1=10,a=10,b=1))
summary(s1fit_dagum)
s1fit_gumbel <- fitdist(samp2$Max_CDAlt, dist="gumbel", start=list(a=1,b=1))
summary(s1fit_gumbel)
s1fit_burr <- fitdist(samp2$Max_CDAlt, dist="burr", start=list(shape1=10, shape2=10, rate=100))
summary(s1fit_burr)

# - prepare for full comparison of goodness-of-fit tests
to.compare <- list(s1fit_g, s1fit_pois, s1fit_logn, s1fit_exp, s1fit_pareto, s1fit_weibull, s1fit_logis, s1fit_cauchy, s1fit_dagum, s1fit_gumbel, s1fit_burr)
fit.names <- c("Gamma", "Poisson", "Log-Normal", "Exponential","Pareto", "Weibull", "Logistic", "Cauchy", "Dagum", "Gumbel", "Burr")
p.names <- c("pgamma", "ppois", "plnorm", "pexp","ppareto", "pweibull", "plogis", "pcauchy", "pdagum", "pgumbel", "pburr")

# - AIC statistics
result <- gofstat(to.compare, fitnames = fit.names)
result

# - Full KS-test
ks.res <- ks.test.vectorised(x=samp2$Max_CDAlt, fitted=to.compare, pnames=p.names)
test.stats.str <- paste0(":\tD = ", format(ks.res$TestStats, digits=4, format="f"))
cat("Kolmogorov-Smirnov test for Goodness-of-fit\n", paste( str_pad(fit.names, width=max(nchar(fit.names)), side="right"),
                                                            str_pad(test.stats.str, max(nchar(test.stats.str)), side="right"), 
                                                            "\tp-value = ", formatC(ks.res$PValues, digits=4, format="f"),
                                                            ":\t", ks.res$Decisions, collapse="\n"))

# - Full AD-test
ad.res <- ad.test.vectorised(x=samp2$Max_CDAlt, fitted=to.compare, pnames=p.names, nullnames=fit.names)
test.stats.str <- paste0(":\tA^2 = ", format(ad.res$TestStats, digits=4, format="f"))
cat("Anderson-Darling test for Goodness-of-fit\n", paste( str_pad(fit.names, width=max(nchar(fit.names)), side="right"),
                                                          str_pad(test.stats.str, max(nchar(test.stats.str)), side="right"), 
                                                          "\tp-value = ", formatC(ad.res$PValues, digits=4, format="f"),
                                                          ":\t", ad.res$Decisions, collapse="\n"))


# - plot theoretical distributions across empirical curves constructed from data
# Do this only for a few candidate distributions, based on previous tests
to.compare <- list(s1fit_dagum, s1fit_logn, s1fit_pareto, s1fit_weibull, s1fit_g, s1fit_exp)
fit.names <- c("Dagum", "Log-Normal", "Pareto", "Weibull",  "Gamma", "Exponential")

plot(x=0, main=bquote({max~g(t)~plain(across)~t==list(1, ..., t[0])}))
chosenFont <- "Times New Roman"

# - compare density curves
s1g1 <- denscomp(to.compare, legendtext = fit.names, plotstyle="ggplot", breaks=50, demp=F, xlab= bquote({max~italic(g[1](t))~plain(across)~italic(t==list(1, ..., t[0]))}),
               fitcol=rev(colorRampPalette(brewer.pal(8, "Dark2"))(length(to.compare))), main="") + theme_minimal() + 
      theme(legend.position = "bottom", text=element_text(size=12,family=chosenFont)) + 
      xlim(-1,50) + coord_cartesian(xlim=c(0,50)) + guides(col=guide_legend(nrow=1, name="", title=""),
                                                           linetype=guide_legend(nrow=1, name="", title="")
                                                           )
s1g1
dpi <- 200
ggsave(s1g1, file="Max_CDAlt (delinquents) - Fitted PDFs vs Histogram.png", width=1200/dpi, height=1100/dpi,dpi=dpi)

# - compare CDF curves
s1g2 <- cdfcomp(to.compare, legendtext=fit.names, plotstyle = "ggplot", xlab = "Maximum of CD", discrete=F, 
              fitcol=rev(colorRampPalette(brewer.pal(8, "Dark2"))(length(to.compare))), main="") + theme_minimal() + theme(legend.position = "bottom")
s1g2
dpi <- 100
ggsave(s1g2, file="Max_CDAlt (delinquents) - Fitted CDFs vs eCDF.png", width=750/dpi, height=600/dpi, dpi=dpi)


# - Final choice
# Exponential
s1prep.arg3 <- split(unname(s1fit_exp$estimate), names(s1fit_exp$estimate))




# ==== S3 (write-offs only)

# - initial
hist(samp3$Max_CDAlt, breaks=50); describe(samp3$Max_CDAlt)

# - Fitting theoretical distributions to the maximum CD observations in S3
fit_g <- fitdist(samp3$Max_CDAlt, dist="gamma")
summary(fit_g)
fit_pois <- fitdist(samp3$Max_CDAlt, dist="pois")
summary(fit_pois)
fit_logn <- fitdist(samp3$Max_CDAlt, dist="lnorm")
summary(fit_logn)
fit_pareto <- fitdist(samp3$Max_CDAlt, dist="pareto", start=list(shape=100, scale=1) )
summary(fit_pareto)
fit_weibull <- fitdist(samp3$Max_CDAlt, dist="weibull")
summary(fit_weibull)
fit_logis <- fitdist(samp3$Max_CDAlt, dist="logis")
summary(fit_logis)
fit_cauchy <- fitdist(samp3$Max_CDAlt, dist="cauchy")
summary(fit_cauchy)
fit_dagum <- fitdist(samp3$Max_CDAlt, dist="dagum",  start=list(p1=1,a=10,b=10))
summary(fit_dagum)
fit_gumbel <- fitdist(samp3$Max_CDAlt, dist="gumbel", start=list(a=1,b=1))
summary(fit_gumbel)
fit_burr <- fitdist(samp3$Max_CDAlt, dist="burr", start=list(shape1=1, shape2=1, rate=1))
summary(fit_burr)

# - prepare for full comparison of goodness-of-fit tests
to.compare <- list(fit_g, fit_pois, fit_logn, fit_pareto,fit_weibull, fit_logis, fit_cauchy, fit_dagum, fit_gumbel, fit_burr)
fit.names <- c("Gamma", "Poisson", "Log-Normal", "Pareto", "Weibull", "Logistic", "Cauchy", "Dagum", "Gumbel", "Burr")
p.names <- c("pgamma", "ppois", "plnorm", "ppareto", "pweibull", "plogis", "pcauchy", "pdagum", "pgumbel", "pburr")

# - AIC statistics
result <- gofstat(to.compare, fitnames = fit.names)
result

# - Full KS-test
ks.res <- ks.test.vectorised(x=samp3$Max_CDAlt, fitted=to.compare, pnames=p.names)
test.stats.str <- paste0(":\tD = ", format(ks.res$TestStats, digits=4, format="f"))
cat("Kolmogorov-Smirnov test for Goodness-of-fit\n", paste( str_pad(fit.names, width=max(nchar(fit.names)), side="right"),
           str_pad(test.stats.str, max(nchar(test.stats.str)), side="right"), 
           "\tp-value = ", formatC(ks.res$PValues, digits=4, format="f"),
           ":\t", ks.res$Decisions, collapse="\n"))

# - Full AD-test
ad.res <- ad.test.vectorised(x=samp3$Max_CDAlt, fitted=to.compare, pnames=p.names, nullnames=fit.names)
test.stats.str <- paste0(":\tA^2 = ", format(ad.res$TestStats, digits=4, format="f"))
cat("Anderson-Darling test for Goodness-of-fit\n", paste( str_pad(fit.names, width=max(nchar(fit.names)), side="right"),
           str_pad(test.stats.str, max(nchar(test.stats.str)), side="right"), 
           "\tp-value = ", formatC(ad.res$PValues, digits=4, format="f"),
           ":\t", ad.res$Decisions, collapse="\n"))

# - Inferences:
# Poisson definitely not a good fit.


# - plot theoretical distributions across empirical curves constructed from data
# Do this only for a few candidate distributions, based on previous tests
to.compare <- list(fit_dagum, fit_burr, fit_weibull, fit_gumbel, fit_g, fit_logis)
fit.names <- c("Dagum", "Burr (Type XII)","Weibull", "Gumbel", "Gamma", "Logistic")

# - compare density curves
g1 <- denscomp(to.compare, legendtext = fit.names, plotstyle="ggplot", breaks=50, demp=F, xlab= bquote({max~italic(g[1](t))~plain(across)~italic(t==list(1, ..., t[0]))}),
         fitcol=rev(colorRampPalette(brewer.pal(8, "Dark2"))(length(to.compare))), main="") + theme_minimal() + 
  theme(legend.position = "bottom", text=element_text(size=12,family=chosenFont), legend.title = element_blank()) + 
  guides(linetype=guide_legend(nrow=1), col=guide_legend(nrow=1))
g1
dpi <- 200
ggsave(g1, file="Max_CDAlt (write-offs) - Fitted PDFs vs Histogram.png", width=1200/dpi, height=1100/dpi,dpi=dpi)

# - compare CDF curves
g2 <- cdfcomp(to.compare, legendtext=fit.names, plotstyle = "ggplot", xlab = "Maximum of CD", discrete=F, 
        fitcol=rev(colorRampPalette(brewer.pal(8, "Dark2"))(length(to.compare))), main="") + theme_minimal() + theme(legend.position = "bottom")
dpi <- 100
ggsave(g2, file="Max_CDAlt (write-offs) - Fitted CDFs vs eCDF.png", width=750/dpi, height=600/dpi, dpi=dpi)


# - Final choice
# Weibull
i <- 3
prep.arg1 <- split(unname(to.compare[[i]]$estimate), names(to.compare[[i]]$estimate))





# ========== Graphical validation tests for use in Paper2/Thesis
# --- Sample 2 (Delinquents) : Exponential

# Test against randomly generated values from Exponential distribution?
test3 <- do.call(rexp, c(list(n=1000000), s1prep.arg3))
plot(density(test3[test3<50])); 
describe(test3); hist(test3[test3<50], breaks=100)
# Mean of ~7.251; median of ~5.02

# compare to historical data
toplot1 <- subset(dat.use.sample, LoanAge == 1 & LoanID %in% which(vec.Del.Ind == 1))[,list(LoanID, Datex, Max_CDAlt, Write.Off.Ind, Sample="Delinquents")]
describe(toplot1$Max_CDAlt); hist(subset(toplot1, Max_CDAlt < 50)$Max_CDAlt, breaks=100)
# Mean of 7.246; median of 3

# --- Sample 3 (Write-offs) : Weibull
prep.arg1

# sample randomly from fitted distribution
test <- do.call(rweibull, c(list(n=1000000), prep.arg1))
plot(density(test)); describe(test)
# Mean of ~21.82; median of ~19.648

# compare to historical data
toplot2 <- subset(dat.use.sample, LoanAge == 1 & LoanID %in% which(vec.Woff == 1))[,list(LoanID, Datex, Max_CDAlt, Write.Off.Ind, Sample="Write-offs")]
describe(toplot2$Max_CDAlt); hist(toplot2$Max_CDAlt, breaks=100)
# Mean of 21.58; median of 20
