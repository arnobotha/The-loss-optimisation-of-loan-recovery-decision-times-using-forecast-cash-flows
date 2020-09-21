# =================== Sample Refinements
# This script completes some ancillary features of loans as if they had progressed to full maturity, only with missing elements.
# Note: this is not the forecasting of receipts-part. This script prepares the scene for forecasting later on.

# Dependencies: 2.1b

# --- Vectors

# - Abstract contractual term only for active loans, use observed maturity as the 'term' for the write-offs/closed accounts.
vec.Term.Treated <- copy(vec.Term)
#vec.Term.Treated[1-vec.Consider] <- vec.Maturity[1-vec.Consider] # old (possibly faulty) version
# vec.Test <- copy(vec.Term.Treated)
vec.Term.Treated[c(1-vec.Consider)==1] <- vec.Maturity[c(1-vec.Consider)==1]



# --- Matrices
# Create empty matrix by which to expand all other matrices to cater for full-length maturities of accounts
mat.empty <- matrix(NA, nrow=period.term-period, ncol=n)

# -- some previously created matrices are recreated below with the full contractual term
# in preparation for forecasting later on

# Receipt matrix
#mat.Receipt <- as.matrix(spread(dat.use.sample[,list(Account,LoanAge,Receipt)], key=Account, value = Receipt)[,-1])
mat.ReceiptAlt <- rbind(as.matrix(spread(dat.use.sample[,list(Account,LoanAge,Receipt_Alt=pmax(Receipt_Alt,0))], 
                                         key=Account, value = Receipt_Alt)[,-1]), mat.empty)
# Instalment matrix
mat.Instal <- rbind(as.matrix(spread(dat.use.sample[,list(Account,LoanAge,Install_Total)], 
                                     key=Account, value = Install_Total)[,-1]), mat.empty)

# Interest Rate Matrix
mat.IntRates <- rbind(as.matrix(spread(dat.use.sample[,list(Account,LoanAge,IntRate)],
                                       key=Account, value = IntRate)[,-1]), mat.empty)
mat.IntRates <- mat.IntRates / 100;

# Balance Matrix
mat.Bal <- rbind(as.matrix(spread(dat.use.sample[,list(Account,LoanAge, Balance)], 
                                  key=Account, value = Balance)[,-1]), mat.empty)

# Transaction-inferred instalment components
mat.Comp.Ins <- rbind(as.matrix(spread(dat.use.sample[,list(Account,LoanAge, Insurance)], 
                                       key=Account, value = Insurance)[,-1]), mat.empty)
mat.Comp.InsPayout <- rbind(as.matrix(spread(dat.use.sample[,list(Account,LoanAge, Receipt_InsPayout)], 
                                             key=Account, value = Receipt_InsPayout)[,-1]), mat.empty)
mat.Comp.Fees <- rbind(as.matrix(spread(dat.use.sample[,list(Account,LoanAge, Fees)], 
                                        key=Account, value = Fees)[,-1]), mat.empty)
mat.Comp.Int <- rbind(as.matrix(spread(dat.use.sample[,list(Account,LoanAge, Interest)], 
                                       key=Account, value = Interest)[,-1]), mat.empty)
mat.Comp.WOff <- rbind(as.matrix(spread(dat.use.sample[,list(Account,LoanAge, WOff)], 
                                        key=Account, value = WOff)[,-1]), mat.empty)



# =========== Calculate CD-measure up to observed maturity only
# -- Calculate CD (g_1: Contractual Delinquency)
mat.CDAlt.obsrvd <- calculate.CD.forData(mat.Instal, mat.ReceiptAlt, sc.Thres, period, n, method="base")
# indicate whether a loan has/had historical delinquency (including write-off, if it went straigth to write-off from performing [data peculiarity]) or not
vec.Del.Ind <- sapply(1:n, function(i,t,cd,woff) {
  # i <- 3; t<-vec.Maturity; cd<-mat.CDAlt.obsrvd; woff<-vec.Woff
  val <- ifelse(any(cd[1:t[i],i] >= 1 | woff[i]==1)==T,1,0)
  return(val)
}, t=vec.Maturity,cd=mat.CDAlt.obsrvd, woff=vec.Woff)
# describe(vec.Del.Ind)

# ======== Merge delinquency measurements back to sample
# - CDAlt (using alternative receipts)
dat.temp <- gather(as.data.frame(mat.CDAlt.obsrvd), key=LoanID, value=CDAlt) %>% na.omit() %>% mutate(LoanID = as.numeric(substr(LoanID, 2, length(LoanID)))) %>%
  group_by(LoanID) %>% mutate(LoanAge = (1:n())-1); setDT(dat.temp, key=c("LoanID","LoanAge"))
dat.use.sample <- merge(dat.use.sample, dat.temp, by=c("LoanID","LoanAge"), all.x=T) #left join (throw away delinquencies at time t=0 for now)

# - cleanup
rm(dat.temp)


# ======== Light Feature Engineering
dat.use.sample[, Max_CDAlt := max(CDAlt),by=list(Account)]
dat.use.sample[, Last_CDAlt := .SD[.N,CDAlt], by=list(Account)]






# ========================== General matrix treatments from observed maturity up to contractual term

# ------- derive other matrices
# Treat interest rate matrix to fill out remaining elements based on last observed interest rate
# This is necessary for conducting MD/DoD calculations that may use mat.Instal.Treated
mat.IntRates.Treated <- matrix(unlist(sapply(1:n, function(i,t,t.max,int) {
  
  # some latest observed interest rates are zero for some reason. Fix this via getIntRate()
  decent.rate <- getIntRate(int[1:t[i],i], t[i])
  
  prep <- c(as.vector(int[1:t[i], i]),
            as.vector(rep(decent.rate, t.max[i]-t[i])),
            as.vector(rep(NA, period.term - t.max[i]))
  )
  # fix zero elements (up to observed maturity however) with decent.rate as well
  prep[which(prep == 0)] <- decent.rate
  
  # sanity check: structure of instalment matrix
  if (NROW(prep) != period.term) cat("\nERROR: ", i)
  return(prep)
}, t=vec.Maturity, t.max=vec.Term, int=mat.IntRates)),
nrow=period.term, byrow=F)


# Treat instalment matrix to fill out remaining elements based on calculated instalment that amortizes to 0
mat.Instal.Treated <- matrix(unlist(sapply(1:n, function(i,t,t.max,instal,b,int) {

  # calculate amortizing instalment - cap at zero in case balance observed is negative
  calc.instal <- max(b[i] / ((1 - (1+int[i]/12)^ (-1* (t.max[i]-t[i]) )  ) / (int[i]/12)), 0)
  # previous instalment: instal[t[i],i]
  # create final instalment structure over time (treated) for ith account
  prep <- c(as.vector(instal[1:t[i],i]),
            as.vector(rep(calc.instal, t.max[i]-t[i])),
            as.vector(rep(NA, period.term - t.max[i]))
            )
  # sanity check: structure of instalment matrix
  if (NROW(prep) != period.term) cat("\nERROR: ", i)
  return(prep)
  }, t=vec.Maturity, t.max=vec.Term.Treated, instal=mat.Instal, 
  # get balances and int rates at last observed period
  b=mat.Bal[cbind( vec.Maturity ,1:n)], int=mat.IntRates.Treated[cbind( vec.Maturity ,1:n)])), 
nrow=period.term, byrow=F)


# ---- more transformations
mat.ReceiptAlt.FBM <- as_FBM(mat.ReceiptAlt)
vec.Instal.last <- mat.Instal.Treated[ cbind( pmin(vec.Maturity.Use+1, period.term),1:n)]
vec.Del.last <- mat.CDAlt.obsrvd[ cbind( pmin(vec.Maturity.Use+1,period.term) ,1:n)  ]




# ========== Some graphs

# Maturity analysis by loan cohort over time vs contractual term ("term gap analysis")
toplot <- subset(dat.use.sample, LoanAge == 1)[,list(LoanID, Datex, Value = Maturity, Max_CDAlt, Write.Off.Ind, NPL_Outcome, Term_Mode)]
toplot[, Month := as.Date(paste0(as.character(Datex),'01'), "%Y%m%d")]
toplot[, Median := as.double( median(Value)), by=list(Month)]
toplot[, Mean := as.double( mean(Value)), by=list(Month)]
toplot.aggr <- unique(toplot[order(Month),list(Month, a_Median_Maturity = Median, b_Mean_Maturity = Mean)])
toplot.aggr[, c_Max_Maturity := round(as.double(difftime(as.Date(date.end.long, "%Y-%m-%d"), Month, units = "days")) / 365.25 * 12) ]
toplot.aggr[, d_Contractual_Maturity := as.double(max(toplot$Term_Mode))]
toplot.aggr2 <- gather(toplot.aggr, key=Type,value=Value, -Month)
toplot <- merge(toplot, toplot.aggr, by="Month", all.x = T); rm(toplot.aggr)
label.vec <- c("Median Tenure", "Mean Tenure","Theoretical Tenure", "Contractual Term")
line.date <- "2010-10-01"
chosenFont <- "Calibri"
chosenFont <- "Times New Roman"

g2 <- ggplot(toplot.aggr2) + theme_minimal() + 
  geom_line(aes(x=Month, y=Value, colour=Type, linetype=Type, size=Type)) + 
  geom_point(aes(x=Month, y=Value, colour=Type, shape=Type, alpha=Type), size=2) + 
  labs(y="Tenure (months)", x="Origination Cohort (ccyymm)") + 
  theme(text=element_text(size=12,family=chosenFont), legend.position="bottom") +
  scale_color_brewer(palette="Spectral", label = label.vec, name="") + 
  scale_size_manual(values = c(0.75, 0.75,1, 1), label = label.vec, name="") + 
  scale_linetype_manual(values=c("dotted", "dotted", "solid", "dotted"), label = label.vec, name="") +
  scale_shape_manual(values=c(1,2,16,4), label = label.vec, name="") +
  scale_alpha_manual(values = c(1,1,0.25,1), label=label.vec, name="") + 
  scale_y_continuous(breaks= pretty_breaks(), labels=scales::comma) + 
  geom_segment(aes(x=as.Date(line.date), xend=as.Date(line.date), y=89, yend=max(toplot$Term_Mode)-2), size=1, arrow = arrow(length = unit(0.5, "cm")),
               colour="gray82")
g2

dpi <- 200
ggsave(g2, file="Maturity (FullSample) Plot by Cohort with overlays.png",width=1200/dpi, height=1000/dpi,dpi=dpi)


# - Important objective of Paper2 (showing in concept the need for forecasting observed maturity for each cohort up to contractual term)
toplot <- subset(dat.use.sample, LoanAge == 1)[,list(LoanID, Datex, Value = Maturity, Max_CDAlt, Write.Off.Ind, NPL_Outcome, Term_Mode)]
toplot[, Month := as.Date(paste0(as.character(Datex),'01'), "%Y%m%d")]
toplot[, Median := as.double( median(Value)), by=list(Month)]
toplot[, Mean := as.double( mean(Value)), by=list(Month)]
toplot.aggr <- unique(toplot[order(Month),list(Month, a_Median_Maturity = Median, b_Mean_Maturity = Mean)])
toplot.aggr[, c_Max_Maturity := round(as.double(difftime(as.Date(date.end.long, "%Y-%m-%d"), Month, units = "days")) / 365.25 * 12) ]
toplot.aggr[, d_Contractual_Maturity := as.double(max(toplot$Term_Mode))]
toplot.aggr[, a_Median_Maturity:=NULL]; toplot.aggr[, b_Mean_Maturity :=NULL]
toplot.aggr2 <- gather(toplot.aggr, key=Type,value=Value, -Month)
toplot <- merge(toplot, toplot.aggr, by="Month", all.x = T); rm(toplot.aggr)
label.vec <- c("Theoretical Tenure", "Contractual Term")
line.date <- "2010-10-01"
chosenFont <- "Calibri"
chosenFont <- "Times New Roman"

g3 <- ggplot(toplot.aggr2) + theme_minimal() + 
  geom_line(aes(x=Month, y=Value, colour=Type, linetype=Type, size=Type)) + 
  #geom_point(aes(x=Month, y=Value, colour=Type, shape=Type, alpha=Type), size=2) + 
  labs(y="Tenure (months)", x="Origination Cohort") + 
  theme(text=element_text(size=12,family=chosenFont), legend.position="bottom") +
  scale_color_manual(values=c("#1c9099","#a6bddb"), label = label.vec, name="") + 
  scale_size_manual(values = c(1, 0.75), label = label.vec, name="") + 
  scale_linetype_manual(values=c("solid", "dashed"), label = label.vec, name="") +
  scale_shape_manual(values=c(16,4), label = label.vec, name="") +
  scale_alpha_manual(values = c(1,1), label=label.vec, name="") + 
  scale_y_continuous(breaks= pretty_breaks(), labels=scales::comma) + 
  geom_segment(aes(x=as.Date(line.date), xend=as.Date(line.date), y=89, yend=max(toplot$Term_Mode)-2), size=1, arrow = arrow(length = unit(0.5, "cm")),
               colour="gray82")
g3

dpi <- 200
ggsave(g3, file="Maturity (Objective) Plot by Cohort with overlays.png",width=1200/dpi, height=1000/dpi,dpi=dpi)




# ---- draw venn diagram of 3 samples used for Article2: S1 (Full), S2 (Delinquents-only), S3 (Write-offs)
dev.off()
venn.plot <- venn.diagram(
  x = list(
    # I = c(1:60, 61:105, 106:140, 141:160, 166:175, 176:180, 181:205,
    #       206:220),
    # II = c(531:605, 476:530, 336:375, 376:405, 181:205, 206:220, 166:175,
    #        176:180),
    # III = c(61:105, 106:140, 181:205, 206:220, 221:285, 286:335, 336:375,
    #         376:405)
    # A = c(1:n),
    # B = c(1:NROW(which(vec.Del.Ind==1))),
    # C = c(1:NROW(which(vec.Woff==1)))
    A = c(1:n),
    B = c((which(vec.Del.Ind==1))),
    C = c((which(vec.Woff==1)))
    
  ),
  category.names = c(
    expression( bold('S'['1']) ),
    expression( bold('S'['2']) ),
    expression( bold('S'['3']) )
  ),
  filename = "VennSamples.png",
  output = TRUE,
  height = 3000,
  width = 3000,
  resolution = 300,
  compression = 'lzw',
  units = 'px',
  lwd = 6,
  lty = 'blank',
  fill = c('#e0ecf4', '#9ebcda', '#8856a7'),
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 2.5,
  cat.fontface = "bold",
  #cat.default.pos = "outer",
  cat.pos = c(0, 0, 180),
  cat.dist = c(-0.025, -0.05, -0.005),
  #cat.fontfamily = "sans",
  rotation = 1,
  alpha=c(0.75,0.75,0.75)
)

grid.draw(venn.plot)



