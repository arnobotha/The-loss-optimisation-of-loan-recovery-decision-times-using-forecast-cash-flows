# =================== Experimental Optimisation Procedure across various scenarios [Parallelized]

# This interactively-run experimental script uses the output from a single forecasting technique
# (specified via a pointer to a particular receipt matrix), as parametrised in a certain way, 
# along with using an accompanying delinquency matrix (computed on the now-completed portfolio),
# to assess the portfolio loss across all specified thresholds, using 3 pre-selected delinquency measures (g1, g2, g3).
# The various ways in which this script was run is recorded in the comments below. 

# Direct dependencies include: 2.1b, 2.1c, 2.2a, 2.2b
# Interactive dependencies include forecast receipt matrices that are all created experimentally in:
#   - Random defaults: 2.2c receipt forecasting
#   - Markovian defaults: 2.3 Estimating MLEs for transition rates and subsequent receipt forecasting

# A NOTE ON vec.Consider: This vector is interactively reset in this script as an artificial proxy for having different loan portfolios
# It was previous set to exclude closed cases when forecasting delinquency.



# ====== 0a. Initializing some parameters

# -- toggle the following matrices to specific versions created from data
mat.Instal.Use <- mat.Instal.Treated
mat.IntRate.Use <- mat.IntRates.Treated
vec.Maturity.Use <- vec.Maturity # observed loan tenures
vec.Mat.Use <- vec.Term.Treated # contractual term

# -- Pointer to a specific forecast matrix, selected from a wider experimental basis
#mat.Receipt.Use <- (mat.ReceiptAlt) # untreated receipt matrix

# - Index of experiments: uncomment one
#### RANDOM DEFAULTS TECHNIQUE [defined in 2.2c ]
####### Using Exp distribution on S2 (delinquents)
mat.Receipt.Use <- mat.ReceiptAlt.Treated7j # for v2_5j(i-iii) (repayment probability estimated from full sample; truncation parameter k drawn randomly from Exponential distribution fitted on Max_CDAlt from delinquents-only sample)
#mat.Receipt.Use <- mat.ReceiptAlt.Treated7k # for v2_5k(i-iii) (repayment probability estimated from delinquents-only sample; truncation parameter k drawn randomly from Exponential distribution fitted on Max_CDAlt from delinquents-only sample)
#mat.Receipt.Use <- mat.ReceiptAlt.Treated7l # for v2_5l(i-iii) (repayment probability estimated from write-offs sample; truncation parameter k drawn randomly from Exponential distribution fitted on Max_CDAlt from delinquents-only sample)
####### Using Weibull distribution on S3 (write-offs)
#mat.Receipt.Use <- mat.ReceiptAlt.Treated7d # for v2_5d(i-vi) (repayment probability estimated from full sample; truncation parameter k drawn randomly from Weibull distribution fitted on Max_CDAlt from write-offs-only sample)
#mat.Receipt.Use <- mat.ReceiptAlt.Treated7e # for v2_5e(i-ix) (repayment probability estimated from delinquents-only sample; truncation parameter k drawn randomly from Weibull distribution fitted on Max_CDAlt from write-offs-only sample)
#mat.Receipt.Use <- mat.ReceiptAlt.Treated7f # for v2_5f(i-ix) (repayment probability estimated from write-offs sample; with truncation parameter k drawn randomly from Weibull distribution fitted on Max_CDAlt from write-offs-only sample)

#### MARKOVIAN DEFAULTS TECHNIQUE [defined in 2.3]
#mat.Receipt.Use <- mat.ReceiptAlt.Treated8a # for v2_6a(i-iii) (treated with multi-state Markovian defaults technique with parameter estimates from full sample)
#mat.Receipt.Use <- mat.ReceiptAlt.Treated8b # for v2_6b(i-iii) (treated with multi-state Markovian defaults technique with parameter estimates from delinquents-only sample)
#mat.Receipt.Use <- mat.ReceiptAlt.Treated8c # for v2_6c(i-iii) (treated with multi-state Markovian defaults technique with parameter estimates from write-offs-only sample)


# ====== 0b. Calculate Delinquency Measures up to full contractual term, using forecast receipts

# -- Calculate CD (g1: Contractual Delinquency)
mat.CD.Use <- calculate.CD.forData(mat.Instal.Use, mat.Receipt.Use, sc.Thres, period.term, n, method="base")

# -- Calculate MD/DoD (g2/g3: Macaulay Duration Index (MD) Measure | Degree of Delinquency (DoD) Measure)
calc.results <- calculate.MDoD.forData(mat.Instal.Use, mat.Receipt.Use, vec.LoanAmt, vec.Mat.Use, 
                                       n, mat.IntRate.Use, vec.DoD.lambda)
mat.MD.Use <- calc.results$MD
mat.DoD.Use <- calc.results$DoD
rm(calc.results) #an optimization, reduces memory usage



# ====== 0c. Loss assessment, risk profile selection, and iteration parameters for subsequent optimisation

# -- toggle sampling for loss optimisation (i: full sample [Lowest risk]; ii; delinquents-only [Medium risk]; iii: write-offs-only [Highest risk])
# This is an artificial proxy for having portfolios on with different risk profiles, on which we may optimise the recovery decision
vec.Consider <- rep(1,n) # -- i: switches to indicate use full sample
 #vec.Consider[which(vec.Del.Ind==0)] <- 0 # -- ii: only consider delinquent (including write-offs) loans, by switching off the rest
 #vec.Consider[which(vec.Woff==0)] <- 0 # -- iii: only consider those written-off loans, by switching off the rest


# - script saving options (may overwrite previously saved data if not careful)
inner.name <- "LossThresh" # experiment theme name
it.name <- "v2_5j(i)-excl_closed" #iteration name
plot.name <- paste0(inner.name, it.name) # full name

# -- Iteration Parameter
num.thresholds <-168; #number of delinquency thresholds (essentially the number of times that loss is calculated)
it.vec <- 1:num.thresholds # iteration vector

# - General parameters
it.max <- NROW(it.vec)
first.iter <- T

# -- Interest Rates & Loss Rates: Specification & Calculation
i.alt <- 0.07; #risk-free rate (effective rate)
i_p.alt <- ((1+i.alt)^(1/12) - 1)*12; #risk-free rate (nominal rate)
Arrears.LossRate <- 0.7;
Outstanding.LossRate <- 0.4;





# ====== 1. Select thresholds for each Delinquency Measure
# NOTE: Some of the ranges of chosen thresholds may need to be tweaked, especially for MD and Dod measures,
# depending on the chosen portfolio on which optimisation is performed, as well as the risk level underyling receipt forecasts.
# Failure to tweak may lead to false conclusions and/or local optima in results.
# Discretionary tweaking itself is currently performed on a trial-and-error basis of running the optimisation multiple times using different ranges
# and trying to isolate the 'neighbourhood' where loss optima seemingly occurs.

# -- CD
vec.k.CD <- seq(0, (num.thresholds-1)) # chose all integer-valued thresholds, no tweaking necessary

# -- MD
max.thres <- max(quantile(mat.MD.Use[!is.na(mat.MD.Use)], 0.985)) + 1
vec.k.MD <- seq(1, ifelse(is.na(max.thres), 5, max(max.thres, 5)),length.out = num.thresholds) # normal selection
vec.k.MD <- c(seq(1, 2.5, length.out = 50),seq(2.51, max.thres, length.out = num.thresholds-50)) # for full sample/delinquents
#vec.k.MD <- c(seq(1, 4.5, length.out = 50),seq(4.51, max.thres, length.out = num.thresholds-50)) # for write-offs
#plot(vec.k.MD)

# -- DoD
max.thres <- max(quantile(mat.DoD.Use[!is.na(mat.DoD.Use)], 0.985)) + 1
vec.k.DoD <- seq(1, ifelse(is.na(max.thres), 5, max(max.thres, 5)),length.out = num.thresholds) # normal selection
vec.k.DoD <- c(seq(1, 2.5, length.out = 50),seq(2.51, max.thres, length.out = num.thresholds-50)) # for full sample/delinquents
#vec.k.DoD <- c(seq(1, 4.5, length.out = 50),seq(4.51, max.thres, length.out = num.thresholds-50)) # for write-offs
#plot(vec.k.DoD)




# ====== 2. LOSS ASSESSMENT: Iterative function definitions (to be run in parallel)

# - main function for assessing the portfolio loss at a specified threshold (one for each of the g1, g2, g3 delinquency measures)
coreJob <- function(vec.Maturity.Use, vec.Mat.Use, mat.Receipt.Use, mat.Instal.Use, mat.IntRate.Use, sc.Thres, period.term, n, vec.LoanAmt, 
                    vec.DoD.lambda, it, num.thresholds, d.CD, d.MD, d.DoD, mat.CD.Use, mat.MD.Use, mat.DoD.Use) {
  
  cat(paste0("\n 1)[",it," of ",num.thresholds,"] Loss assessments .. "),
      file="assesslog.txt", append=TRUE)
  
  # ---- Total Loss across given Threshold (d.CD, d.MD, d.DoD)
    
    # - get default start times of first episode (if multiple exist), given threshold d, otherwise return -1 to indicate a performing loan
    # uses custom function default.start.first.v2()
    vec.default.start_first.CD <- sapply(1:n, default.start.first.v2, thres.d=d.CD, del.mat=mat.CD.Use, t=vec.Mat.Use)
    vec.default.start_first.MD <- sapply(1:n, default.start.first.v2, thres.d=d.MD, del.mat=mat.MD.Use, t=vec.Mat.Use)
    vec.default.start_first.DoD <- sapply(1:n, default.start.first.v2, thres.d=d.DoD, del.mat=mat.DoD.Use, t=vec.Mat.Use)
    
    # - get (g,d)-defaulting account indices across measure, given current thresholds
    def.CD <- which(vec.default.start_first.CD >= 0 & vec.Consider==1)
    def.MD <- which(vec.default.start_first.MD >= 0 & vec.Consider==1)
    def.DoD <- which(vec.default.start_first.DoD >= 0 & vec.Consider==1)
    
    # - get (g,d)-performing account indices across measure, given current thresholds
    perf.CD <- which(vec.default.start_first.CD < 0 & vec.Consider==1)
    perf.MD <- which(vec.default.start_first.MD < 0 & vec.Consider==1)
    perf.DoD <- which(vec.default.start_first.DoD < 0 & vec.Consider==1)
    
    # - calculate final maturity as either contractual term / maturity or default time, given (g,d)-default 
    # for use in discounting and other loss calculations
    vec.maturity.CD <- copy(vec.Mat.Use)
    vec.maturity.CD[def.CD] <- vec.default.start_first.CD[def.CD]
    vec.maturity.MD <- copy(vec.Mat.Use)
    vec.maturity.MD[def.MD] <- vec.default.start_first.MD[def.MD]
    vec.maturity.DoD <- copy(vec.Mat.Use)
    vec.maturity.DoD[def.DoD] <- vec.default.start_first.DoD[def.DoD]
    
    
    # - Calculate NPV of receipts, given maturity and relevant receipts
    vec.ReceiptsPV.CD <- sapply(1:n, function(i,r,t) {
      if (t[i] > 0) {
        val <- sum( r[1:t[i], i] * (1+i_p.alt/12)^(-1*1:(t[i]) ) )
      } else {
        val <- 0
      }
      return (val)
    }, r=mat.Receipt.Use, t=vec.maturity.CD)
    vec.ReceiptsPV.MD <- sapply(1:n, function(i,r,t) {
      if (t[i] > 0) {
        val <- sum( r[1:t[i], i] * (1+i_p.alt/12)^(-1*1:(t[i]) ) )
      } else {
        val <- 0
      }
      return (val)
    }, r=mat.Receipt.Use, t=vec.maturity.MD)
    vec.ReceiptsPV.DoD <- sapply(1:n, function(i,r,t) {
      if (t[i] > 0) {
        val <- sum( r[1:t[i], i] * (1+i_p.alt/12)^(-1*1:(t[i]) ) )
      } else {
        val <- 0
      }
      return (val)
    }, r=mat.Receipt.Use, t=vec.maturity.DoD)
    
    
    # - calculate NPV of arrears, given maturity, relevant instalments and relevant receipts
    vec.ArrearsPV.CD <- sapply(1:n, function(i,ins,r,t) {
      if (t[i] > 0) {
        val <- sum( ins[1:t[i],i] * (1+i_p.alt/12)^(-1*1:(t[i]) ) ) - r[i]
      } else {
        val <- 0
      }
      return (val)
    }, ins=mat.Instal.Use, r=vec.ReceiptsPV.CD, t=vec.maturity.CD)
    vec.ArrearsPV.MD <- sapply(1:n, function(i,ins,r,t) {
      if (t[i] > 0) {
        val <- sum( ins[1:t[i],i] * (1+i_p.alt/12)^(-1*1:(t[i]) ) ) - r[i]
      } else {
        val <- 0
      }
      return (val)
    }, ins=mat.Instal.Use, r=vec.ReceiptsPV.MD, t=vec.maturity.MD)
    vec.ArrearsPV.DoD <- sapply(1:n, function(i,ins,r,t) {
      if (t[i] > 0) {
        val <- sum( ins[1:t[i],i] * (1+i_p.alt/12)^(-1*1:(t[i]) ) ) - r[i]
      } else {
        val <- 0
      }
      return (val)
    }, ins=mat.Instal.Use, r=vec.ReceiptsPV.DoD, t=vec.maturity.DoD)
    
    
    # - calculate expected balance, given tenure at (g,d)-default, resulting remaining tenure, instalments, and interest rates
    vec.ExpBalance.CD <- sapply(1:n, function(i,ins,intr,t,tt) {
      if (t[i] < tt[i]) {
        val <- sum( ins[((t[i]+1):tt[i]),i] * (1+intr[((t[i]+1):tt[i]),i]/12)^(-1*1:(tt[i] - t[i]) ) ) ;
      } else {
        val <- 0
      }
      # discount to origination
      val <- val *  (1+i_p.alt/12)^(-1*t[i] )
      return (val)
    }, ins=mat.Instal.Use, intr=mat.IntRate.Use, t=vec.maturity.CD, tt=vec.Mat.Use)
    vec.ExpBalance.MD <- sapply(1:n, function(i,ins,intr,t,tt) {
      if (t[i] < tt[i]) {
        val <- sum( ins[((t[i]+1):tt[i]),i] * (1+intr[((t[i]+1):tt[i]),i]/12)^(-1*1:(tt[i] - t[i]) ) ) ;
      } else {
        val <- 0
      }
      # discount to origination
      val <- val *  (1+i_p.alt/12)^(-1*t[i] )
      return (val)
    }, ins=mat.Instal.Use, intr=mat.IntRate.Use, t=vec.maturity.MD, tt=vec.Mat.Use)
    vec.ExpBalance.DoD <- sapply(1:n, function(i,ins,intr,t,tt) {
      if (t[i] < tt[i]) {
        val <- sum( ins[((t[i]+1):tt[i]),i] * (1+intr[((t[i]+1):tt[i]),i]/12)^(-1*1:(tt[i] - t[i]) ) ) ;
      } else {
        val <- 0
      }
      # discount to origination
      val <- val *  (1+i_p.alt/12)^(-1*t[i] )
      return (val)
    }, ins=mat.Instal.Use, intr=mat.IntRate.Use, t=vec.maturity.DoD, tt=vec.Mat.Use)
    
    
    # - calculate losses as weighted combination between arrears and expected balance, and associated loss rates
    vec.Losses.CD <- pmax(vec.ArrearsPV.CD*Arrears.LossRate + vec.ExpBalance.CD*Outstanding.LossRate, 0)
    vec.Losses.MD <- pmax(vec.ArrearsPV.MD*Arrears.LossRate + vec.ExpBalance.MD*Outstanding.LossRate, 0)
    vec.Losses.DoD <- pmax(vec.ArrearsPV.DoD*Arrears.LossRate + vec.ExpBalance.DoD*Outstanding.LossRate, 0)
    
    
    # - calculate actual balance [ancillary information]
    vec.bal.CD <- pmax(vec.ArrearsPV.CD + vec.ExpBalance.CD, 0)
    vec.bal.MD <- pmax(vec.ArrearsPV.MD + vec.ExpBalance.MD, 0)
    vec.bal.DoD <- pmax(vec.ArrearsPV.DoD + vec.ExpBalance.DoD, 0)
    
    
    # --- curate some vectors for reporting/graphing purposes - remove accounts not to be considered
    {
      vec.bal.CD[which(vec.Consider==0)] <- NA
      vec.bal.MD[which(vec.Consider==0)] <- NA
      vec.bal.DoD[which(vec.Consider==0)] <- NA
      
      # - zero the loss if a particular account is not to be considered [sampling]
      vec.Losses.CD[which(vec.Consider==0)] <- 0
      vec.Losses.MD[which(vec.Consider==0)] <- 0
      vec.Losses.DoD[which(vec.Consider==0)] <- 0
    }
    
    
    # ---------- Concatenate relevant information, including profit/loss calculations for optimisation
    dat.EL.core <- rbind(data.table(Measure="CD", MeasureName ="g1: CD", Threshold=d.CD,
                                       Vol_Perf=length(perf.CD),Vol_Def=length(def.CD),
                                       Bal_Perf = sum(vec.bal.CD[perf.CD], na.rm = T), Bal_Def = sum(vec.bal.CD[def.CD], na.rm = T),
                                       Loss=sum(vec.Losses.CD, na.rm = T)),
                            data.table(Measure="MD", MeasureName ="g2: MD", Threshold=d.MD,
                                       Vol_Perf=length(perf.MD),Vol_Def=length(def.MD),
                                       Bal_Perf = sum(vec.bal.MD[perf.MD], na.rm = T), Bal_Def = sum(vec.bal.MD[def.MD], na.rm = T),
                                       Loss=sum(vec.Losses.MD, na.rm = T)),
                            data.table(Measure="DoD", MeasureName ="g3: DoD", Threshold=d.DoD,
                                       Vol_Perf=length(perf.DoD),Vol_Def=length(def.DoD),
                                       Bal_Perf = sum(vec.bal.DoD[perf.DoD], na.rm = T), Bal_Def = sum(vec.bal.DoD[def.DoD], na.rm = T),
                                       Loss=sum(vec.Losses.DoD, na.rm = T))
    )
    
    cat(paste0("\n\t 2)[",it," of ",num.thresholds,"] Loss assessed! "),
        file="assesslog.txt", append=TRUE)
  
  return (dat.EL.core)
}




# ====== 3. OPTIMISATION PROCEDURE: Iterating across chosen thresholds (g1, g2, g3 delinquency measures)

ptm <- proc.time() #IGNORE: for computation time calculation

cat(paste("New Job: Assessing delinquency and profitability of given portfolio across various thresholds",sep=''),
    file="assesslog.txt", append=FALSE)

# -- parallelization parameters
cl.port<-makeCluster(6) # number of threads to register in the OS for this procedure
registerDoParallel(cl.port)

# using foreach() from foreach package for distributing loop iterations across registered threads: remarkable improvement in run time
# Note: need to import all custom functions used within the loss assessment.
dat.EL.par <- foreach(it=1:it.max, .combine='rbind', .verbose=T, .inorder=F, .packages ='data.table',
                   .export=c('default.start.first.v2', 'coreJob')) %dopar%
                   {
                     
                     dat.EL.core <- coreJob(vec.Maturity.Use=vec.Maturity.Use, vec.Mat.Use=vec.Mat.Use,
                                            mat.Receipt.Use=mat.Receipt.Use, mat.Instal.Use=mat.Instal.Use, mat.IntRate.Use=mat.IntRate.Use,
                                            sc.Thres=sc.Thres, period.term=period.term, n=n, vec.LoanAmt=vec.LoanAmt, vec.DoD.lambda=vec.DoD.lambda, 
                                            it=it, num.thresholds=num.thresholds,
                                            d.CD=vec.k.CD[it], d.MD=vec.k.MD[it], d.DoD=vec.k.DoD[it],
                                            mat.CD.Use=mat.CD.Use, mat.MD.Use=mat.MD.Use, mat.DoD.Use=mat.DoD.Use)
                   }

stopCluster(cl.port) # release threads back to the OS

# - last data preparation
setDT(dat.EL.par, key=c("Measure","Threshold"))
# - zip and save optimisation results to disk
pack.ffdf(paste0("EL",it.name),dat.EL.par)




# =========== OPTIMISATION RESULTS: Isolated optima and graphs

# -- Balance and Volume graphs across chosen threshold range, using g1-measure [CD]
# just for inspection purposes
ex.g <- "CD"
toplot <- gather(dat.EL.par[Measure==ex.g, list(Threshold, Bal_Perf, Bal_Def)], key=Type, value=Value, -Threshold)
label.vec <- c("Defaulting Balance", "Performing Balance")
ggplot(toplot, aes(x=Threshold, group=Type)) + theme_minimal() + 
  geom_bar(aes(x = Threshold, y=Value, fill = Type), position="fill", stat="identity") + 
  scale_y_continuous(breaks=pretty_breaks(), labels=percent) + theme(legend.position = "bottom") + 
  labs(y="Proportionate Balances (%)",x=paste0("Default Threshold (",ex.g,")")) + 
  scale_fill_manual(values=c("paleturquoise4","paleturquoise"),labels=label.vec, name="")

toplot <- gather(dat.EL.par[Measure==ex.g,list(Threshold, Vol_Perf, Vol_Def)], key=Type, value=Value, -Threshold)
label.vec <- c("Defaulting Volume", "Performing Volume")
ggplot(toplot, aes(x=Threshold, group=Type)) + theme_minimal() + 
  geom_bar(aes(x = Threshold, y=Value, fill = Type), position="fill", stat="identity") + 
  scale_y_continuous(breaks=pretty_breaks(), labels=percent) + theme(legend.position = "bottom") + 
  labs(y="Proportionate Volume (%)",x=paste0("Default Threshold (",ex.g,")")) + 
  scale_fill_manual(values=c("paleturquoise4","paleturquoise"),labels=label.vec, name="")

# save graphs to disk
dpi <- 110
ggsave(g.Bal, file=paste0("EL-Balances-",it.name,".png"),width=600/dpi, height=450/dpi,dpi=dpi)
ggsave(g.Vol, file=paste0("EL-Volumes-",it.name,".png"),width=600/dpi, height=450/dpi,dpi=dpi)


# -- Optimisation results across chosen threshold range, by delinquency measure
g <- ggplot(dat.EL.par, aes(x=Threshold, y=Loss)) + 
  geom_point(aes(x=Threshold,y=Loss, color=MeasureName, shape=MeasureName), size=1.5) +
  geom_line(aes(x=Threshold, y=Loss, color=MeasureName), size = 1) + 
  labs(y="PV of Losses (R)",x="Thresholds (w)") + theme_minimal() + 
  theme(text=element_text(family="Calibri", size=12),
        legend.position="bottom") + 
  scale_color_economist(name="Delinquency Measure",guide = guide_legend(ncol=2)) +
  scale_shape_manual(values=c(1,16,8), 
                     name="Delinquency Measure",guide = guide_legend(ncol=2)) +
  scale_y_continuous(breaks= pretty_breaks(), labels=unit_format(unit="b", scale=0.000000001))
g


# -- Optimisation results across chosen threshold range, for the g1_measure [CD] only
toplot <- dat.EL.par[Measure==ex.g, list(MeasureName, Threshold, Loss)]
g.CDOnly <- ggplot(toplot, aes(x=Threshold, y=Loss)) + 
  geom_point(aes(x=Threshold,y=Loss, color=MeasureName, shape=MeasureName), size=1.5) +
  geom_line(aes(x=Threshold, y=Loss, color=MeasureName), size = 1) + 
  labs(y="PV of Losses (R)",x="Thresholds (w)") + theme_minimal() + 
  theme(text=element_text(family="Calibri", size=12),
        legend.position="bottom") + 
  scale_color_economist(name="Delinquency Measure",guide = guide_legend(ncol=2)) +
  scale_shape_manual(values=c(1,16,8), 
                     name="Delinquency Measure",guide = guide_legend(ncol=2)) +
  scale_y_continuous(breaks= pretty_breaks(), labels=unit_format(unit="b", scale=0.000000001))
g.CDOnly


# save graphs to disk
ggsave(g, file=paste0(plot.name,".png"),width=600/100, height=450/100,dpi=100)
ggsave(g.CDOnly, file=paste0(plot.name,"-CD.png"),width=600/100, height=450/100,dpi=100)


# -- Loss-optimal thresholds found across chosen threshold range, by delinquency measure
minima <- function(given) {
  dat.min <- given[Measure=="CD", list(MeasureName, Threshold, Loss)]
  min.pos <- Position(function(x) x==min(dat.min$Loss), dat.min$Loss)
  cat("CD: Minimum Loss of", comma(min(dat.min$Loss)),"at threshold d =", dat.min[min.pos, Threshold], "at position", min.pos, "in threshold vector")
  dat.min <- given[Measure=="MD", list(MeasureName, Threshold, Loss)]
  min.pos <- Position(function(x) x==min(dat.min$Loss), dat.min$Loss)
  cat("\nMD: Minimum Loss of", comma(min(dat.min$Loss)),"at threshold d =", dat.min[min.pos, Threshold], "at position", min.pos, "in threshold vector")
  dat.min <- given[Measure=="DoD", list(MeasureName, Threshold, Loss)]
  min.pos <- Position(function(x) x==min(dat.min$Loss), dat.min$Loss)
  cat("\nDoD: Minimum Loss of", comma(min(dat.min$Loss)),"at threshold d =", dat.min[min.pos, Threshold], "at position", min.pos, "in threshold vector")  
}
minima(dat.EL.par)

proc.time() - ptm #IGNORE: computation time taken



