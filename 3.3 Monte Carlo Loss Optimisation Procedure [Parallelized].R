# =================== MONTE CARLO VARIANCE DIAGNOSTICS OF LOSS OPTIMISATION PROCEDURE
# Similar to 3.2 this interactive script also performs the loss optimisation procedure for a pre-specified set of settings, 
# but it does so iteratively within a Monte Carlo framework by re-estimating forecasts repeatedly using different seed values. 
# This script is also run experimentaly across a few toggable settings (recorded in comments)
# The goal is to examine the variance underlying the loss curves, using a Monte Carlo-based method to estimate the variance.
# This is achieved by producing a graph as the end-result wherein we show the "mean" loss curves, accompanied by estimated confidence bands

# Three main experimental settings are used (s_ij).
# Note that the same sample is used for both forecast training (i) and the optimisation procedure (j):
# - s11-full: Random defaults and Markovian defaults techniques.
# - s22-delinquents: Random defaults and Markovian defaults techniques.
# - s33-WOff: Random defaults and Markovian defaults techniques.

# Direct dependencies include: 2.1b, 2.1c


# ====== 1. Initializing some global parameters

# -- Pointers
mat.ReceiptAlt.FBM <- as_FBM(mat.ReceiptAlt)
mat.Instal.Use <- (mat.Instal.Treated)
mat.IntRate.Use <- (mat.IntRates.Treated)
vec.Maturity.Use <- vec.Maturity # observed loan tenures
vec.Mat.Use <- vec.Term.Treated # contractual term

# -- Iteration Parameter
num.thresholds <-75; #number of thresholds (# of times loss calculations, given threshold d)
it.vec <- 1:num.thresholds

# -- General iteration parameters
it.max <- NROW(it.vec)
first.iter <- T

# -- Interest Rates & Loss Rates: Specification & Calculation
i.alt <- 0.07; #risk-free rate (effective rate)
i_p.alt <- ((1+i.alt)^(1/12) - 1)*12; #risk-free rate (nominal rate)
Arrears.LossRate <- 0.7;
Outstanding.LossRate <- 0.4;

# -- Select thresholds for g1-measure
vec.k.CD <- seq(0, (num.thresholds-1))



# ====== 2. LOSS ASSESSMENT: Iterative function definitions (to be run in parallel)

# - main function for assessing the portfolio loss at a specified threshold (using only the g1 delinquency measure)
coreJob_CD <- function(vec.Maturity.Use, vec.Mat.Use, mat.Receipt.Use, mat.Instal.Use, mat.IntRate.Use, sc.Thres, period.term, n, vec.LoanAmt, 
                    it, num.thresholds, d.CD, mat.CD.Use, vec.Consider,i_p.alt, Arrears.LossRate, Outstanding.LossRate) {
  
  
  cat(paste0("\n1) Threshold [",it," of ",num.thresholds,"]: (", Sys.time(),") Assessing portfolio loss at threshold .. "),
      file="assesslog.txt", append=TRUE)
  
  # ---- Total Loss across given threshold (d.CD)
  
  # - get default start times of first episode (if multiple exist), given threshold d, otherwise return -1 to indicate a performing loan
  vec.default.start_first.CD <- sapply(1:n, default.start.first.v2, thres.d=d.CD, del.mat=mat.CD.Use, t=vec.Mat.Use)
  
  # - get (g,d)-defaulting account indices across measure, given current thresholds
  def.CD <- which(vec.default.start_first.CD >= 0 & vec.Consider==1)
  
  # - get (g,d)-performing account indices across measure, given current thresholds
  perf.CD <- which(vec.default.start_first.CD < 0 & vec.Consider==1)
  
  # - calculate final maturity as either contractual term / maturity or default time, given (g,d)-default 
  # for use in discounting and other loss calculations
  vec.maturity.CD <- copy(vec.Mat.Use)
  vec.maturity.CD[def.CD] <- vec.default.start_first.CD[def.CD]
  
  # - Calculate NPV of receipts, given maturity and relevant receipts
  vec.ReceiptsPV.CD <- sapply(1:n, function(i,r,t) {
    if (t[i] > 0) {
      val <- sum( r[1:t[i], i] * (1+i_p.alt/12)^(-1*1:(t[i]) ) )
    } else {
      val <- 0
    }
    return (val)
  }, r=mat.Receipt.Use, t=vec.maturity.CD)
  
  # - calculate NPV of arrears, given maturity, relevant instalments and relevant receipts
  vec.ArrearsPV.CD <- sapply(1:n, function(i,ins,r,t) {
    if (t[i] > 0) {
      val <- sum( ins[1:t[i],i] * (1+i_p.alt/12)^(-1*1:(t[i]) ) ) - r[i]
    } else {
      val <- 0
    }
    return (val)
  }, ins=mat.Instal.Use, r=vec.ReceiptsPV.CD, t=vec.maturity.CD)
  
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
  
  # - calculate losses as weighted combination between arrears and expected balance, and associated loss rates
  vec.Losses.CD <- pmax(vec.ArrearsPV.CD*Arrears.LossRate + vec.ExpBalance.CD*Outstanding.LossRate, 0)
  
  # - calculate actual balance [ancillary information]
  vec.bal.CD <- pmax(vec.ArrearsPV.CD + vec.ExpBalance.CD, 0)
  
  # --- curate some vectors for reporting/graphing purposes - remove accounts not to be considered
  {
    vec.bal.CD[which(vec.Consider==0)] <- NA
    
    # - zero the loss if a particular account is not to be considered [sampling]
    vec.Losses.CD[which(vec.Consider==0)] <- 0
  }
  
  
  # ---------- Concatenate relevant information, including profit/loss calculations for optimisation
  dat.EL.core <- rbind(data.table(Iteration=it, Measure="CD", MeasureName ="g1: CD", Threshold=d.CD,
                                  Vol_Perf=length(perf.CD),Vol_Def=length(def.CD),
                                  Bal_Perf = sum(vec.bal.CD[perf.CD], na.rm = T), Bal_Def = sum(vec.bal.CD[def.CD], na.rm = T),
                                  Loss=sum(vec.Losses.CD, na.rm = T))
  )
  
  cat(paste0("\n2) Threshold [",it," of ",num.thresholds,"]: (", Sys.time(),") Loss assessed! "),
      file="assesslog.txt", append=TRUE)
  
  return (dat.EL.core)
}





# ====== 3. OUTER FUNCTION TO BE CALLED WITHIN MONTE CARLO FRAMEWORK 

# - Function call scenarios, including possible values for input variables or other pointers:
# forecast.tech: 'random', 'Markov'
# forecast.tech='random': iterate for prob.b: [prob.b.full], [prob.b.del], [prob.b.woff]
# forecast.tech='random': iterate for fitted.params: [prep.arg1] (Weibull), [s1prep.arg3] (Exponential)
# forecast.tech='random': iterate for trunc.dist: 'exp', 'weibull'
# forecast.tech='Markov': iterate for p.trans: [p.DelStates], [p.DelStates.del], [p.DelStates.woff]
outerJob <- function(forecast.tech='random', prob.b, fitted.params, trunc.dist, p.trans, states, seed.value, n, 
                     vec.Maturity.Use, vec.Mat.Use, mat.ReceiptAlt, mat.ReceiptAlt.FBM, mat.Instal.Use, period.term, mat.CDAlt.obsrvd, sc.Thres,
                     mat.IntRate.Use, vec.LoanAmt, num.thresholds, it.max, vec.k.CD, outer.iter, outer.iter.max, vec.Consider,
                     i_p.alt, Arrears.LossRate, Outstanding.LossRate, vec.Instal.last, vec.Del.last) {
  
  ptm <- proc.time() #IGNORE: for computation time calculation
  
  # ===================  Forecast cash flows to full contractual term
    
  cat(paste0("\n\n1)[",outer.iter," of ",outer.iter.max,"] Monte Carlo simulations (", Sys.time(),"). Forecasting .. "),
      file="montecarlo_log.txt", append=TRUE)
  
  if (forecast.tech == "random") {
    
    # --- Random Defaults
    mat.Receipt.Use <- forecastJob_RandomDefaults(seed.value=seed.value*(n-1), n=n, fitted.params=fitted.params, prob.b=prob.b, trunc.dist=trunc.dist,
                                                  vec.Maturity.Use=vec.Maturity.Use, vec.Mat.Use=vec.Mat.Use, mat.ReceiptAlt=mat.ReceiptAlt, 
                                                  mat.Instal.Use=mat.Instal.Use, period.term=period.term, mat.CDAlt.obsrvd=mat.CDAlt.obsrvd, sc.Thres=sc.Thres)
  } else if (forecast.tech == "Markov") {

    # --- Markovian Defaults
    mat.Receipt.Use <- forecastJob_MarkovDefaults_p(seed.value=seed.value*(n-1), p.trans=p.trans, states=states, n=n, vec.Maturity.Use=vec.Maturity.Use, 
                               vec.Mat.Use=vec.Mat.Use, mat.ReceiptAlt=mat.ReceiptAlt.FBM, vec.Instal.last=vec.Instal.last, 
                               period.term=period.term, vec.Del.last=vec.Del.last, sc.Thres)
  }
  
  # ---- Calculate g1-measure up to full contractual term using forecasts
  mat.CD.Use <- calculate.CD.forData(mat.Instal.Use, mat.Receipt.Use, sc.Thres, period.term, n, method="base")
  

  # ===================  Portfolio loss assessment across thresholds: parallelized loop
  
  cat(paste0("\n2)[",outer.iter," of ",outer.iter.max,"] Monte Carlo simulations (", Sys.time(),"): Receipts forecasted and delinquency assessed. Assessing portfolio losses across thresholds .. "),
      file="montecarlo_log.txt", append=TRUE)
  
  cat(paste("New Job (", Sys.time(),"): Assessing delinquency and profitability of given portfolio across various thresholds",sep=''),
      file="assesslog.txt", append=FALSE)
  
  
  # using foreach() from foreach package for distributing loop iterations across registered threads: remarkable improvement in run time
  # Note: need to import all custom functions used within the loss assessment.
  dat.EL.par <- foreach(it=1:it.max, .combine='rbind', .verbose=F, .inorder=F, .packages ='data.table',
                        .export=c('default.start.first.v2', 'coreJob_CD')) %dopar%
    {
      
      dat.EL.core <- coreJob_CD(vec.Maturity.Use=vec.Maturity.Use, vec.Mat.Use=vec.Mat.Use,
                             mat.Receipt.Use=mat.Receipt.Use, mat.Instal.Use=mat.Instal.Use,
                             mat.IntRate.Use=mat.IntRate.Use, sc.Thres=sc.Thres, period.term=period.term, n=n, 
                             vec.LoanAmt=vec.LoanAmt, it=it, num.thresholds=num.thresholds, d.CD=vec.k.CD[it], mat.CD.Use=mat.CD.Use,
                             vec.Consider=vec.Consider, i_p.alt=i_p.alt, Arrears.LossRate=Arrears.LossRate, Outstanding.LossRate=Outstanding.LossRate)
    }
  
  elapsed <- proc.time() - ptm
  
  cat(paste0("\n3)[",outer.iter," of ",outer.iter.max,"] Monte Carlo simulations (", Sys.time(),"): Done! Elapsed time: ", round(elapsed['elapsed'],digits=0), " seconds"),
      file="montecarlo_log.txt", append=TRUE)
  
  # - last data preparation
  dat.EL.par[, Iteration := outer.iter]
  setDT(dat.EL.par, key=c("Iteration","Measure","Threshold"))

  return(dat.EL.par)
  
}







# ====== 4. MONTE CARLO FRAMEWORK

# ---- Parametrisation history:
# 1a -- Random defaults, trained from s1, applied on s1 [v2_5j(i)]: prob.b.full, s1prep.arg3, "exp", vec.Consider <- rep(1,n)
# 1b -- Markovian defaults, trained from s1, applied on s1 [v2_6a(i)]: p.DelStates, vec.Consider <- rep(1,n)
# 2a -- Random defaults, trained from s2, applied on s2 [v2_5k(ii)]: prob.b.del, s1prep.arg3, "exp", vec.Consider[which(vec.Del.Ind==0)] <- 0
# 2b -- Markovian defaults, trained from s2, applied on s2 [v2_6b(ii)]: p.DelStates.del, vec.Consider[which(vec.Del.Ind==0)] <- 0
# 3a -- Random defaults, trained from s3, applied on s3 [v2_5f(ix)]: prob.b.woff, prep.arg1, "weibull", vec.Consider[which(vec.Woff==0)] <- 0
# 3b -- Markovian defaults, trained from s3, applied on s3 [v2_6c(iii)]: p.DelStates.woff, vec.Consider[which(vec.Woff==0)] <- 0


# ---- Parametrisation

# -- Monte Carlo parameters
outer.iter.max <- 500
cpu.threads <- 6

# -- Forecasting technique parameters
forecast.tech.use <- "random" # random, Markov
# - Random Defaults
prob.b.use <- prob.b.full # prob.b.full, prob.b.del, prob.b.woff
fitted.params.use <- s1prep.arg3 # prep.arg1 (Weibull for write-offs), s1prep.arg3 (Exponential for full or delinquents-only)
trunc.dist.use <- "exp" # exp, weibull
# - Markovian Defaults
p.trans.use <- p.DelStates # p.DelStates, p.DelStates.del, p.DelStates.woff

# -- File name settings for storing results
it.name <- "v1_1a"

# -- toggle sampling for loss optimisation (i: full sample [Lowest risk]; ii; delinquents-only [Medium risk]; iii: write-offs-only [Highest risk])
# This is an artificial proxy for having portfolios on with different risk profiles, on which we may optimise the recovery decision
vec.Consider <- rep(1,n) # -- i: switches to indicate use full sample
#vec.Consider[which(vec.Del.Ind==0)] <- 0 # -- ii: only consider delinquent (including write-offs) loans, by switching off the rest
#vec.Consider[which(vec.Woff==0)] <- 0 # -- iii: only consider those written-off loans, by switching off the rest


# ---- Main Loop

ptm <- proc.time() #IGNORE: for computation time calculation

cl.port <- makeCluster(cpu.threads) # number of threads to register in the OS for this procedure (used in outerJob)
registerDoParallel(cl.port)

cat(paste("New Job (", Sys.time(),"): ", outer.iter.max, " Monte Carlo simulations. Experiment series ", it.name, ".",sep=''),
    file="montecarlo_log.txt", append=FALSE)



# Note: need to import all custom functions used within the loss assessment.
dat.EL.MC <- foreach(outer.iter=1:outer.iter.max, .combine='rbind', .verbose=T, .inorder=F, .packages=c('data.table','foreach'),
                      .export=c('default.start.first.v2', 'coreJob_CD', 'forecastJob_RandomDefaults', 
                                'forecastJob_MarkovDefaults_p','outerJob')) %do%
  {
    dat.EL.outercore <- outerJob(forecast.tech=forecast.tech.use, prob.b=prob.b.use, fitted.params=fitted.params.use, trunc.dist=trunc.dist.use,
             p.trans=p.trans.use, states=states, seed.value=(outer.iter), n=n, vec.Maturity.Use=vec.Maturity.Use, 
             vec.Mat.Use=vec.Mat.Use, mat.ReceiptAlt=mat.ReceiptAlt, mat.ReceiptAlt.FBM=mat.ReceiptAlt.FBM, mat.Instal.Use=mat.Instal.Use, period.term=period.term,
             mat.CDAlt.obsrvd=mat.CDAlt.obsrvd, sc.Thres=sc.Thres, mat.IntRate.Use=mat.IntRate.Use, vec.LoanAmt=vec.LoanAmt,
             num.thresholds=num.thresholds, it.max=it.max, vec.k.CD=vec.k.CD, outer.iter=outer.iter, outer.iter.max=outer.iter.max,
             vec.Consider=vec.Consider, i_p.alt=i_p.alt, Arrears.LossRate=Arrears.LossRate, Outstanding.LossRate=Outstanding.LossRate,
             vec.Instal.last, vec.Del.last)
  }

# - last data preparation
setDT(dat.EL.MC, key=c("Iteration","Measure","Threshold"))
# - zip and save optimisation results to disk
pack.ffdf(paste0("MonteCarlo_EL", outer.iter.max, "-",it.name),dat.EL.MC)

elapsed <- proc.time() - ptm
elapsed

stopCluster(cl.port) # release threads back to the OS

cat(paste0("\n\nEnd of Job (", Sys.time(),"): ", outer.iter.max, " Monte Carlo simulations done. Elapsed time: ", round(elapsed['elapsed'],digits=0), " seconds"),
    file="montecarlo_log.txt", append=TRUE)

