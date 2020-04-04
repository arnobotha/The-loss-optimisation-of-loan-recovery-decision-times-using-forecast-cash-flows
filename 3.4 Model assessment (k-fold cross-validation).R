# =================== FORECAST ERROR ACCURACY ASSESSMENT: k-fold Cross-validation approach
# A script for implemneting a few measures of forecast error in a systematic way.
# Accordingly, forecasts are retrained and assessed on validation dataset, following a k-fold cross-validation approach.
# For expediency, only estimates from S1 (Full sample) are tested here, with S2 and S3 disregarded.

# Direct dependencies include scripts: 2.1b, 2.1c, 2.2b




# ====== 1a. Initializing some global parameters

# -- Pointers
mat.ReceiptAlt.FBM <- as_FBM(mat.ReceiptAlt)
mat.Instal.Use <- mat.Instal.Treated 
mat.IntRate.Use <- mat.IntRates.Treated
vec.Maturity.Use <- vec.Maturity # observed loan tenures
vec.Mat.Use <- vec.Term.Treated # contractual term


# -- Forecasting Technique parameters

# - Random Defaults (using S1 settings)
fitted.params.use <- s1prep.arg3 #  Exponential for full or delinquents-only
trunc.dist.use <- "exp" # exp, weibull

# - Create payment indicator function for Random Defaults' estimation part
mat.repay.ind.T <- matrix(unlist(sapply(1:n, function(i,t,ins,r) {
  val <- c(ifelse(r[1:t[i],i] >= ins[1:t[i],i], 1, 0),
           as.vector(rep(NA, period.term - t[i]))
  )
  return(val)
}, t=vec.Maturity, ins=mat.Instal.Use, r=mat.ReceiptAlt.FBM)),
nrow=period.term, byrow=F)


# -- toggle sampling for loss optimisation (i: full sample; ii; delinquents-only; iii: write-offs-only)
vec.Consider <- rep(1,n) # -- switches to indicate use full sample


# -- k-fold cross validation parameters
k <- 5
consider.n <- sum(vec.Consider)
# - for those data points under consideration, draw a randomized order of indexes (same size as 'full data')
set.seed(123)
ind.random <- sample(x=which(vec.Consider==1),size=as.integer(consider.n))

# -- filename
fil.name <- paste0("k",k,"-crossValidation-tests")




# ====== 1b. Preliminary parameter stability calculation

# Prepare estimates from the full dataset (no training/validation scheme)
# This is to enable comparisons with estimates trained from each kth fold later on
# in assessing parameter sensitivity

# --- Random defaults
prob.b.AllData <- sum(sapply(which(vec.Consider==1), function(i,m,t) { sum(m[1:t[i],i]) }, m=mat.repay.ind.T, t=vec.Maturity)) / 
  sum(vec.Maturity[which(vec.Consider==1)])

# --- Markovian defaults
# - prepare base matrix
mat.DelStates.T.AllData <- matrix(unlist(sapply(which(vec.Consider==1), function(i,t,del,woff) {
  # i<-3; ins<-mat.Instal; t<-vec.Maturity; r<- mat.ReceiptAlt
  prep.val <- case_when(
    del[1:t[i], i] == 0 ~ 1, # "Active/CD0"
    del[1:t[i], i] == 1 ~ 2, # CD1
    del[1:t[i], i] == 2 ~ 3, # CD2
    del[1:t[i], i] == 3 ~ 4, # CD3
    del[1:t[i], i] == 4 ~ 5, # CD4
    del[1:t[i], i] == 5 ~ 6, # CD5
    del[1:t[i], i] >= 6 ~ 7  # CD6+
  )
  # now curate last element as state 8 when written-off 
  prep.val[t[i]] <- woff[i] * 8 + (1-woff[i]) * prep.val[t[i]] # write-off is state 8
  
  val <- c(prep.val, 
           as.vector(rep(NA, period.term - t[i]))
  )
  return(val)
}, t=vec.Maturity.Use, del=mat.CDAlt.obsrvd, woff=vec.Woff)),
nrow=period.term, byrow=F)

# ---- MLE Estimates for transition matrix
p.DelStates.T.AllData <- matrix(nrow=NROW(states), ncol=NROW(states), 0)
colnames(p.DelStates.T.AllData) <- states.T; rownames(p.DelStates.T.AllData) <- states.T
for (t in 1:(period.term-1)) {
  # t <- 1
  for (state.i in 1:NROW(states.T)) {
    for (state.j in 1:NROW(states.T)) {
      p.DelStates.T.AllData[state.i, state.j] <- p.DelStates.T.AllData[state.i, state.j] + NROW(which(mat.DelStates.T.AllData[t,] ==state.i & mat.DelStates.T.AllData[t+1,] ==state.j))   
    }
  }
}
for (state.i in 1:NROW(states.T)) p.DelStates.T.AllData[state.i, ] <- p.DelStates.T.AllData[state.i, ] / sum(p.DelStates.T.AllData[state.i, ])




# ====== 2. Panel of model asssessment tests

# -- custom function that implements one run of model training, forecasting, and final assessment - all within
# a particular kth fold within a bigger k-fold cross-validation setup
# Assessment contains a broad panel of measures, including forecast error on the observation-level and portfolio-level,
# and stability of parameter estimates
forecastAssess.kcrossvalid <- function(it.k=1) {
  
  # -- Set training indices depending on current k
   # first, determine k-th swathe of indices to use as validation observations
  ind.fold <- ((it.k-1)*as.integer(consider.n/k) + 1) : ((it.k)*as.integer(consider.n/k))
  # secondly, set the complement hereof as training observations
  ind.train <- ind.random[-ind.fold]
  
  # - create validation flag vector, pointing to indices of records in original dataset meant for validation
  vec.Validation <- rep(0, consider.n)
  vec.Validation[setdiff(ind.random, ind.train)] <- 1
  # - unset validation flags in vec.Consider as training flag vector
  vec.Consider.train <- copy(vec.Consider)
  vec.Consider.train[-ind.train] <- 0
  
  # -- Artifically set the loan ages to 1 of those validation records we wish to forecast.
  vec.Maturity.Valid <- copy(vec.Maturity.Use)
  vec.Maturity.Valid[which(vec.Validation==1)] <- 1
  
  # - Also reset the 'last observed' vectors (instalment and g1-delinquency) to values appropriate in vec.Maturity.Valid
  vec.Instal.Valid <- mat.Instal.Use[1, 1:n]
  vec.Del.Valid <- rep(0,n)
  
  
  # ====== Training Forecasts ..
  
  # --- Random Defaults
  {
    # ==== Training .. 
    prob.b.use <- sum(sapply(which(vec.Consider.train==1), function(i,m,t) { sum(m[1:t[i],i]) }, m=mat.repay.ind.T, t=vec.Maturity)) / 
      sum(vec.Maturity[which(vec.Consider.train==1)])
    prob.b.use
    
    # ==== Forecasting .. 
    # use vec.Maturity.Valid instead of vec.Maturity.Use, thereby forcing forecasts from loan age 1 up to the full contractual term
    mat.Receipt.Test.Random <- forecastJob_RandomDefaults(seed.value=1, n=n, fitted.params=fitted.params.use, prob.b=prob.b.use, trunc.dist=trunc.dist.use,
                                                          vec.Maturity.Use=vec.Maturity.Valid, vec.Mat.Use=vec.Mat.Use, mat.ReceiptAlt=mat.ReceiptAlt, 
                                                          mat.Instal.Use=mat.Instal.Use, period.term=period.term, mat.CDAlt.obsrvd=mat.CDAlt.obsrvd, sc.Thres=sc.Thres)
    
    # ==== Delinquency calcualtion
    mat.CD.Test.Random <- calculate.CD.forData(mat.Instal.Use, mat.Receipt.Test.Random, sc.Thres, period.term, n, method="base")
  }
  
  
  # --- Markovian Defaults
  {  
    # ==== Training .. 
    
    # ---- prepare base matrix
    mat.DelStates.T <- matrix(unlist(sapply(which(vec.Consider.train==1), function(i,t,del,woff) {
      # i<-3; ins<-mat.Instal; t<-vec.Maturity; r<- mat.ReceiptAlt
      prep.val <- case_when(
        del[1:t[i], i] == 0 ~ 1, # "Active/CD0"
        del[1:t[i], i] == 1 ~ 2, # CD1
        del[1:t[i], i] == 2 ~ 3, # CD2
        del[1:t[i], i] == 3 ~ 4, # CD3
        del[1:t[i], i] == 4 ~ 5, # CD4
        del[1:t[i], i] == 5 ~ 6, # CD5
        del[1:t[i], i] >= 6 ~ 7  # CD6+
      )
      # now curate last element as state 8 when written-off 
      prep.val[t[i]] <- woff[i] * 8 + (1-woff[i]) * prep.val[t[i]] # write-off is state 8
      
      val <- c(prep.val, 
               as.vector(rep(NA, period.term - t[i]))
      )
      return(val)
    }, t=vec.Maturity.Use, del=mat.CDAlt.obsrvd, woff=vec.Woff)),
    nrow=period.term, byrow=F)
    
    # ---- MLE Estimates for transition matrix
    states.T <- c("1.UpToDate", "2.CD1", "3.CD2", "4.CD3", "5.CD4", "6.CD5", "7.CD6+", "8.WOff")
    p.DelStates.T <- matrix(nrow=NROW(states), ncol=NROW(states), 0)
    colnames(p.DelStates.T) <- states.T; rownames(p.DelStates.T) <- states.T
    for (t in 1:(period.term-1)) {
      # t <- 1
      for (state.i in 1:NROW(states.T)) {
        for (state.j in 1:NROW(states.T)) {
          p.DelStates.T[state.i, state.j] <- p.DelStates.T[state.i, state.j] + NROW(which(mat.DelStates.T[t,] ==state.i & mat.DelStates.T[t+1,] ==state.j))   
        }
      }
    }
    for (state.i in 1:NROW(states.T)) p.DelStates.T[state.i, ] <- p.DelStates.T[state.i, ] / sum(p.DelStates.T[state.i, ])
    
    # ==== Forecasting .. 
    # use vec.Maturity.Valid instead of vec.Maturity.Use, thereby forcing forecasts from loan age 1 up to the full contractual term  
    mat.Receipt.Test.Markov <- forecastJob_MarkovDefaults_p(seed.value=1, p.trans=p.DelStates.T, states=states.T, n=n, 
                                                            vec.Maturity.Use=vec.Maturity.Valid, vec.Mat.Use=vec.Mat.Use, 
                                                            mat.ReceiptAlt=mat.ReceiptAlt.FBM, vec.Instal.last=vec.Instal.Valid, 
                                                            period.term=period.term, vec.Del.last=vec.Del.Valid, sc.Thres, createOwnPar = T)
    
    # ==== Delinquency calcualtion
    mat.CD.Test.Markov <- calculate.CD.forData(mat.Instal.Use, mat.Receipt.Test.Markov, sc.Thres, period.term, n, method="base")  
  }
  
  
  
  # ========= Assessing Accuracy of Forecasts ..
  
  # ====== Accuracy Test 1: Mean Absolute Error (MAE)
  # https://en.wikipedia.org/wiki/Mean_absolute_error
  
  # MAE: portfolio-level sum of summed absolute errors
  T.MAE <- function(mat.Receipt.Given) {
    errs <- unlist(sapply(which(vec.Validation==1), function(loan) {
      # absolute forecast errors
      out <- abs(mat.ReceiptAlt[2:vec.Maturity.Use[loan],loan] - mat.Receipt.Given[2:vec.Maturity.Use[loan],loan])
      # data error correction (cases that only have 1-month observed history should be excluded, therefore we zero the forecast error)
      if (any(is.na(out))) out <- NA
      return(out)
    }))
    mean(errs, na.rm = T)
  }
  
  # - MAE
  vec.Test1a <- c(T.MAE(mat.Receipt.Test.Random), T.MAE(mat.Receipt.Test.Markov))
  # value closer to 0 is better
  
  # - expressed as a proportion of the mean instalment
  # get all instalments across all periods for validation records only
  installs <- unlist(sapply(which(vec.Validation==1), function(i) {
    mat.Instal.Use[2:vec.Maturity.Use[i], i ]
  }))
  n.eff <- length(installs) - sum(is.na(installs))
  # mean
  instal.mean <- mean(installs, na.rm=T)
  # variance
  instal.var <- var(installs, na.rm=T)
  
  # 95% confidence bands around estimated mean instalment
  instal.CI <- c( instal.mean - 1.96* sqrt(instal.var/ n.eff ),
     instal.mean + 1.96*sqrt(instal.var/ n.eff ) )
  # quite narrow bands, so we are quite confidence about our estimate of the mean.
  
  # - Scaled MAE (Normalised)
  vec.Test1b <- vec.Test1a / instal.mean
  # value closer to 0 is better, but this shows scale as well.
  
  
  
  # ====== Accuracy Test 2: Root Mean Squared Error (RMSE) / Mean Squared Prediction Error (MSPE) actually in our hold-out case
  # https://en.wikipedia.org/wiki/Root-mean-square_deviation
  
  # RMSE: portfolio-level sum of summed squared errors
  T.RMSE <- function(mat.Receipt.Given) {
    sq.errs <-  unlist(sapply(which(vec.Validation==1), function(loan) {
      # Loan-level squared errors
      # start at 2 since that is when forecasting actually happens
      ( mat.ReceiptAlt[2:vec.Maturity.Use[loan],loan] - mat.Receipt.Given[2:vec.Maturity.Use[loan],loan] )^2
    }))
    
    sqrt(mean(sq.errs,na.rm=T))
  }
  
  # - RMSE
  vec.Test2a <- c(T.RMSE(mat.Receipt.Test.Random), T.RMSE(mat.Receipt.Test.Markov))
  # value closer to 0 is better
  
  # - Normalized RMSE
  vec.Test2b <- vec.Test2a / instal.mean
  
  
  
  # ====== Accuracy Test 3: Delinquency Forecast Error (DFE) (own measure)
  
  # Mean or Median
  T.AFE <- function(mat.CD.Test.Given, cent.meas="Mean") {
    errs <- unlist(sapply(which(vec.Validation==1), function(loan, d, d.f, t) {
      d.f[2:t[loan],loan] - d[2:t[loan],loan]
    }, d=mat.CDAlt.obsrvd, d.f=mat.CD.Test.Given, t=vec.Maturity.Use
    ))
    if (cent.meas=="Mean") {
      mean(errs, na.rm=T)
    } else if (cent.meas=="Median") {
      median(errs, na.rm=T)    
    } else { # just output general summary statistics
      describe(errs)
    }
  }
  
  vec.Test3a <- c(T.AFE(mat.CD.Test.Random), T.AFE(mat.CD.Test.Markov))
  vec.Test3b <- c(T.AFE(mat.CD.Test.Random, cent.meas="Median"), T.AFE(mat.CD.Test.Markov, cent.meas="Median"))
  
  
  
  # ====== Accuracy Test 4: Portfolio Arrears Rate (PAR) (own measure)
  
  T.PAR <- function(mat.Receipt.Given, disc.rate=((1.07)^(1/12) - 1)*12 ) {
    
    PVs <- unlist(sapply(which(vec.Validation==1), function(loan, r, t, d, ins) {
      (ins[2:t[loan],loan] - r[2:t[loan],loan] ) * (1+d/12)^(-1*1:(t[loan]-1) )
    }, r=mat.Receipt.Given, t=vec.Maturity.Use, d=disc.rate, ins=mat.Instal.Use
    )) 
    
    sum(PVs, na.rm=T) / sum(vec.LoanAmt[which(vec.Consider==1)], na.rm=T)
  }
  
  vec.Test4a <- c(T.PAR(mat.Receipt.Test.Random), T.PAR(mat.Receipt.Test.Markov))
  Port.Arr.Rate.Train <- T.PAR(mat.ReceiptAlt) # Arrears rate of the whole portfolio, for comparison:
  
  
  
  # ====== Parameter Variance Test 1:
  # --- Random Defaults
  # Percentage difference in parameters
  Params.Random <- c( prob.b.AllData, prob.b.use)
  Params.Random.Diff <- (Params.Random[2] - Params.Random[1]) / Params.Random[1]
  
  # --- Markovian Defaults
  Params.Markov <- list(AllData=p.DelStates.T.AllData, Test=p.DelStates.T)
  Params.Markov.Diff <- mean(as.vector((Params.Markov$Test[1:7,] - Params.Markov$AllData[1:7,]) / 
                                         Params.Markov$AllData[1:7,]), na.rm=T)
  
  
  # ====== Combining test results
  
  dat.accur.tests.core <- data.table(k=it.k, Model=c("Random Defaults","Markovian Defaults"), MAE=vec.Test1a, MAE.Norm=vec.Test1b,
               RMSE=vec.Test2a, RMSE.Norm=vec.Test2b, DFE.Mean=vec.Test3a, DFE.Median=vec.Test3b,
               PAR=vec.Test4a, PAR.Train=Port.Arr.Rate.Train, Mean.Param.Diff=c(Params.Random.Diff,Params.Markov.Diff))

  return (dat.accur.tests.core)

}



# ====== 3. k-fold Cross-validation framework

# - loop over each kth fold, and assess
dat.accur.tests <- foreach(it.k=1:k, .combine='rbind', .verbose=T, .inorder=T, .packages=c('data.table','foreach'),
                     .export=c('forecastJob_RandomDefaults', 'forecastJob_MarkovDefaults_p', 'calculate.CD.forData',
                               'forecastAssess.kcrossvalid')) %do%
  {
    dat.accur.tests.core <- forecastAssess.kcrossvalid(it.k=it.k)
  }

# - zip and store performance results to disk
pack.ffdf(paste0("Data/",fil.name),dat.accur.tests)


# - aggregate metrics across folds, for presentation
dat.accur.tests.overall <- dat.accur.tests[,list(Avg.MAE= comma(mean(MAE, na.rm=T)),
                      Avg.MAE.Norm=paste0( round(mean(MAE.Norm, na.rm=T)*100, digits=1), "%"),
                      Avg.RMSE=mean(RMSE, na.rm=T),
                      Avg.RMSE.NORM=percent(mean(RMSE.Norm, na.rm=T)),
                      Avg.DFE.Mean = round(mean(DFE.Mean, na.rm=T),digits=1),
                      Avg.DFE.Median = round(mean(DFE.Median, na.rm=T),digits=1),
                      Avg.PAR = paste0( round(mean(PAR, na.rm=T)*100, digits=3), "%"),
                      AVG.PAR.Train = percent(mean(PAR.Train, na.rm=T)),
                      AVG.Mean.Param.Diff = percent(mean(Mean.Param.Diff, na.rm=T))
                      ), by=list(Model)]
dat.accur.tests.overall # final results
