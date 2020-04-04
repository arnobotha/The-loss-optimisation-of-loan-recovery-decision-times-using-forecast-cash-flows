# =============== DELINQUENCY MEASURES
# This script contains function definitions for constructing the g_1, g_2, and g_3 delinquency measures
# It also contains variants hereof designed to work with real-world data



# ========= Function Declarations

# ==== Function: Calculate CD (g1: Contractual Delinquency)
# - Inputs the following:
# 1) vec.Instal: a vector of fixed instalments, one for each loan account within a loan portfolio
# 2) mat.Receipt: a matrix of time-indexed cashflows/receipts, each column corresponding to each account
# 3) sc.Thres: a scalar threshold above which the repayment ratio is considered current and beneath which it is considered delinquent
# 4) period: a scalar indicating the contractual period of all loans
# 5) n: the number of loans within the portfolio
# 6) method: arguments include either 'base' or 'simple':
#    'base' -> construct the robust g_1 measure
#    'simple' -> construct the simpler arrears / instalment ratio and take the ceiling thereof (not used)
calculate.CD <- function(vec.Instal, mat.Receipt, sc.Thres, period, n, method="base") {
  
  # - prepare various working matrices, to be filled later
  mat.CD <- matrix(-1, nrow=period+1, ncol=n); #include time 0
  mat.CD.c <- mat.CD;
  mat.CD.d <- mat.CD;
  mat.CD.m <- mat.CD;
  
  # - two methods are implemented: 'base' (default) and 'simple'
  if (method == "base") {
    
    # - calculate repayment ratio h_t at each time point for each loan
    mat.RepayRatio <- sapply(1:n, function(i,R,I) {
      val <- c(0,R[1:period,i] / I[i])
      return(val)
    }, R=mat.Receipt, I=vec.Instal)
    
    # - create a matrix containing Boolean-values of the decision function d_1(t,i), at every period t of each loan i
    mat.CD.d <- ifelse(mat.RepayRatio < sc.Thres, 1, 0)
    
    # - create a matrix containing integer-values of the m(t,i) function, i.e., the reduction in accrued delinquency at every period t of each loan i
    mat.CD.m <- floor(mat.RepayRatio / sc.Thres)*(1-mat.CD.d) - 1
    
    # - let g_1(t) at t=0 (row=1) be equal to zero, by mathematical design
    mat.CD[1,] <- 0
    
    # - finally, create the g_1 measure at each period t simultaneously across all loans (vectorised approach)
    for (ii in 2:(period+1)) {
      # - fill in at period t=ii the Boolean-values of the decision function d_2, across all loans
      mat.CD.c[ii,] <- ifelse(mat.CD[ii-1,] == 0, 1, 0)
      # - now construct g_1 at period t=11 simultaneously across all loans
      mat.CD[ii,] <-  pmax(0, mat.CD.d[ii,]*mat.CD.c[ii,] + (1-mat.CD.c[ii,])*(mat.CD[ii-1,] - mat.CD.m[ii,]))
    }
    
  } else if (method == "simple") {
    
    # - convert the previous instalment vector into a time-indexed matrix of instalments (repeat the value across all periods)
    mat.Instal <- sapply(1:n, function(i,I) {
      return(c(0,rep(I[i],period)))
    },I=vec.Instal);
    
    # - compute the differences between instalments and receipts at each period for each loan
    mat.Diff <- mat.Instal - rbind(rep(0,n),mat.Receipt);
    
    # - create the cumulative sum of these differences, i.e., the arrears balance
    mat.Arrears <- sapply(1:n, function(i,y) {
      bal <- cumsum(y[,i])
      return(bal)
    }, y=mat.Diff)
    
    # - simply divide the accumulated arrears with the fixed instalment, and take the ceiling hereof as the number of payments in arrears
    mat.CD <- ceiling(mat.Arrears %*% diag(1/vec.Instal))
  }
  
  return (mat.CD)
}



# ==== Function: Calculate CD (g1: Contractual Delinquency), assuming variable instalments from to real data
# - Inputs the following:
# 1) mat.Instal: a matrix of time-indexed variable instalments, each column corresponding to each account within a wider loan portfolio
# 2) mat.Receipt: a matrix of time-indexed cashflows/receipts, each column corresponding to each account
# 3) sc.Thres: a scalar threshold above which the repayment ratio is considered current and beneath which it is considered delinquent
# 4) period: a scalar indicating the contractual period of all loans
# 5) n: the number of loans within the portfolio
# 6) method: arguments include either 'base' or 'simple':
#    'base' -> construct the robust g_1 measure
#    'simple' -> construct the simpler arrears / instalment ratio and take the ceiling thereof (not used)
calculate.CD.forData <- function(mat.Instal, mat.Receipt, sc.Thres, period, n, method="base") {
  
  # - prepare various working matrices, to be filled later
  mat.CD <- matrix(-1, nrow=period+1, ncol=n); #include time 0
  mat.CD.c <- mat.CD;
  mat.CD.d <- mat.CD;
  mat.CD.m <- mat.CD;
  
  # - two methods are implemented: 'base' (default) and 'simple'
  if (method == "base") {
    
    # - calculate repayment ratio h_t at each time point for each loan
    mat.RepayRatio <- sapply(1:n, function(i,R,I) {
      
      ratios <- ifelse(I[1:period,i] > 0,  R[1:period,i] / I[1:period,i], 
                       ifelse(R[1:period,i] > 0, R[1:period,i] / 0.0001, sc.Thres)
                       );
      val <- c(0,ratios) # pad with zeros at origination (t=0)
      return(val)
    }, R=mat.Receipt, I=mat.Instal)
    
    # - create a matrix containing Boolean-values of the decision function d_1(t,i), at every period t of each loan i
    mat.CD.d <- ifelse(mat.RepayRatio < sc.Thres, 1, 0)
    
    # - create a matrix containing integer-values of the m(t,i) function, i.e., the reduction in accrued delinquency at every period t of each loan i
    mat.CD.m <- floor(mat.RepayRatio / sc.Thres)*(1-mat.CD.d) - 1
    
    # - let g_1(t) at t=0 (row=1) be equal to zero, by mathematical design
    mat.CD[1,] <- 0
    
    # - finally, create the g_1 measure at each period t simultaneously across all loans (vectorised approach)
    for (ii in 2:(period+1)) {
      # - fill in at period t=ii the Boolean-values of the decision function d_2, across all loans
      mat.CD.c[ii,] <- ifelse(mat.CD[ii-1,] == 0, 1, 0)
      # - now construct g_1 at period t=11 simultaneously across all loans
      mat.CD[ii,] <-  pmax(0, mat.CD.d[ii,]*mat.CD.c[ii,] + (1-mat.CD.c[ii,])*(mat.CD[ii-1,] - mat.CD.m[ii,]))
    }
    
  } else if (method == "simple") {
    
    # - compute the differences between instalments and receipts at each period for each loan
    mat.Diff <- mat.Instal - rbind(rep(0,n),mat.Receipt);
    
    # - create the cumulative sum of these differences, i.e., the arrears balance
    mat.Arrears <- sapply(1:n, function(i,y) {
      bal <- cumsum(y[,i])
      return(bal)
    }, y=mat.Diff)
    
    # - simply divide the accumulated arrears with the fixed instalment, and take the ceiling hereof as the number of payments in arrears
    mat.CD <- ceiling(mat.Arrears %*% diag(1/vec.Instal))
  }
  
  return (mat.CD)
}



# ==== Function: Calculates g2: Macaulay Duration Index (MD-measure), g3: Degree of Delinquency (DoD-measure)
# - Inputs the following:
# 1) vec.Instal: a vector of fixed instalments, one for each loan account within a loan portfolio
# 2) mat.Receipt: a matrix of time-indexed cashflows/receipts (starting at t=0), each column corresponding to the history of each account
# 3) vec.Principal: a vector of loan amounts/principals, constituting individual accounts within a loan portfolio
# 4) period: a scalar indicating the contractual period of all loans
# 5) n: the number of loans within the portfolio
# 6) i.rate: either a single effective rate (to be repeated across the portfolio) or a vector of effective rates per annum, each element corresponding to the rate of a loan account
# 7) vec.DoD.lambda: a vector of account-level multipliers by which to stress delinquency in g3. This is the lambda-function defined in equation 22
calculate.MDoD <- function(vec.Instal, mat.Receipt, vec.Principal, period, n, i.rate, vec.DoD.lambda) {
  
  # - if a vector of interest rates are given, then use that, otherwise repeat the given interest rate across the portfolio
  if (NROW(i.rate) == 1) {
    vec.i <- rep(i.rate, n)
  } else {
    vec.i <- i.rate
  }
  
  # - convert given effective rates per annum to nominal rate per period (convertibly monthly per annum, in this case)
  i_p.rate <- ((1+vec.i)^(1/12) - 1)*12
  
  # - transform the given instalment vector to a matrix by repeating the level instalment across the history of each loan account, which 
  # is mapped to each column in the matrix (loan term, and therefore number of rows, is the same for each account by design)
  # Also, include time 0 as the first row
  mat.Instal <- sapply(1:n, function(i,I) {
    return(c(0,rep(I[i],period))) # instalment at time t=0 is 0 by design
  },I=vec.Instal);
  
  # - calculate the present value (at time 0) of every instalment of every loan account: row=period, column=account
  # Also, include time 0 as the first row
  mat.InstalPV <- sapply(1:n, function(i,I) {
    return(c(0, (1+i_p.rate[i]/12)^(-1*(1:period)) * I[i]))
  },I=vec.Instal);
  
  # - calculate \delta_t as the difference between I_t and R_t for each loan account at each period t=0,...,T
  mat.Diff <- mat.Instal - rbind(rep(0,n), mat.Receipt);
  
  # - Calculate expected duration: Part 1 (standard Macaulay Duration, but without summing)
  # 1) Weigh each discounted instalment of an account by the loan principal, at each period
  # 2) Multiply this with time elapsed (loan age) per annual period
  MDoD.ExpWeightTime <- (mat.InstalPV %*% diag(1/vec.Principal)) * (0:period / 12);
  
  # - create various matrices to be populated later
  MDoD.ActWeightTime <- MDoD.ExpWeightTime # to be modified recursively later, based on actual experience (receipts)
  MDoD.CashPV <- mat.InstalPV
  vec.term <- rep(period,n) # assume all accounts have the same contractual term
  
  # - create a copy of (undiscounted) instalments, such that the last instalment of each account
  # can be updated with arrears recursively later
  MDoD.CashFlow.Star <- mat.Instal
  
  # - create matrices containing the eventual expected and actual duration values that are recursively calculated later
  MDoD.ExpDuration <- matrix(-1.00, nrow=period+1, ncol=n)
  MDoD.ActDuration <- MDoD.ExpDuration # basically serves line 9's purpose in g3's algorithm
  
  # - create a vew n-sized vectors with one value per account within a portfolio
  vec.Term.Star <- vec.term # a vector of the behavioural term, set to the contractual term, only to be incremented gradually once an account becomes out-of-contract
  vec.InContract <- rep(1,n) # a vector of Boolean flags indicating whether an account is still within its contractual tenure (1) or not (0)
  vec.saved.Arrears <- rep(0.00, n) # a vector of accumulated arrears to be added to the last contractually expected instalment
  
  # - main loop: iterate by period, progressing through the history of a loan
  for (i in 1:(period+1)) {
    
    # -- g3-specific
    # - retrieve the last instalment (previously modified or not) of each account for modification (if necessary)
    # This is line 4 in g3's algorithm, calculating \alpha
    vec.saved.Arrears <- MDoD.CashFlow.Star[matrix(c(vec.Term.Star+1,1:n),nrow=n)]
    
    # -- g3-specific
    # - update the behaviour term of each loan (if necessary)
    # If an account is out-of-contract, then increment the behavioural term, otherwise, simply return the contractual term
    # This is line 5 in g3's algorithm, calculating \Tau
    vec.Term.Star <- sapply(1:n, function(k,y) { 
      val <- period
      if (y[k] == 0) val <- i -1
      return(val) },y=vec.InContract)
    
    # -- g3-specific
    # - determine if an account is out-of-contract: if yes, return 0, otherwise return 1
    # The current period i is compared to the contractual term, whilst adjusting for t=0 included in all matrices as row 1
    # This is implementing the decision function \delta_3() used in constructing g3
    vec.InContract <- ifelse(i > (vec.term + 1), 0, 1)
    
    # -- g3-specific
    # - only execute for t>=1 (ignore origination time point)
    if (i > 1) {
      
      # - If in-contract, then add arrears (if any) to last contractually expected instalment per account, after accumulating these arrears for one period with interest.
      # - If out-of-contract, then accumulate the previously-modified last instalment for one period with interest 
      # This is line 6 in g3's algorithm, calculating I'_(\Tau).
      MDoD.CashFlow.Star[matrix(c(vec.Term.Star+1,1:n),nrow=n)] <- mat.Diff[i,] * (1+vec.delta_pp)^(vec.Term.Star-i+1) +
        MDoD.CashFlow.Star[matrix(c(vec.Term.Star+1,1:n),nrow=n)] * vec.InContract + 
        vec.saved.Arrears * (1-vec.InContract) * (1+vec.delta_pp)
      
    }
    
    # - iterate for each account (across columns), at this particular time period i
    for (j in 1:n) {  
      
      # - only execute for t>=1 (ignore origination time point)
      if (i>1) {
        
        # - calculate remaining discounting periods for jth loan as at time i, based on in-contract Boolean flags, whilst adjusting for t=0 included in all matrices as row 1
        # This is line 7 in g3's algorithm, calculating \beta(m)
        periods <- (0:(vec.Term.Star[j]-i+1)) + (1-vec.InContract[j])
        
        # - calculate discount factors for jth loan, corresponding to the remaining discounting periods, as at time i
        discount <- (1+i_p.rate[j]/12)^(-periods)
        
        # - Recalculate the expected Macaulay Duration as at time i for the jth's loan remaining expected (and unmodified) instalments
        # This is line 8 in g3's algorithm, calculating f_{ED}(t=i) - before summation
        MDoD.ExpWeightTime[i:(vec.Term.Star[j]+1),j] <- mat.Instal[i:(vec.Term.Star[j]+1),j] * discount / vec.Principal[j] * (periods)/12
        
        # - Recalculate the present value of each element within the jth loan's remaining and future cash flows, as at time i
        # Note that this may include the arrears-modified last instalment (i.e., line 6's effect in g3's algorithm)
        MDoD.CashPV[i:(vec.Term.Star[j]+1),j] <- discount * MDoD.CashFlow.Star[i:(vec.Term.Star[j]+1),j]
        
        # - Recalcualte the actual Macaualy Duration as at time i for the jth's loan remaining future cash flows
        # This is line 10 in g3's algorithm, calculating f_{AD}(t=i) - before summation
        MDoD.ActWeightTime[i:(vec.Term.Star[j]+1),j] <- MDoD.CashPV[i:(vec.Term.Star[j]+1),j] / vec.Principal[j] * (periods)/12
      }
      
      # -- Calculate the Expected/Actual Duration vectors for kth loan
      # - sum across remaining periods as at time i, for final step in calculating Macaulay Duration
      MDoD.ExpDuration[i,j] <- sum(MDoD.ExpWeightTime[i:(vec.Term.Star[j]+1),j])
      MDoD.ActDuration[i,j] <- sum(MDoD.ActWeightTime[i:(vec.Term.Star[j]+1),j])
      
    } # exit inner loop
  } # exit main loop
  
  # - create a final matrix containing measurements from the g2 function: (row=period, column=loan)
  mat.MD <- MDoD.ActDuration / MDoD.ExpDuration
  
  # - create Boolean-valued decision matrix, filled with values from the decision function d4
  MDoD.ExceedsExp <- ifelse(MDoD.ActDuration > MDoD.ExpDuration, 1, 0)
  
  # - create a final matrix containing measurements from the g3 function: (row-period, column=loan)
  mat.DoD <- MDoD.ActDuration / MDoD.ExpDuration * (MDoD.ExceedsExp %*% diag(vec.DoD.lambda) + 1)
  
  # - Failsafe: eliminate negative values in the MD-matrix and DoD-matrix
  # This may happen due to bad data input (e.g., negative instalments)
  mat.MD <- ifelse(mat.MD >= 0, 1, 0)*mat.MD
  mat.DoD <- ifelse(mat.DoD >= 0, 1, 0)*mat.DoD
  
  # - Failsafe: eliminate NaNs in the MD-matrix and DoD-matrix by replacing them with NAs (missing values)
  mat.MD[is.na(mat.MD)] <- NA
  mat.DoD[is.na(mat.DoD)] <- NA
  
  return(list(MD=mat.MD,DoD=mat.DoD, ActDur=MDoD.ActDuration, ExpDur=MDoD.ExpDuration))
}




# -- Assuming variable instalments and interest rates (nominal)
# ==== Function: Calculates g2: Macaulay Duration Index (MD-measure), g3: Degree of Delinquency (DoD-measure)
# - Inputs the following:
# 1) vec.Instal: a matrix of time-indexed instalments (starting at t=1), row=period and column=loan
# 2) mat.Receipt: a matrix of time-indexed cashflows/receipts (starting at t=1), row=period and column=loan
# 3) vec.LoanAmt: a vector of loan amounts/principals, constituting individual accounts within a loan portfolio
# 4) vec.Mat: a vector of loan terms
# 5) n: the number of loans within the portfolio
# 6) mat.IntRates.g: a matrix of time-indexed nominal interest rates (starting at t=1), row=period and column=loan
# 7) vec.DoD.lambda: a vector of account-level multipliers by which to stress delinquency in g3. This is the lambda-function defined in equation 22 of paper1
# 8) silent: a flag to indicate silent (T) or verbose (F) operation
calculate.MDoD.forData <- function(mat.Instal.g, mat.Receipt.g, vec.LoanAmt, vec.Mat, n, mat.IntRates.g, vec.DoD.lambda, silent=F ) {
  
  # - for testing purposes
  # mat.Instal.g <- mat.Instal.Treated
  # mat.Receipt.g <- mat.ReceiptAlt.Treated7l
  # vec.Mat <- vec.Term.Treated
  # mat.IntRates.g <- mat.IntRates.Treated
  
  # - append the instalment and receipt matrices with one zero-valued row at the start, to indicate origination (t=0)
  mat.Instal.t <- rbind(rep(0, n), mat.Instal.g);
  mat.Receipt.t <- rbind(rep(0, n), mat.Receipt.g);
  
  # - get max maturity/loan tenure/contractual term (depending on what is input)
  period.g <- max(vec.Mat);
  
  # - create a corresponding time-indexed matrix of continuously compounded interest rates per annum, whilst adding an extra zero-valued row at the start
  # to cater for t=0 (origination)
  mat.delta_pp <- rbind(rep(0,n), (exp(log(1 + (1 + mat.IntRates.g/12)^(12)-1)/12) - 1)) #used in MD/DoD delinquency calculation
  
  # - calculate the present value (at time 0) of every instalment of every loan account: row=period, column=account
  # Also, include time 0 as the first row
  mat.InstalPV <- sapply(1:n, function(i,I, int) {
    return(c(0, (1+int[1:period.g,i]/12)^(-1*(1:period.g)) * I[1:period.g,i]))
  },I=mat.Instal.g,int=mat.IntRates.g);
  
  # - cater for missing values
  mat.InstalPV[is.na(mat.InstalPV)] <- 0
  
  # - calculate \delta_t as the difference between I_t and R_t for each loan account at each period t=0,...,T
  mat.Diff <- mat.Instal.t - mat.Receipt.t;
  
  #WeightTime is: weight (PV of instalment at every time point, divided by constant principal) x time (granular compound-on-book)
  
  # - Calculate expected duration: Part 1 (standard Macaulay Duration, but without summing)
  # 1) Weigh each discounted instalment of an account by the loan principal, at each period
  # 2) Multiply this with time elapsed (loan age) per annual period
  MDoD.ExpWeightTime <- sapply(1:n, function(i,I,P) {
      return(I[,i]/P[i] * (0:period.g / 12) )
    },I=mat.InstalPV, P=vec.LoanAmt);
  
  # - create various matrices to be populated later
  MDoD.ActWeightTime <- MDoD.ExpWeightTime # to be modified recursively later, based on actual experience (receipts)
  MDoD.CashPV <- mat.InstalPV
  vec.term <- vec.Mat # as input
  
  # - create a copy of (undiscounted) instalments, such that the last instalment of each account
  # can be updated with arrears recursively later
  MDoD.CashFlow.Star <- mat.Instal.t
  
  # - create matrices containing the eventual expected and actual duration values that are recursively calculated later
  MDoD.ExpDuration <- matrix(-1.00, nrow=period+1, ncol=n)
  MDoD.ActDuration <- MDoD.ExpDuration # basically serves line 9's purpose in g3's algorithm
  
  # - create a few n-sized vectors with one value per account within a portfolio
  vec.Term.Star <- vec.term # a vector of the behavioural term, set to the input contractual term, only to be incremented gradually once an account becomes out-of-contract
  vec.InContract <- rep(1,n) # a vector of Boolean flags indicating whether an account is still within its contractual tenure (1) or not (0)
  vec.saved.Arrears <- rep(0.00, n) # a vector of accumulated arrears to be added to the last contractually expected instalment
  
  
  ptm <- proc.time()
  
  # - testing purposes
  #period.g <- 165 #error occurs when period.g=166, so run loop for period.g=165 and manually step through subsequent loop until error found.
  # i <- 1
  # i <- 2
  
  # - main loop: iterate by period (row), progressing through the history of each loan (column)
  for (i in 1:(period.g+1)) {
    # i <- 167
    
    # -- g3-specific
    # - retrieve the last instalment (previously modified or not) of each account for modification (if necessary)
    # This is line 4 in g3's algorithm, calculating \alpha    
    vec.saved.Arrears <- MDoD.CashFlow.Star[matrix(c(vec.Term.Star+1,1:n),nrow=n)]
    
    # -- g3-specific
    # - update the 'behavioural' term of each loan (if necessary)
    # If an account is out-of-contract, then increment the behavioural term, otherwise, simply return the contractual term
    # This step will only effectively execute if vec.InContract has been changed for any particular loan in a previous loop-iteration
    # This is line 5 in g3's algorithm, calculating \Tau
    vec.Term.Star <- sapply(1:n, function(k,y,per) { 
      val <- per[k]
      if (y[k] == 0) val <- i -1
      return(val) },y=vec.InContract, per=vec.Mat)
    
    # -- g3-specific
    # - determine if an account is out-of-contract: if yes, return 0, otherwise return 1
    # The current period i is compared to the contractual term, whilst adjusting for t=0 included in all matrices as row 1
    # This is implementing the decision function \delta_3() used in constructing g3
    vec.InContract <- ifelse(i > (vec.term + 1 ) & !is.na(mat.Receipt.t[i,]), 0, 1)
    
    # -- g3-specific
    # - only execute for t>=1 (ignore origination time point)
    if (i > 1) {
      
      # - If in-contract, then add arrears (if any) to last contractually expected instalment per account, after accumulating these arrears by one period's interest.
      # - If out-of-contract, then accumulate the (possibly previously-modified) LAST contractual instalment by one period's interest 
      # This is line 6 in g3's algorithm, calculating I'_(\Tau).
      MDoD.CashFlow.Star[matrix(c(vec.Term.Star+1,1:n),nrow=n)] <- mat.Diff[i,] * (1+mat.delta_pp[i,])^(vec.Term.Star-i+1) +
        MDoD.CashFlow.Star[matrix(c(vec.Term.Star+1,1:n),nrow=n)] * vec.InContract + 
        vec.saved.Arrears * (1-vec.InContract) * (1+mat.delta_pp[i,])
      
    }
    
    # - iterate for each account (across columns), at this particular time period i

    for (j in 1:n) {  
      #j <- 2; i <- 2
      # j <- 1
      
      # Only compute the following if history has been observed for the jth loan.
      # This disallows unobservable time periods for certain loans, which only happens because of loans having different terms
      if (i>1 & !is.na(mat.Receipt.t[i,j])) {
        
        # - calculate remaining discounting periods for jth loan as at time i, based on in-contract Boolean flags, whilst adjusting for t=0 included in all matrices as row 1
        # This is line 7 in g3's algorithm, calculating \beta(m)
        periods <- (0:(vec.Term.Star[j]-i+1)) + (1-vec.InContract[j])
        
        # - calculate discount factors for jth loan, corresponding to the remaining discounting periods, as at time i
        discount <- (1+mat.IntRates.g[i-1,j]/12)^(-periods)
        
        # - Recalculate the expected Macaulay Duration as at time i for the jth's loan remaining expected (and unmodified) instalments
        # This is line 8 in g3's algorithm, calculating f_{ED}(t=i) - before summation
        MDoD.ExpWeightTime[i:(vec.Term.Star[j]+1),j] <- mat.Instal.t[i:(vec.Term.Star[j]+1),j] * discount / vec.LoanAmt[j] * (periods)/12
        
        # - Recalculate the present value of each element within the jth loan's remaining and future cash flows, as at time i
        # Note that this may include the arrears-modified last instalment (i.e., line 6's effect in g3's algorithm)
        MDoD.CashPV[i:(vec.Term.Star[j]+1),j] <- discount * MDoD.CashFlow.Star[i:(vec.Term.Star[j]+1),j]
        
        # - Recalculate the actual Macaulay Duration as at time i for the jth's loan remaining future cash flows
        # This is line 10 in g3's algorithm, calculating f_{AD}(t=i) - before summation
        MDoD.ActWeightTime[i:(vec.Term.Star[j]+1),j] <- MDoD.CashPV[i:(vec.Term.Star[j]+1),j] / vec.LoanAmt[j] * (periods/12)
        
      }
      
      # -- Calculate the Expected/Actual Duration vectors for kth loan
      # only execute for observable histories
      if (!is.na(mat.Receipt.t[i,j])) {

        # - sum across remaining periods as at time i, for final step in calculating Macaulay Duration
        MDoD.ExpDuration[i,j] <- sum(MDoD.ExpWeightTime[i:(vec.Term.Star[j]+1),j])
        MDoD.ActDuration[i,j] <- sum(MDoD.ActWeightTime[i:(vec.Term.Star[j]+1),j])
        
      } else {
        MDoD.ExpDuration[i,j] <- NA
        MDoD.ActDuration[i,j] <- NA
      }
    } # exit inner loop
    
  } # exit main loop
  
  if (silent==F) {
    print(proc.time() - ptm) #IGNORE: computation time taken
  }
  
  # - create a final matrix containing measurements from the g2 function: (row=period, column=loan)
  mat.MD <- MDoD.ActDuration / MDoD.ExpDuration
  
  # - create Boolean-valued decision matrix, filled with values from the decision function d4
  MDoD.ExceedsExp <- ifelse(MDoD.ActDuration > MDoD.ExpDuration, 1, 0)
  MDoD.ExceedsExp[is.na(MDoD.ExceedsExp)] <- 0 # cater for missing values
  
  # - create a final matrix containing measurements from the g3 function: (row-period, column=loan)
  mat.DoD <- sapply(1:n, function(i,Act,Ex,Exceeds,l) {
    # i<-4
    # Act <- MDoD.ActDuration; Ex <- MDoD.ExpDuration; Exceeds<-MDoD.ExceedsExp; l<-vec.DoD.lambda
    return (Act[,i] / Ex[,i] * (Exceeds[,i] * l[i] + 1))
  },Act=MDoD.ActDuration, Ex=MDoD.ExpDuration, Exceeds=MDoD.ExceedsExp, l=vec.DoD.lambda);
  
  
  # - Failsafe: eliminate NaNs in the MD-matrix and DoD-matrix by replacing them with NAs (missing values)
  mat.MD[is.na(mat.MD)] <- 0 #used to be NA
  mat.DoD[is.na(mat.DoD)] <- 0 #used to be NA
  
  # tests for infinite values
  # vec.Test <- sapply(1:n, function(i) {
  #   return(Position(function(x) is.infinite(x)==T, mat.MD[,i], nomatch=-1))
  # })
  # describe(vec.Test)
  # 2 infinite values found
  # POST HOC: no inf values found.
  # find.accs <- which( vec.Test > -1)
  # mat.MD[, find.accs[1]]
  # Infinite values are always preceded by one NA and sometimes succeeded by one NA.
  # Always occur at time 184/183.
  # vec.Maturity[find.accs]
  # no comon pattern in maturities, except that they are all early-ish (max of 13)
  # vec.term[find.accs]
  # all of them have a term of 182 (first 4) and 181 (next 3). Seems connected with positioning of infinite values
  
  # tests for NaN/NA values
  # vec.Test <- sapply(1:n, function(i) {
  #   return(Position(function(x) is.na(x)==T, mat.MD[,i], nomatch=-1))
  # })
  # describe(vec.Test)
  # NA's seem ubiquiteous
  # POST HOC: NaNs replaced with NA, by design.
  # find.accs <- which( vec.Test > -1)
  # mat.MD[, find.accs[1]]
  
  # - test specific loanID: 180
  # testID <- 4
  # mat.Bal[,testID]
  # mat.Receipt.g[, testID]
  # mat.Instal.g[,testID]
  # mat.MD[,testID]
  
  return(list(MD=mat.MD,DoD=mat.DoD, ActDur=MDoD.ActDuration, ExpDur=MDoD.ExpDuration))
}

