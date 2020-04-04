# =================== Random Defaults Technique: Experimental Forecasting
# Implement forecasting in "completing the portfolio" across various scenarios

# Dependencies: 2.1b, 2.1c, 2.2a, 2.2b




# ======== Draw truncation parameter repeatedly from fitted distributions

# We have 3 estimates for b (probability of repayment), estimated from 3 different samples:
prob.b.full # 1) Full sample (S1) 
prob.b.del  # 2) Delinquents-only (include write-offs) (S2)
prob.b.woff # 3) Write-offs-only (S3)

# for Weibull curve fitted to write-offs only (S3) # prep.arg1
trunc.k <- do.call(rweibull, c(list(n=n), prep.arg1))
#plot(density(trunc.k)) # like expected shape
#plot(density(trunc.k[which(vec.Del.Ind==1)])) # like expected shape
#plot(density(trunc.k[which(vec.Woff==1)])) # like expected shape

# for exponential curve fitted to delinquents-only (S2) # s1prep.arg3
s1trunc.k3 <- do.call(rexp, c(list(n=n), s1prep.arg3))
#plot(density(s1trunc.k3[which(vec.Del.Ind==1)]), xlim=c(0,50)) # like expected shape




# ======== Forecast receipts experimentally using random defaults model across various scenarios

# ------ Index of experiments:
#### RANDOM DEFAULTS TECHNIQUE
####### Using Weibull distribution on S3 (write-offs)
#mat.Receipt.Use <- mat.ReceiptAlt.Treated7d # for v2_5d(i-vi) (repayment probability estimated from full sample; truncation parameter k drawn randomly from Weibull distribution fitted on Max_CDAlt from write-offs-only sample)
#mat.Receipt.Use <- mat.ReceiptAlt.Treated7e # for v2_5e(i-ix) (repayment probability estimated from delinquents-only sample; truncation parameter k drawn randomly from Weibull distribution fitted on Max_CDAlt from write-offs-only sample)
#mat.Receipt.Use <- mat.ReceiptAlt.Treated7f # for v2_5f(i-ix) (repayment probability estimated from write-offs sample; with truncation parameter k drawn randomly from Weibull distribution fitted on Max_CDAlt from write-offs-only sample)
####### Using Exp distribution on S2 (delinquents)
#mat.Receipt.Use <- mat.ReceiptAlt.Treated7j # for v2_5j(i-iii) (repayment probability estimated from full sample; truncation parameter k drawn randomly from Exponential distribution fitted on Max_CDAlt from delinquents-only sample)
#mat.Receipt.Use <- mat.ReceiptAlt.Treated7k # for v2_5k(i-iii) (repayment probability estimated from delinquents-only sample; truncation parameter k drawn randomly from Exponential distribution fitted on Max_CDAlt from delinquents-only sample)
#mat.Receipt.Use <- mat.ReceiptAlt.Treated7l # for v2_5l(i-iii) (repayment probability estimated from write-offs sample; truncation parameter k drawn randomly from Exponential distribution fitted on Max_CDAlt from delinquents-only sample)

# Treat Receipts matrix to fill out remaining elements calculated instalment, treated with random defaults technique
# using probability of repayment: full sample
# Randomized truncation drawn from Weibull fitted on Max_CDAlt from write-offs-only sample
mat.ReceiptAlt.Treated7d <- matrix(unlist(sapply(1:n, function(i,t,t.max,r,ins,prob=0.95, del, trunck) {
  #i<-1
  rem.receipts <- sample(x=c(0,1), size=(t.max[i]-t[i]), prob=c(1-prob,prob), replace=T) * ins[i]
  
  # - treat simulated receipts for truncation, but only if we should forecast
  if (length(rem.receipts) > 0) {
    # obtain quick-and-dirty CD count
    CD.vec <- 0:(length(which(rem.receipts == 0))-1) + del[t[i], i]
    if (length(which(CD.vec >= trunck[i])) > 0) {
      start.truncation <- which(CD.vec >= trunck[i])[1]
      # now truncate with zeros
      rem.receipts[start.truncation:length(rem.receipts)] <- 0
    }
  }
  
  prep <- c(as.vector(r[1:t[i], i]),
            as.vector(rem.receipts),
            as.vector(rep(NA, period.term - t.max[i]))
  )
  # sanity check: structure of instalment matrix
  if (NROW(prep) != period.term) cat("\nERROR: ", i)
  return(prep)
}, t=vec.Maturity, t.max=vec.Term.Treated, r=mat.ReceiptAlt, ins=mat.Instal.Treated[ cbind( pmin(vec.Maturity+1,period.term) ,1:n)],
prob=prob.b.full, del=mat.CDAlt.obsrvd, trunck=trunc.k)), 
nrow=period.term, byrow=F)


# Treat Receipts matrix to fill out remaining elements calculated instalment, treated with random defaults technique
# using probability of repayment: defaults-only sample
# Randomized truncation drawn from Weibull fitted on Max_CDAlt from write-offs-only sample
mat.ReceiptAlt.Treated7e <- matrix(unlist(sapply(1:n, function(i,t,t.max,r,ins,prob=0.95, del, trunck) {
  #i<-1
  rem.receipts <- sample(x=c(0,1), size=(t.max[i]-t[i]), prob=c(1-prob,prob), replace=T) * ins[i]
  
  # - treat simulated receipts for truncation, but only if we should forecast
  if (length(rem.receipts) > 0) {
    # obtain quick-and-dirty CD count
    CD.vec <- 0:(length(which(rem.receipts == 0))-1) + del[t[i], i]
    if (length(which(CD.vec >= trunck[i])) > 0) {
      start.truncation <- which(CD.vec >= trunck[i])[1]
      # now truncate with zeros
      rem.receipts[start.truncation:length(rem.receipts)] <- 0
    }
  }
  prep <- c(as.vector(r[1:t[i], i]),
            as.vector(rem.receipts),
            as.vector(rep(NA, period.term - t.max[i]))
  )
  # sanity check: structure of instalment matrix
  if (NROW(prep) != period.term) cat("\nERROR: ", i)
  return(prep)
}, t=vec.Maturity, t.max=vec.Term.Treated, r=mat.ReceiptAlt, ins=mat.Instal.Treated[ cbind( pmin(vec.Maturity+1,period.term) ,1:n)],
prob=prob.b.del, del=mat.CDAlt.obsrvd, trunck=trunc.k)), 
nrow=period.term, byrow=F)


# Treat Receipts matrix to fill out remaining elements calculated instalment, treated with random defaults technique
# using probability of repayment: write-offs-only sample
# Randomized truncation drawn from Weibull fitted on Max_CDAlt from write-offs-only sample
mat.ReceiptAlt.Treated7f <- matrix(unlist(sapply(1:n, function(i,t,t.max,r,ins,prob=0.95, del, trunck) {
  #i<-1
  rem.receipts <- sample(x=c(0,1), size=(t.max[i]-t[i]), prob=c(1-prob,prob), replace=T) * ins[i]
  
  # - treat simulated receipts for truncation, but only if we should forecast
  if (length(rem.receipts) > 0) {
    # obtain quick-and-dirty CD count
    CD.vec <- 0:(length(which(rem.receipts == 0))-1) + del[t[i], i]
    if (length(which(CD.vec >= trunck[i])) > 0) {
      start.truncation <- which(CD.vec >= trunck[i])[1]
      # now truncate with zeros
      rem.receipts[start.truncation:length(rem.receipts)] <- 0
    }
  }
  
  prep <- c(as.vector(r[1:t[i], i]),
            as.vector(rem.receipts),
            as.vector(rep(NA, period.term - t.max[i]))
  )
  # sanity check: structure of instalment matrix
  if (NROW(prep) != period.term) cat("\nERROR: ", i)
  return(prep)
}, t=vec.Maturity, t.max=vec.Term.Treated, r=mat.ReceiptAlt, ins=mat.Instal.Treated[ cbind( pmin(vec.Maturity+1,period.term) ,1:n)],
prob=prob.b.woff, del=mat.CDAlt.obsrvd, trunck=trunc.k)), 
nrow=period.term, byrow=F)


# Treat Receipts matrix to fill out remaining elements calculated instalment, treated with random defaults technique
# using probability of repayment: full sample
# Randomized truncation drawn from Exponential distribution fitted on Max_CDAlt from delinquents-only sample
mat.ReceiptAlt.Treated7j <- matrix(unlist(sapply(1:n, function(i,t,t.max,r,ins,prob=0.95, del, trunck) {
  #i<-1
  rem.receipts <- sample(x=c(0,1), size=(t.max[i]-t[i]), prob=c(1-prob,prob), replace=T) * ins[i]
  
  # - treat simulated receipts for truncation, but only if we should forecast
  if (length(rem.receipts) > 0) {
    # obtain quick-and-dirty CD count
    CD.vec <- 0:(length(which(rem.receipts == 0))-1) + del[t[i], i]
    if (length(which(CD.vec >= trunck[i])) > 0) {
      start.truncation <- which(CD.vec >= trunck[i])[1]
      # now truncate with zeros
      rem.receipts[start.truncation:length(rem.receipts)] <- 0
    }
  }
  
  prep <- c(as.vector(r[1:t[i], i]),
            as.vector(rem.receipts),
            as.vector(rep(NA, period.term - t.max[i]))
  )
  # sanity check: structure of instalment matrix
  if (NROW(prep) != period.term) cat("\nERROR: ", i)
  return(prep)
}, t=vec.Maturity, t.max=vec.Term.Treated, r=mat.ReceiptAlt, ins=mat.Instal.Treated[ cbind( pmin(vec.Maturity+1,period.term) ,1:n)],
prob=prob.b.full, del=mat.CDAlt.obsrvd, trunck=s1trunc.k3)), 
nrow=period.term, byrow=F)


# Treat Receipts matrix to fill out remaining elements calculated instalment, treated with random defaults technique
# using probability of repayment: delinquents-only sample
# Randomized truncation drawn from Exponential distribution fitted on Max_CDAlt from delinquents-only sample
mat.ReceiptAlt.Treated7k <- matrix(unlist(sapply(1:n, function(i,t,t.max,r,ins,prob=0.95, del, trunck) {
  #i<-1
  rem.receipts <- sample(x=c(0,1), size=(t.max[i]-t[i]), prob=c(1-prob,prob), replace=T) * ins[i]
  
  # - treat simulated receipts for truncation, but only if we should forecast
  if (length(rem.receipts) > 0) {
    # obtain quick-and-dirty CD count
    CD.vec <- 0:(length(which(rem.receipts == 0))-1) + del[t[i], i]
    if (length(which(CD.vec >= trunck[i])) > 0) {
      start.truncation <- which(CD.vec >= trunck[i])[1]
      # now truncate with zeros
      rem.receipts[start.truncation:length(rem.receipts)] <- 0
    }
  }
  
  prep <- c(as.vector(r[1:t[i], i]),
            as.vector(rem.receipts),
            as.vector(rep(NA, period.term - t.max[i]))
  )
  # sanity check: structure of instalment matrix
  if (NROW(prep) != period.term) cat("\nERROR: ", i)
  return(prep)
}, t=vec.Maturity, t.max=vec.Term.Treated, r=mat.ReceiptAlt, ins=mat.Instal.Treated[ cbind( pmin(vec.Maturity+1,period.term) ,1:n)],
prob=prob.b.del, del=mat.CDAlt.obsrvd, trunck=s1trunc.k3)), 
nrow=period.term, byrow=F)


# Treat Receipts matrix to fill out remaining elements calculated instalment, treated with random defaults technique
# using probability of repayment: write-offs-only sample
# Randomized truncation drawn from Exponential distribution fitted on Max_CDAlt from delinquents-only sample
mat.ReceiptAlt.Treated7l <- matrix(unlist(sapply(1:n, function(i,t,t.max,r,ins,prob=0.95, del, trunck) {
  #i<-1
  rem.receipts <- sample(x=c(0,1), size=(t.max[i]-t[i]), prob=c(1-prob,prob), replace=T) * ins[i]
  
  # - treat simulated receipts for truncation, but only if we should forecast
  if (length(rem.receipts) > 0) {
    # obtain quick-and-dirty CD count
    CD.vec <- 0:(length(which(rem.receipts == 0))-1) + del[t[i], i]
    if (length(which(CD.vec >= trunck[i])) > 0) {
      start.truncation <- which(CD.vec >= trunck[i])[1]
      # now truncate with zeros
      rem.receipts[start.truncation:length(rem.receipts)] <- 0
    }
  }
  
  prep <- c(as.vector(r[1:t[i], i]),
            as.vector(rem.receipts),
            as.vector(rep(NA, period.term - t.max[i]))
  )
  # sanity check: structure of instalment matrix
  if (NROW(prep) != period.term) cat("\nERROR: ", i)
  return(prep)
}, t=vec.Maturity, t.max=vec.Term.Treated, r=mat.ReceiptAlt, ins=mat.Instal.Treated[ cbind( pmin(vec.Maturity+1,period.term) ,1:n)],
prob=prob.b.woff, del=mat.CDAlt.obsrvd, trunck=s1trunc.k3)), 
nrow=period.term, byrow=F)

