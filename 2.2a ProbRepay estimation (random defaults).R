# =================== Random Defaults Technique: Parameter Estimation
# simple probability of repayment estimation (for random defaults technique in treating receipts)

# Dependencies: 2.1b, 2.1c


# ========= Estimate probability of payment, b

# - create indicator matrix, by incorporating overpayments into "P" (Paid) definition
mat.repay.ind2 <- matrix(unlist(sapply(1:n, function(i,t,ins,r) {
  val <- c(ifelse(r[1:t[i],i] >= ins[1:t[i],i], 1, 0),
           as.vector(rep(NA, period.term - t[i]))
  )
  return(val)
}, t=vec.Maturity, ins=mat.Instal, r=mat.ReceiptAlt)),
nrow=period.term, byrow=F)


# -- Estimate probability of payment for next period
# more appropriate method and closer to intended use case in "rolling the dice" for the next receipt in treating the unobserved remainder of each loan account

# MLE: full sample (S1)
prob.b.full <- sum(sapply(1:n, function(i,m,t) { sum(m[1:t[i],i]) }, m=mat.repay.ind2, t=vec.Maturity)) / sum(vec.Maturity)
# 86.7%

# MLE: Delinquents-only sample (S2)
prob.b.del <- sum(sapply(1:length(which(vec.Del.Ind==1)), function(i,m,t,ind) { 
  ii <- which(ind==1)[i]
  return( sum(m[1:t[ii],ii]) )
  }, m=mat.repay.ind2, t=vec.Maturity, ind=vec.Del.Ind)) / sum(vec.Maturity[which(vec.Del.Ind==1)])
# 81%

# MLE: Write-offs sample (S3)
prob.b.woff <- sum(sapply(1:length(which(vec.Woff==1)), function(i,m,t,ind) { 
  ii <- which(ind==1)[i]
  return( sum(m[1:t[ii],ii]) )
}, m=mat.repay.ind2, t=vec.Maturity, ind=vec.Woff)) / sum(vec.Maturity[which(vec.Woff==1)])
# 45%

  