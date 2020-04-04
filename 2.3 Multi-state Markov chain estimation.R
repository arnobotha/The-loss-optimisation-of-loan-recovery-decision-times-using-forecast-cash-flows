# =================== Markovian Defaults Technique: Estimation + Experimental Forecasting
# 1) Estimation of the parameters of a more sophisticated Markov chain-based receipt forecasting technique
# 2) Implement forecasting in "completing the portfolio" across various scenarios

# Dependencies: 2.1b, 2.1c



# ======== Observe delinquency states over time for each account

# - prepare base matrix into which history is categorised into the states of the Markov chain
ptm <- proc.time()
mat.DelStates <- matrix(unlist(sapply(1:n, function(i,t,del,woff) {
  prep.val <- case_when(
    del[1:t[i], i] == 0 ~ 1, # "Active/CD0"
    del[1:t[i], i] == 1 ~ 2, # CD1
    del[1:t[i], i] == 2 ~ 3, # CD2
    del[1:t[i], i] == 3 ~ 4, # CD3
    del[1:t[i], i] == 4 ~ 5, # CD4
    del[1:t[i], i] == 5 ~ 6, # CD5
    del[1:t[i], i] >= 6 ~ 7  # CD6+
  )
  # now truncate last element when written-off 
  prep.val[t[i]] <- woff[i] * 8 + (1-woff[i]) * prep.val[t[i]] # write-off is state 8
  
  val <- c(prep.val, 
           as.vector(rep(NA, period.term - t[i]))
  )
  return(val)
}, t=vec.Maturity, del=mat.CDAlt.obsrvd, woff=vec.Woff)),
nrow=period.term, byrow=F)
proc.time() - ptm #IGNORE: computation time taken



# ======== Estimate MLEs of resulting transition matrices


# -- 1) MLE of transition rates of Markov chain: S1 (Full sample)

states <- c("1.UpToDate", "2.CD1", "3.CD2", "4.CD3", "5.CD4", "6.CD5", "7.CD6+", "8.WOff")
p.DelStates <- matrix(nrow=NROW(states), ncol=NROW(states), 0)
colnames(p.DelStates) <- states; rownames(p.DelStates) <- states
for (t in 1:(period.term-1)) {
  for (state.i in 1:NROW(states)) {
    for (state.j in 1:NROW(states)) {
      p.DelStates[state.i, state.j] <- p.DelStates[state.i, state.j] + NROW(which(mat.DelStates[t,] ==state.i & mat.DelStates[t+1,] ==state.j))   
    }
  }
}
for (state.i in 1:NROW(states)) p.DelStates[state.i, ] <- p.DelStates[state.i, ] / sum(p.DelStates[state.i, ])
p.DelStates

# MLE of initial rates
p.DelStates.init <- c(NROW(which(mat.DelStates[1,] == 1)), NROW(which(mat.DelStates[1,] == 2)), NROW(which(mat.DelStates[1,] == 3)), NROW(which(mat.DelStates[1,] == 4)),
                      NROW(which(mat.DelStates[1,] == 5)), NROW(which(mat.DelStates[1,] == 6)), NROW(which(mat.DelStates[1,] == 7)), NROW(which(mat.DelStates[1,] == 8))
                      ) / n
# save rates to file
write.csv(x=matrix(sprintf("%1.2f%%", p.DelStates*100), ncol=NROW(states), byrow=F, dimnames=list(states,states)), "6a-estimates.csv")

# however, initial state probabilities are irrelevant since:
#   1) they are not sensible (very near 100% at state 1 - which is expected)
#   2) will introduce a sudden "break" in the deliquency development of a loan (when completing it from last observed point), which is jarring
# Rather, we will observe the last CD count of each loan and set that as the starting state in treating the remaining unobserved receipts




# -- 2) MLE of transition rates of Markov chain: S2 (delinquents-only sample)

states <- c("1.UpToDate", "2.CD1", "3.CD2", "4.CD3", "5.CD4", "6.CD5", "7.CD6+", "8.WOff")
p.DelStates.del <- matrix(nrow=NROW(states), ncol=NROW(states), 0)
colnames(p.DelStates.del) <- states; rownames(p.DelStates.del) <- states
for (t in 1:(period.term-1)) {
  for (state.i in 1:NROW(states)) {
    for (state.j in 1:NROW(states)) {
      p.DelStates.del[state.i, state.j] <- p.DelStates.del[state.i, state.j] + NROW(which(mat.DelStates[t,which(vec.Del.Ind==1)] ==state.i & mat.DelStates[t+1,which(vec.Del.Ind==1)] ==state.j))   
    }
  }
}
for (state.i in 1:NROW(states)) p.DelStates.del[state.i, ] <- p.DelStates.del[state.i, ] / sum(p.DelStates.del[state.i, ])
p.DelStates.del

# MLE of initial rates
p.DelStates.del.init <- c(NROW(which(mat.DelStates[1,which(vec.Del.Ind==1)] == 1)), NROW(which(mat.DelStates[1,which(vec.Del.Ind==1)] == 2)),
                          NROW(which(mat.DelStates[1,which(vec.Del.Ind==1)] == 3)), NROW(which(mat.DelStates[1,which(vec.Del.Ind==1)] == 4)),
                          NROW(which(mat.DelStates[1,which(vec.Del.Ind==1)] == 5)), NROW(which(mat.DelStates[1,which(vec.Del.Ind==1)] == 6)),
                          NROW(which(mat.DelStates[1,which(vec.Del.Ind==1)] == 7)), NROW(which(mat.DelStates[1,which(vec.Del.Ind==1)] == 8))
) / NROW(cbind((which(mat.DelStates[1,which(vec.Del.Ind==1)] == 1)), (which(mat.DelStates[1,which(vec.Del.Ind==1)] == 2)),
               NROW((mat.DelStates[1,which(vec.Del.Ind==1)] == 3)), (which(mat.DelStates[1,which(vec.Del.Ind==1)] == 4)),
               NROW((mat.DelStates[1,which(vec.Del.Ind==1)] == 5)), (which(mat.DelStates[1,which(vec.Del.Ind==1)] == 6)),
               NROW((mat.DelStates[1,which(vec.Del.Ind==1)] == 7)), (which(mat.DelStates[1,which(vec.Del.Ind==1)] == 8))
))

# save rates to file
write.csv(x=matrix(sprintf("%1.2f%%", p.DelStates.del*100), ncol=NROW(states), byrow=F, dimnames=list(states,states)), "6b-estimates.csv")




# -- 3) MLE of transition rates of Markov chain: S3 (write-offs sample)

states <- c("1.UpToDate", "2.CD1", "3.CD2", "4.CD3", "5.CD4", "6.CD5", "7.CD6+", "8.WOff")
p.DelStates.woff <- matrix(nrow=NROW(states), ncol=NROW(states), 0)
colnames(p.DelStates.woff) <- states; rownames(p.DelStates.woff) <- states
for (t in 1:(period.term-1)) {
  for (state.i in 1:NROW(states)) {
    for (state.j in 1:NROW(states)) {
      p.DelStates.woff[state.i, state.j] <- p.DelStates.woff[state.i, state.j] + NROW(which(mat.DelStates[t,which(vec.Woff==1)] ==state.i & mat.DelStates[t+1,which(vec.Woff==1)] ==state.j))   
    }
  }
}
for (state.i in 1:NROW(states)) p.DelStates.woff[state.i, ] <- p.DelStates.woff[state.i, ] / sum(p.DelStates.woff[state.i, ])
p.DelStates.woff

# MLE of initial rates
p.DelStates.woff.init <- c(NROW(which(mat.DelStates[1,which(vec.Woff==1)] == 1)), NROW(which(mat.DelStates[1,which(vec.Woff==1)] == 2)),
                          NROW(which(mat.DelStates[1,which(vec.Woff==1)] == 3)), NROW(which(mat.DelStates[1,which(vec.Woff==1)] == 4)),
                          NROW(which(mat.DelStates[1,which(vec.Woff==1)] == 5)), NROW(which(mat.DelStates[1,which(vec.Woff==1)] == 6)),
                          NROW(which(mat.DelStates[1,which(vec.Woff==1)] == 7)), NROW(which(mat.DelStates[1,which(vec.Woff==1)] == 8))
) / NROW(cbind((which(mat.DelStates[1,which(vec.Woff==1)] == 1)), (which(mat.DelStates[1,which(vec.Woff==1)] == 2)),
               NROW((mat.DelStates[1,which(vec.Woff==1)] == 3)), (which(mat.DelStates[1,which(vec.Woff==1)] == 4)),
               NROW((mat.DelStates[1,which(vec.Woff==1)] == 5)), (which(mat.DelStates[1,which(vec.Woff==1)] == 6)),
               NROW((mat.DelStates[1,which(vec.Woff==1)] == 7)), (which(mat.DelStates[1,which(vec.Woff==1)] == 8))
))

# save rates to file
write.csv(x=matrix(sprintf("%1.2f%%", p.DelStates.woff*100), ncol=NROW(states), byrow=F, dimnames=list(states,states)), "6c-estimates.csv")








# ======== Forecast receipts experimentally using Markovian defaults model across various scenarios

# ------ Index of experiments:
#### MARKOVIAN DEFAULTS TECHNIQUE 
# mat.ReceiptAlt.Treated8a # for v2_6a(i-iii) (treated with multi-state Markovian defaults technique with parameter estimates from full sample)
# mat.ReceiptAlt.Treated8b # for v2_6b(i-iii) (treated with multi-state Markovian defaults technique with parameter estimates from delinquents-only sample)
# mat.ReceiptAlt.Treated8c # for v2_6c(i-iii) (treated with multi-state Markovian defaults technique with parameter estimates from write-offs-only sample)


# Treat Receipts matrix to fill out remaining elements calculated instalment, treated with markovian defaults technique
# using estimates from data (full sample)
ptm <- proc.time()
mat.ReceiptAlt.Treated8a <- matrix(unlist(sapply(1:n, function(i,t,t.max,r,ins, p.trans, del, n.states) {
  
  # prepare outcome vector to be predicted by markov chain
  pred.outcomes <- rep(0, (t.max[i]-t[i]))
  # prepare receipt vector to be simulated
  rec.receipts <- rep(0, (t.max[i]-t[i])) 
  
  if (t.max[i]-t[i] > 0) {
    
    # sample initial states
    last.observed <- case_when(
      del[i] == 0 ~ 1, # "Active/CD0"
      del[i] == 1 ~ 2, # CD1
      del[i] == 2 ~ 3, # CD2
      del[i] == 3 ~ 4, # CD3
      del[i] == 4 ~ 5, # CD4
      del[i] == 5 ~ 6, # CD5
      del[i] >= 6 ~ 7  # CD6+
    )
    # We purposefully ignore the write-off outcome.
    # We're continuing forecasting as if write-off didn't happen, just to complete the account's remainding receipts
    # This is done just for consistency with the greater experimental scenario setup of our study
    # We don't foresee that a bank will actually do this in practice.
    
    # simulate "initial" outcome based on last observed state
    pred.outcomes[1] <- sample(x=1:n.states, size=1, replace=T, prob=p.trans[last.observed,])
    
    # simulate first receipt based on "initial" outcome
    if (pred.outcomes[1] == 1) {
      rec.receipts[1] <- ins[i]
    } else if (pred.outcomes[1] > 1 & pred.outcomes[1] < n.states) {
      rec.receipts[1] <- 0
    } else {
      rec.receipts[1:(t.max[i]-t[i])] <- 0
    }
    
    # complete the chain
    # and map chain outcomes to payment outcomes in the same loop - unless truncation/write-off is triggered from the "initial" outcome
    if (t.max[i]-t[i] > 1 & pred.outcomes[1] < n.states) {
      
      for (tt in 2:(t.max[i]-t[i]) ) {
        
        pred.outcomes[tt] <- sample(x=1:n.states, size=1, replace=T, prob=p.trans[pred.outcomes[tt-1],])
        
        # cater for truncation/woff (last state, i.e., n.states), if it was sampled above
        if (pred.outcomes[tt] == n.states) {
          
          # fill all remaining outcomes with truncation outcome (zero payments) and exit loop
          pred.outcomes[tt:(t.max[i]-t[i])] <- n.states
          rec.receipts[tt:(t.max[i]-t[i])] <- 0 # no-payments
          break
          
        } else {
          
          trans <- diff(pred.outcomes[(tt-1):tt])
          
          if (trans < 0) { # CD decreased
            
            # therefore, at least 2 full instalmnet was paid
            rec.receipts[tt] <- -(trans - 1) * ins[i]
            
            
          } else if (trans == 0) { # CD stayed the same 
            
            if (pred.outcomes[tt] == (n.states-1)) {
              
              # remember that penultimate state is semi-absorbing "6+"
              # this actually means that if it "stayed the same", we ASSUME
              # that delinquency has actually worsed since CD6 and CD7 will technically be in the same state.
              # therefore, in this case, roll a non-payment
              rec.receipts[tt] <- 0
			  
      			  # NB! The mathematical formulation of this IF..THEN (and indeed, the definition of trans) is different in Paper2 than to this.
      			  # However, it amounts to the very same effect. Paper2 has a \delta_t (which is the difference in g_1), while this code 
      			  # has a [trans], which is the difference in state number. But state is derived from g_1, therefore, the same effect.
			  
            } else {
              
              # therefore, a payment must have been paid since CD stayed the same (but is also not at CD6+)
              rec.receipts[tt] <- -(trans - 1) * ins[i]
              
            }
          } else { # CD increased by 1 (can only ever be one)
            
            # therefore, a zero is effectively simulated
            rec.receipts[tt] <- 0
          }
        }
      }
    }
  }
  
  # stitch together
  prep <- c(as.vector(r[1:t[i], i]),
            as.vector(rec.receipts),
            as.vector(rep(NA, period.term - t.max[i]))
  )
  
  # sanity check: structure of instalment matrix
  if (NROW(prep) != period.term) cat("\nERROR: ", i)
  return(prep)
  
}, t=vec.Maturity, t.max=vec.Term.Treated, r=mat.ReceiptAlt, ins=mat.Instal.Treated[ cbind( pmin(vec.Maturity+1,period.term) ,1:n)],
p.trans=p.DelStates, del=mat.CDAlt.obsrvd[ cbind( pmin(vec.Maturity+1,period.term) ,1:n)  ], n.states=NROW(states))),
nrow=period.term, byrow=F)
proc.time() - ptm #IGNORE: computation time taken




# Treat Receipts matrix to fill out remaining elements calculated instalment, treated with markovian defaults technique
# using estimates from data (delinquents-only sample)
ptm <- proc.time()
mat.ReceiptAlt.Treated8b <- matrix(unlist(sapply(1:n, function(i,t,t.max,r,ins, p.trans, del, n.states) {

  # prepare outcome vector to be predicted by markov chain
  pred.outcomes <- rep(0, (t.max[i]-t[i]))
  # prepare receipt vector to be simulated
  rec.receipts <- rep(0, (t.max[i]-t[i])) 
  
  if (t.max[i]-t[i] > 0) {
    
    # sample initial states
    last.observed <- case_when(
      del[i] == 0 ~ 1, # "Active/CD0"
      del[i] == 1 ~ 2, # CD1
      del[i] == 2 ~ 3, # CD2
      del[i] == 3 ~ 4, # CD3
      del[i] == 4 ~ 5, # CD4
      del[i] == 5 ~ 6, # CD5
      del[i] >= 6 ~ 7  # CD6+
    )
    
    # simulate "initial" outcome based on last observed state
    pred.outcomes[1] <- sample(x=1:n.states, size=1, replace=T, prob=p.trans[last.observed,])
    
    # simulate first receipt based on "initial" outcome
    if (pred.outcomes[1] == 1) {
      rec.receipts[1] <- ins[i]
    } else if (pred.outcomes[1] > 1 & pred.outcomes[1] < n.states) {
      rec.receipts[1] <- 0
    } else {
      rec.receipts[1:(t.max[i]-t[i])] <- 0
    }
    
    # complete the chain
    # and map chain outcomes to payment outcomes in the same loop - unless truncation/write-off is triggered from the "initial" outcome
    if (t.max[i]-t[i] > 1 & pred.outcomes[1] < n.states) {
      
      for (tt in 2:(t.max[i]-t[i]) ) {
        
        pred.outcomes[tt] <- sample(x=1:n.states, size=1, replace=T, prob=p.trans[pred.outcomes[tt-1],])
        
        # cater for truncation/woff (last state, i.e., n.states), if it was sampled above
        if (pred.outcomes[tt] == n.states) {
          
          # fill all remaining outcomes with truncation outcome (zero payments) and exit loop
          pred.outcomes[tt:(t.max[i]-t[i])] <- n.states
          rec.receipts[tt:(t.max[i]-t[i])] <- 0 # no-payments
          break
          
        } else {
          
          trans <- diff(pred.outcomes[(tt-1):tt])
          
          if (trans < 0) { # CD decreased
            
            # therefore, at least 2 full instalmnet was paid
            rec.receipts[tt] <- -(trans - 1) * ins[i]
            
            
          } else if (trans == 0) { # CD stayed the same 
            
            if (pred.outcomes[tt] == (n.states-1)) {
              
              # remember that penultimate state is semi-absorbing "6+"
              # this actually means that if it "stayed the same", we ASSUME
              # that delinquency has actually worsed since CD6 and CD7 will technically be in the same state.
              # therefore, in this case, roll a non-payment
              rec.receipts[tt] <- 0
              
            } else {
              
              # therefore, a payment must have been paid since CD stayed the same (but is also not at CD6+)
              rec.receipts[tt] <- -(trans - 1) * ins[i]
              
            }
          } else { # CD increased by 1 (can only ever be one)
            
            # therefore, a zero is effectively simulated
            rec.receipts[tt] <- 0
            
          }
        }
      }
    }
  }
  
  # stitch together
  prep <- c(as.vector(r[1:t[i], i]),
            as.vector(rec.receipts),
            as.vector(rep(NA, period.term - t.max[i]))
  )
  
  # sanity check: structure of instalment matrix
  if (NROW(prep) != period.term) cat("\nERROR: ", i)
  return(prep)
  
}, t=vec.Maturity, t.max=vec.Term.Treated, r=mat.ReceiptAlt, ins=mat.Instal.Treated[ cbind( pmin(vec.Maturity+1,period.term) ,1:n)],
p.trans=p.DelStates.del, del=mat.CDAlt.obsrvd[ cbind( pmin(vec.Maturity+1,period.term) ,1:n)  ], n.states=NROW(states))),
nrow=period.term, byrow=F)
proc.time() - ptm #IGNORE: computation time taken




# Treat Receipts matrix to fill out remaining elements calculated instalment, treated with markovian defaults technique
# using estimates from data (write-offs-only sample)
ptm <- proc.time()
mat.ReceiptAlt.Treated8c <- matrix(unlist(sapply(1:n, function(i,t,t.max,r,ins, p.trans, del, n.states) {

  # prepare outcome vector to be predicted by markov chain
  pred.outcomes <- rep(0, (t.max[i]-t[i]))
  # prepare receipt vector to be simulated
  rec.receipts <- rep(0, (t.max[i]-t[i])) 
  
  
  if (t.max[i]-t[i] > 0) {
    
    # sample initial states
    last.observed <- case_when(
      del[i] == 0 ~ 1, # "Active/CD0"
      del[i] == 1 ~ 2, # CD1
      del[i] == 2 ~ 3, # CD2
      del[i] == 3 ~ 4, # CD3
      del[i] == 4 ~ 5, # CD4
      del[i] == 5 ~ 6, # CD5
      del[i] >= 6 ~ 7  # CD6+
    )
    
    # simulate "initial" outcome based on last observed state
    pred.outcomes[1] <- sample(x=1:n.states, size=1, replace=T, prob=p.trans[last.observed,])
    
    # simulate first receipt based on "initial" outcome
    if (pred.outcomes[1] == 1) {
      rec.receipts[1] <- ins[i]
    } else if (pred.outcomes[1] > 1 & pred.outcomes[1] < n.states) {
      rec.receipts[1] <- 0
    } else {
      rec.receipts[1:(t.max[i]-t[i])] <- 0
    }
    
    # complete the chain
    # and map chain outcomes to payment outcomes in the same loop - unless truncation/write-off is triggered from the "initial" outcome
    if (t.max[i]-t[i] > 1 & pred.outcomes[1] < n.states) {
      
      for (tt in 2:(t.max[i]-t[i]) ) {
        
        pred.outcomes[tt] <- sample(x=1:n.states, size=1, replace=T, prob=p.trans[pred.outcomes[tt-1],])
        
        # cater for truncation/woff (last state, i.e., n.states), if it was sampled above
        if (pred.outcomes[tt] == n.states) {
          
          # fill all remaining outcomes with truncation outcome (zero payments) and exit loop
          pred.outcomes[tt:(t.max[i]-t[i])] <- n.states
          rec.receipts[tt:(t.max[i]-t[i])] <- 0 # no-payments
          break
          
        } else {

          trans <- diff(pred.outcomes[(tt-1):tt])
          
          if (trans < 0) { # CD decreased
            
            # therefore, at least 2 full instalmnet was paid
            rec.receipts[tt] <- -(trans - 1) * ins[i]
            
          } else if (trans == 0) { # CD stayed the same 
            
            if (pred.outcomes[tt] == (n.states-1)) {
              
              # remember that penultimate state is semi-absorbing "6+"
              # this actually means that if it "stayed the same", we ASSUME
              # that delinquency has actually worsed since CD6 and CD7 will technically be in the same state.
              # therefore, in this case, roll a non-payment
              rec.receipts[tt] <- 0
              
            } else {
              
              # therefore, a payment must have been paid since CD stayed the same (but is also not at CD6+)
              rec.receipts[tt] <- -(trans - 1) * ins[i]
              
            }
          } else { # CD increased by 1 (can only ever be one)
            
            # therefore, a zero is effectively simulated
            rec.receipts[tt] <- 0
            
          }
        }
      }
    }
  }
  
  # stitch together
  prep <- c(as.vector(r[1:t[i], i]),
            as.vector(rec.receipts),
            as.vector(rep(NA, period.term - t.max[i]))
  )
  
  # sanity check: structure of instalment matrix
  if (NROW(prep) != period.term) cat("\nERROR: ", i)
  return(prep)
  
}, t=vec.Maturity, t.max=vec.Term.Treated, r=mat.ReceiptAlt, ins=mat.Instal.Treated[ cbind( pmin(vec.Maturity+1,period.term) ,1:n)],
p.trans=p.DelStates.woff, del=mat.CDAlt.obsrvd[ cbind( pmin(vec.Maturity+1,period.term) ,1:n)  ], n.states=NROW(states))),
nrow=period.term, byrow=F)
proc.time() - ptm #IGNORE: computation time taken

