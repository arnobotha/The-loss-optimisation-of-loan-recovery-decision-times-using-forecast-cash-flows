# =================== Forecast function definitions
# This is a necessary script since it contain various background function definitions.
# Used primarily for Monte Carlo validation later and for model accuracy tests later


# ---------- Random Defaults forecasting technique

# -- function for forecasting cash flows using Random Defaults, given probability of repayment value [prob.b],
# fitted truncation parameters [fitted.params], and truncation distribution [trunc.dist].
# The rest of the input arguments are simply copies of those objects within the outer data scope passed to within
forecastJob_RandomDefaults <- function(seed.value, n, fitted.params=s1prep.arg3, prob.b=prob.b.full, trunc.dist='exp', vec.Maturity.Use,
                                       vec.Mat.Use, mat.ReceiptAlt, mat.Instal.Use, period.term, mat.CDAlt.obsrvd, sc.Thres) {
  
  # ---- Sample random truncation points
  # -- draw truncation parameters from chosen fitted distribution: one for each loan
  set.seed(seed.value)
  if (trunc.dist == 'exp') {
    # for exponential curve fitted to delinquents-only (S2): used for both S1 and S2
    trunc.k <- do.call(rexp, c(list(n=n), fitted.params))
    #plot(density(s1trunc.k3[which(vec.Del.Ind==1)]), xlim=c(0,50)) # like expected shape    
  } else if (trunc.dist == 'weibull') {
    # for Weibull curve fitted to write-offs only (S3) [USED]
    trunc.k <- do.call(rweibull, c(list(n=n), fitted.params))
  }
  
  # ---- Forecast receipts
  mat.Receipt.output <- matrix(unlist(sapply(1:n, function(i,t,t.max,r,ins,prob=0.95, del, trunck, seedval) {
    #i<-2
    set.seed(seedval + i)
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
    if (NROW(prep) != period.term) cat("\nERROR: ", i, ". Row length: ", NROW(prep))
    return(prep)
  }, t=vec.Maturity.Use, t.max=vec.Mat.Use, r=mat.ReceiptAlt, ins=mat.Instal.Use[ cbind( pmin(vec.Maturity.Use+1,period.term) ,1:n)],
  prob=prob.b, del=mat.CDAlt.obsrvd, trunck=trunc.k, seedval=seed.value)), 
  nrow=period.term, byrow=F)
  
  return (mat.Receipt.output)
}



# ---------- Markovian defaults technique

# -- function for forecasting cash flows using Markovian Defaults, given transition rates [p.trans] and state spaces [states].
# The rest of the input arguments are simply copies of those objects within the outer data scope passed to within
forecastJob_MarkovDefaults <- function(seed.value, p.trans, states, n, vec.Maturity.Use, vec.Mat.Use, mat.ReceiptAlt, mat.Instal.Use, period.term,
                                       mat.CDAlt.obsrvd, sc.Thres) {
  
  mat.Receipt.output.old <- matrix(unlist(sapply(1:n, function(i,t,t.max,r,ins, p.trans, del, n.states, seedval) {
    
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
      set.seed(seedval+i)
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
          
          set.seed(seedval+i+(tt-1))
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
                
                # NB!!!! The mathematical formulation of this IF..THEN (and indeed, the definition of trans is different in Paper2 than to this.
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
    
  }, t=vec.Maturity.Use, t.max=vec.Mat.Use, r=mat.ReceiptAlt, ins=mat.Instal.Use[ cbind( pmin(vec.Maturity.Use+1,period.term) ,1:n)],
  p.trans=p.trans, del=mat.CDAlt.obsrvd[ cbind( pmin(vec.Maturity.Use+1,period.term) ,1:n)  ], n.states=NROW(states), seedval=seed.value)),
  nrow=period.term, byrow=F)
  
  return (mat.Receipt.output.old)
  
}


# -- parallelized function for forecasting cash flows using Markovian Defaults, given transition rates [p.trans] and state spaces [states].
# The rest of the input arguments are simply copies of those objects within the outer data scope passed to within
forecastJob_MarkovDefaults_p <- function(seed.value, p.trans, states, n, vec.Maturity.Use, vec.Mat.Use, mat.ReceiptAlt, 
                                         vec.Instal.last, period.term, vec.Del.last, sc.Thres, createOwnPar=F, cpu.threads=6) {
  
  # parallelized inner function that iterates once for every account
  innerJob <- function(i,t,t.max,r,ins, p.trans, del, n.states, seedval) {
    
    # prepare outcome vector to be predicted by markov chain
    pred.outcomes <- rep(0, (t.max-t))
    # prepare receipt vector to be simulated
    rec.receipts <- copy(pred.outcomes) 
    
    if (t.max-t > 0) {
      
      # sample initial states
      last.observed <- case_when(
        del == 0 ~ 1, # "Active/CD0"
        del == 1 ~ 2, # CD1
        del == 2 ~ 3, # CD2
        del == 3 ~ 4, # CD3
        del == 4 ~ 5, # CD4
        del == 5 ~ 6, # CD5
        del >= 6 ~ 7  # CD6+
      )
      # We purposefully ignore the write-off outcome.
      # We're continuing forecasting as if write-off didn't happen, just to complete the account's remainding receipts
      # This is done just for consistency with the greater experimental scenario setup of our study
      # We don't foresee that a bank will actually do this in practice.
      
      # simulate "initial" outcome based on last observed state
      set.seed(seedval+i)
      pred.outcomes[1] <- sample(x=1:n.states, size=1, replace=T, prob=p.trans[last.observed,])
      
      # simulate first receipt based on "initial" outcome
      if (pred.outcomes[1] == 1) {
        rec.receipts[1] <- ins
      } else if (pred.outcomes[1] > 1 & pred.outcomes[1] < n.states) {
        rec.receipts[1] <- 0
      } else {
        rec.receipts[1:(t.max-t)] <- 0
      }
      
      # complete the chain
      # and map chain outcomes to payment outcomes in the same loop - unless truncation/write-off is triggered from the "initial" outcome
      if (t.max-t > 1 & pred.outcomes[1] < n.states) {
        
        for (tt in 2:(t.max-t) ) {
          
          set.seed(seedval+i+(tt-1))
          pred.outcomes[tt] <- sample(x=1:n.states, size=1, replace=T, prob=p.trans[pred.outcomes[tt-1],])
          
          # cater for truncation/woff (last state, i.e., n.states), if it was sampled above
          if (pred.outcomes[tt] == n.states) {
            
            # fill all remaining outcomes with truncation outcome (zero payments) and exit loop
            pred.outcomes[tt:(t.max-t)] <- n.states
            rec.receipts[tt:(t.max-t)] <- 0 # no-payments
            break
            
          } else {

            trans <- diff(pred.outcomes[(tt-1):tt])
            
            if (trans < 0) { # CD decreased
              
              # therefore, at least 2 full instalmnet was paid
              rec.receipts[tt] <- -(trans - 1) * ins
              
            } else if (trans == 0) { # CD stayed the same 
              
              if (pred.outcomes[tt] == (n.states-1)) {
                
                # remember that penultimate state is semi-absorbing "6+"
                # this actually means that if it "stayed the same", we ASSUME
                # that delinquency has actually worsed since CD6 and CD7 will technically be in the same state.
                # therefore, in this case, roll a non-payment
                rec.receipts[tt] <- 0

              } else {
                
                # therefore, a payment must have been paid since CD stayed the same (but is also not at CD6+)
                rec.receipts[tt] <- -(trans - 1) * ins
                
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
    prep <- c(as.vector(r[1:t]),
              as.vector(rec.receipts),
              as.vector(rep(NA, period.term - t.max))
    )
    
    return(prep)    
  }
  
  
  #ptm <- proc.time() #IGNORE: for computation time calculation
  if (createOwnPar==T) {
    cl.port <- makeCluster(cpu.threads)
    registerDoParallel(cl.port)    
  }
  
  # using foreach() from foreach package for advanced looping (using %do%), including parallelization (using %dopar%)
  mat.Receipt.output <- foreach(it=1:n, .combine='cbind', .verbose=F, .inorder=T, .packages =c('dplyr','data.table'))
    %dopar%
    {
      #it <- 1
      receipt.temp <- innerJob(i=it, t=vec.Maturity.Use[it], t.max=vec.Mat.Use[it], r=mat.ReceiptAlt[1:vec.Maturity.Use[it],it], ins=vec.Instal.last[it],
                               p.trans=p.trans, del=vec.Del.last[it], n.states=NROW(states), seedval=seed.value)
    }  
  
  if (createOwnPar==T) {
    stopCluster(cl.port)
  }
  
  #proc.time() - ptm
  # unparalellized variant took 300s
  # parallelized variant $dopar$ took 111s
  # parallelized variant $dopar$ took 47s
  
  return (mat.Receipt.output)
  
}