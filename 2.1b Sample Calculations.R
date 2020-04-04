# =================== Sample Calculations
# This script creates a set of scalars, vectors, and matrices from a given sample, for easier calculation later on


# ====== Set sample pointer
dat.use.sample <- copy(dat.sample6)


# --- Scalars
n <- NROW(unique(dat.use.sample[,Account]))
period <- max(dat.use.sample[,Maturity])

# -- Parameters for Delinquency Measures
sc.Thres <- 0.9; # repayment ratio: g_1-measure
sc.DelinqSens <- 1; # Delinquency Sensitivity: g_3-measure
sc.maxLoan <- max(dat.use.sample[,LoanAmt]) # Delinquency Sensitivity: g_3-measure


# --- Vectors

# Maturity vector (number of monthly observations per account)
vec.Maturity <- as.vector(unique(dat.use.sample[order(LoanID),list(Account, Maturity)])$Maturity)

# Loan amount vector
vec.LoanAmt <- as.vector(unique(dat.use.sample[order(LoanID),list(Account, LoanAmt)])$LoanAmt)
# Term vector (use a slightly different aggregation method since account 60016933001 had a term change)
vec.Term <- as.vector(unique(dat.use.sample[order(LoanID),list(TermOriginal = TermOriginal[1]),by=list(Account)])$TermOriginal)
#   describe(vec.Term)
#   results: mean of 238, median of 240. Lowest of 180, highest of 240. Don't worry too much about this.

# Closed Vector
vec.Closed <- as.vector(unique(dat.use.sample[order(LoanID),list(Account, Closed.Ind)])$Closed.Ind)

# Write-off Vector
vec.Woff <- as.vector(unique(dat.use.sample[order(LoanID),list(Account, Write.Off.Ind)])$Write.Off.Ind)

# Last date observed vector
vec.LastObs <-  as.Date(paste0(dat.use.sample[,list(LastDate = max(Datex)),by=list(LoanID)]$LastDate,"01"), "%Y%m%d")

# Number of periods from last date observed until end of observation period
vec.RemToPerfWindow <- round(as.double(difftime(
  as.Date(paste0(date.end,"01"), "%Y%m%d"), vec.LastObs, units = "days")) / 365.25 * 12)

# Vector to indicate active accounts (excludes settled, closed, written-off)
vec.Consider <- rep(1, n)
vec.Consider[which(vec.Closed==1)] <- 0
#vec.Consider[which(vec.Woff==1)] <- 0

#describe(vec.Consider)
# results: abut 79% of accounts in sample6 are considered. The remainder have been closed.
# This only impacts the vec.Term.Treated vector accordingly, which affects the extent of forecasting in the various receipt matrices in future.
# For example, loan 1 (in sample 6) was closed. It will therefore retain its behavioural term, instead of switching 
# to the contractual term that enables forecasting forwards.
# It makes no sense in our study's context to forecast loans for closed accounts (closed good without accrued delinquency), when
# we want to optimise the recovery decision' timing on a delinquency basis.
# Note that vec.Consider is repeatedly reset during the loss optimisation process itself (after having forecast using its setting here)
# as an artificial proxy for having different loan portfolios with different risk profiles.


# To enable g_3-measure calculation (Delinquency Sensitivity) later on
vec.DoD.lambda <- sc.DelinqSens * (1-((sc.maxLoan-vec.LoanAmt)/sc.maxLoan));


# --- Matrices: rows are periods (up to 240 by sample design), columns are loan accounts
# Receipt matrix
mat.ReceiptAlt <- as.matrix(spread(dat.use.sample[,list(Account,LoanAge,Receipt_Alt)], key=Account, value = Receipt_Alt)[,-1])

# Instalment matrix
mat.Instal <- as.matrix(spread(dat.use.sample[,list(Account,LoanAge,Install_Total)], key=Account, value = Install_Total)[,-1])

# Interest Rate Matrix
mat.IntRates <- as.matrix(spread(dat.use.sample[,list(Account,LoanAge,IntRate)], key=Account, value = IntRate)[,-1])
mat.IntRates <- mat.IntRates / 100;

# Balance Matrix
mat.Bal <- as.matrix(spread(dat.use.sample[,list(Account,LoanAge, Balance)], key=Account, value = Balance)[,-1])

# Transaction-inferred instalment components
mat.Comp.Ins <- as.matrix(spread(dat.use.sample[,list(Account,LoanAge, Insurance)], key=Account, value = Insurance)[,-1])
mat.Comp.InsPayout <- as.matrix(spread(dat.use.sample[,list(Account,LoanAge, Receipt_InsPayout)], key=Account, value = Receipt_InsPayout)[,-1])
mat.Comp.Fees <- as.matrix(spread(dat.use.sample[,list(Account,LoanAge, Fees)], key=Account, value = Fees)[,-1])
mat.Comp.Int <- as.matrix(spread(dat.use.sample[,list(Account,LoanAge, Interest)], key=Account, value = Interest)[,-1])
mat.Comp.WOff <- as.matrix(spread(dat.use.sample[,list(Account,LoanAge, WOff)], key=Account, value = WOff)[,-1])



# ===================  Calculate Delinquency Measures for available loan performance history

# -- Calculate CD (g_1: Contractual Delinquency)
mat.CDAlt <- calculate.CD.forData(mat.Instal, mat.ReceiptAlt, sc.Thres, period, n, method="base")
# -- Calculate MD/DoD (g_2/g_3: Macaulay Duration Index (MD) Measure | Degree of Delinquency (DoD) Measure)
calc.results <- calculate.MDoD.forData(mat.Instal, mat.ReceiptAlt, vec.LoanAmt, vec.Maturity, n, mat.IntRates, vec.DoD.lambda)
mat.MDAlt <- calc.results$MD
mat.DoDAlt <- calc.results$DoD
rm(calc.results) # an optimization, reduces memory usage

# ======== Merge delinquency measurements back to sample

# - CDAlt (using alternative receipts)
dat.temp <- gather(as.data.frame(mat.CDAlt), key=LoanID, value=CDAlt) %>% na.omit() %>% mutate(LoanID = as.numeric(substr(LoanID, 2, length(LoanID)))) %>% 
  group_by(LoanID) %>% mutate(LoanAge = (1:n())-1); setDT(dat.temp, key=c("LoanID","LoanAge"))
dat.use.sample <- merge(dat.use.sample, dat.temp, by=c("LoanID","LoanAge"), all.x=T) #left join (throw away delinquencies at time t=0 for now)

# - cleanup
rm(dat.temp)


# ======== Light Feature Engineering
dat.use.sample[,Max_CDAlt := max(CDAlt),by=list(Account)]
dat.use.sample[, Last_CDAlt := .SD[.N,CDAlt], by=list(Account)] #.SD is just a subset itself of an account's history, here
# we return the last record (referenced via '.N', which refers to the number of rows for each Account) within the [CDAlt] field.




