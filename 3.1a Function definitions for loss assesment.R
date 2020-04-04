# =================== Critical Loss assessment function definitions
# This is a necessary script since it contain various background function definitions.


# - function that returns the first starting point of the first default episode within a given loan's history,
# given each loans's term (to restrict search space)
# [MAIN FUNCTION]: Used in 5.1d
default.start.first.v2 <- function(i, thres.d, del.mat, t) {
  # testing purposes
  # i <- 57; del.mat <- mat.CD.Use; thres.d <- d.CD; i <- 129; i <- 7; t.offset<-0
  
  # ------ tests:
  #--times-------------------------------0,1,2, 3,4, 5,6,7,8,9,  10,11,12,13,  14,15,16,17,18,19,20
  # test case (assume CD): test.del <- c(3,3,2, 3,2, 4,5,2,1,2,   3, 4, 3, 2,   4, 5, 6, 5, 4, 3, 2)
  #   with d=3, first default episode ought to be time=0
  #--times---------------------------------0,1,2,3,4,5, 6,7,8
  # test case 2 (assume CD): test.del <- c(0,0,1,2,1,0, 0,1,0)
  #   with d=1, first default episode ought to be time=2
  
  # find positions (times) in delinquency matrix where threshold is/was reached to be (g,d)-defaulting at t
  #   - this should include from time t=0 (at origination)
  vec.found <- which(del.mat[1:(t[i]+1),i] >= thres.d)
  # test case: vec.found <- which(test.del >= 3)
  # test case 2: vec.found <- which(test.del >= 1)
  
  if(length(vec.found) == 1) {
    # only one index found
    episodes.start <- vec.found
  } else {
    # 1. Find positions in these positions where the lagged difference is greater than 1.
    #   - these incidate 'breaks' between episodes 
    # 2. Add 1 to these found positions to move to the 'initial starting points' of the next episode in succession
    # 3. Pre-fix this vector with '1' to re-include the first 'episode' that was deselected previously
    # 4. Given this vector of indices, return starting positions again
    episodes.start <-  vec.found[c(1, which(diff(vec.found) > 1) + 1 )]
  }
  
  # return starting period of first episode (if it exists)
  #   - if it doesn't exist, return -1
  #   - subtract 1 to account for a skewed index that includes in its range t=0 as index=1, 
  #     as a result of using the entire range of {del.mat}
  first.start <- ifelse(length(vec.found) == 0, -1, episodes.start[1] - 1)
  return(first.start)
  #rm(del.mat);
}

