NeedlemanWunsch <- function(seq_1, seq_2, match, mismatch, gap){
  
  seq_1 <- strsplit(seq_1, "")[[1]]
  seq_2 <- strsplit(seq_2, "")[[1]]
  
  m <- 1 + length(seq_1)
  n <- 1 + length(seq_2)
  S <- (0:(n-1))*gap
  
  for (i in 2:m){
    s <- S[1]
    c <- S[1] + gap
    S[1] <- c
    for (j in 2:n){
      if (seq_1[i-1] == seq_2[j-1]){
        a <- match
      }else{
        a <- mismatch
      }
      c <- max(c(S[j]+gap, c+gap, s+a))
      s <- S[j]
      S[j] <- c
    }
  }
  return(S)
}