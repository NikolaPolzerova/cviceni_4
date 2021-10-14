Hirschberg <- function(seq_1, seq_2, match, mismatch, gap){
   
  seq_1 <- strsplit(seq_1, "")[[1]]
 
  s1 <- floor(length(seq_1)/2)
  seq_1Half <- seq_1[1:s1] # 1. polovina sekvence 1
  seq_2Half <- seq_1[(s1+1):length(seq_1)] # 2. polovina sekvence 2
  
  seq_1Half <- paste(seq_1Half, collapse='')
  
  prvniRadek <- NeedlemanWunsch(seq_1Half, seq_2, match, mismatch, gap)
  
  seq_2 <- strsplit(seq_2, "")[[1]]
  otocena_seq2 <- paste(rev(seq_2), collapse= '') # otoceni druhe sekvence
  otocena_seq2Half <- paste(rev(seq_2Half), collapse= '') # otoceni druhe poloviny prvni sekvence
  
  druhyRadek <- NeedlemanWunsch(otocena_seq2Half, otocena_seq2, match, mismatch, gap)
  
  otoceny_2radek <- rev(druhyRadek) # otoceni druheho vysledku k souctu
  
  konecnyRadek <- prvniRadek + otoceny_2radek # soucet
  bodDeleni <- which.max(konecnyRadek) # bod deleni pro sekvenci 2
  print(bodDeleni)
  
  return(konecnyRadek)
}