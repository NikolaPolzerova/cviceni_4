# Hirschberguv algoritmus - vlastni cestou, OPRAVIT
Hirschberg <- function(seq_1, seq_2, match, mismatch, gap, align){
  
  seq_1 <- strsplit(seq_1, "")[[1]]
 
  s1 <- floor(length(seq_1)/2)
  seq_1Half1 <- seq_1[1:s1] # 1. polovina sekvence 1
  seq_1Half2 <- seq_1[(s1+1):length(seq_1)] # 2. polovina sekvence 2
  
  seq_1Half1 <- paste(seq_1Half1, collapse='')
  
  prvniRadek <- NeedlemanWunsch(seq_1Half1, seq_2, match, mismatch, gap)
  
  seq_2 <- strsplit(seq_2, "")[[1]]
  otocena_seq2 <- paste(rev(seq_2), collapse= '') # otoceni druhe sekvence
  otocena_seq2Half <- paste(rev(seq_1Half2), collapse= '') # otoceni druhe poloviny prvni sekvence
  
  druhyRadek <- NeedlemanWunsch(otocena_seq2Half, otocena_seq2, match, mismatch, gap)
  
  otoceny_2radek <- rev(druhyRadek) # otoceni druheho vysledku k souctu
  
  konecnyRadek <- prvniRadek + otoceny_2radek # soucet
  bodDeleni <- which.max(konecnyRadek) # bod deleni pro sekvenci 2
  
  print(bodDeleni) # podminka podle bodu deleni
  if (bodDeleni < 1){
    seq_2Half1 <- seq_2[1:bodDeleni-1]
    seq_2Half2 <- seq_2[bodDeleni:length(seq_2)]
    
    Hirschberg(seq_1Half1, seq_2Half1, match, mismatch, gap)
  }else{
    align <- c(align, seq_1Half1, paste(replicate(nchar(seq_1Half1), '-'), collapse = ''))
  }

  
  return(align)
}