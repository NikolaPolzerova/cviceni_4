library(Biostrings)
library(stringi)


NeedlemanWunsch <- function(seq_1, seq_2, match, mismatch, gap){
  
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


Hirschberg_sablona <- function(X,Y, Align, match,mismatch,gap){
  
  Z <- Align[[1]] #inicializace prvniho radku zarovnani, typ DNAStringSet
  W <- Align[[2]] #inicializace druheho radku zarovnani
  
  if (length(X) == 0) #delka X je nula
  {                                                            
   for (i in length(Y)) #podle poctu znaku v Y
    {
      
      Z <- c(Z,DNAString('-')) #prida se mezera                  #ZDE POUŽÍT pøetypování v rámci DNAString .. celé to mùžete beztrestnì spojit pomocí c()
      W <- c(W,Y[i]) #prida se znak z Y
      
    }
    Align <- c(DNAStringSet(Z),DNAStringSet(W))                  #ZDE si znovu pohlídat správný formát použitím DNAString
    print(Align)
  }
  else if (length(Y) == 0) #delka Y je nula
  {                                                              #ZDE for cyklus se píše jinak!!
    for (j in length(X)) #podle poctu znaku v X
    {

      Z <- c(Z,X[j]) #prida se znak z X
      W <- c(W,DNAString('-'))#prida se mezera
      
      Align <- c(DNAStringSet(Z),DNAStringSet(W))
      print(Align)
    }
  }
  else if ((length(X) == 1 & length(Y)) == 1) #delka X a Y je jedna
  {
    Z <- c(Z,X[1])  #prida se znak z X                                                   
    W <- c(W,Y[1]) #prida se znak z Y
    Align <- c(DNAStringSet(Z),DNAStringSet(W))
    print(Align)
  }
  else
  {
      xlen <- length(X) #delka X
      xmid <- floor(xlen/2) #pulka delky X
      ylen <- length(Y) #delka Y
      
      ScoreL <- NeedlemanWunsch(X[1:xmid], Y, match, mismatch, gap) #NW skore pro prvni polovinu X a cele Y
      
      ScoreR <- NeedlemanWunsch(reverse(X[(xmid+1):xlen]),reverse(Y),match,mismatch,gap)    #ZDE se spíše použije funkce reverse() na DNAString type
                                                                                            #NW skore pro druhou polovinu X a cele Y (obe sek jsou obracene)
      
      ymid <- ScoreL + rev(ScoreR)                                      #ZDE BY MÌLO BYT SCORER v reverzní formì.. 
                                                                        #index deleni Y
      ymid <- which.max(ymid)-1                                         #ZDE chybí -1  .. abyste nepøesáhla indexování dále
      
      
      #pro prvni cast
      if (ymid == 0) #index deleni Y je nula  
      ### ZDE JSEM MÍSTO c().. doplnila DNAString
    {
      #volani Hirschberg pro prvni polovinu X a prazdny vektor Y
      Align <- Hirschberg_sablona (X[1:xmid],DNAString(),Align,match,mismatch,gap)
      }
    else
    {
      Align <- Hirschberg_sablona(X[1:xmid], Y[1:ymid], Align, match, mismatch, gap) #volani Hirschberg pro prvni polovinu X a prvni cast Y
    }
    
    #pro druhou cast
    if ((xmid+1)>xlen) #X jiz nelze pulit
    {
      #volani Hirschberg pro prazdny vektor X a druhou cast Y
      Align <-  Hirschberg_sablona (DNAString(),Y[(ymid+1):length(Y)],Align,match,mismatch,gap)
    }
    else if ((ymid+1)>ylen) #Y jiz nejde delit
    {
      #volani Hirschberg pro druhou polovinu X a prazdny vektor Y
      Align <-  Hirschberg_sablona (X[(xmid+1):length(X)],DNAString(),Align,match,mismatch,gap)
      }
    else 
    {
      Align <- Hirschberg_sablona(X[(xmid+1):length(X)], Y[(ymid+1):length(Y)], Align, match, mismatch, gap) #volani Hirschberg pro druhou polovinu X a druhou cast Y
    }
  }
  return(Align)
}



s1<-DNAString('AGTACGCA')
s2<-DNAString('TATGC')
match<-2; mismatch<-1; gap<-2;
Hirschberg_sablona(s1,s2,c(),2,-1,-2)
