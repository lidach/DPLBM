## put raw lengths into L_a


lengths <- rnorm(40, mL, 5)

ages <- 0:18
L_a <- 64.6*(1-exp(-0.1*(ages - t0)))

lengths2 <- lengths
for (i in 1:length(lengths)){
  for(j in 1:length(L_a)){  
    ifelse(lengths[i] > L_a[j] & lengths[i] < L_a[j+1], lengths2[i] <- L_a[j],0)
  }
  ifelse(lengths[i] > L_a[length(L_a)], lengths2[i] <- L_a[length(L_a)], 0)
  
  xx <- count(lengths2)
}

for(i in 1:length(L_a)){
  if(xx[i,1] != L_a[i]) xx[i+nrow(xx),] <- c(L_a[i],0)
  xx <- na.omit(xx)
}

xx <- tapply(xx$freq, xx$x, FUN = sum)


bheq
