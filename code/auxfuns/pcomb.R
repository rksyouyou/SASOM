
tippet = function(x){
    1- (1-min(x))^2
}

fisher = function(x){
    q = -2*(sum(log(x)))
    1- pchisq(q,4)
}


DAPC <- function(x){
  if(any(is.na(x))) return(NA)
  if(max(x) > 0.5) qobs = min(x) else qobs = 2*prod(x)
  if(qobs<=0.5) p = (3/2)*qobs - (1/2)*qobs*log(2*qobs) else p=1-(1-qobs)^2
  return(p)
}
