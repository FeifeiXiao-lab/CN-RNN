CNVcluster_new=function (Y, cp, L, type) {
  bic.v <- vector()
  EM <- vector("list", 3)
  
  st1 = 2
  p1 = rep(1/st1, st1)
  sigma1 = rep(1, st1)
  mu1 = c(-1.8, 0.001) # assume all Deletion and Normal
  priors1 = list(p = p1, mu = mu1, sigma = sigma1)
  EM[[1]] = gausianMixture(Y, cp, priors1, L, st1)
  
  st2 = 2
  p2 = rep(1/st2, st2)
  sigma2 = rep(1, st2)
  mu2 = c(0.001, 1.5) # assume all Duplication and Normal
  priors2 = list(p = p2, mu = mu2, sigma = sigma2)
  EM[[2]] = gausianMixture(Y, cp, priors2, L, st2)
  EM[[2]]$state.new[which(EM[[2]]$state.new == 2)] = 0
  EM[[2]]$state.new[which(EM[[2]]$state.new == 1)] = 2
  EM[[2]]$state.new[which(EM[[2]]$state.new == 0)] = 1
  
  st3 = 3
  mu3 = c(-1.8, 0.001, 1.5)
  p3 = rep(1/st3, st3)
  sigma3 = rep(1, st3) # assume all Deletion, Duplication and Normal
  priors3 <- list(p = p3, mu = mu3, sigma = sigma3)
  EM[[3]] = gausianMixture(Y, cp, priors3, L, st3)
  for (i in 1:3) {
    bic.v[i] <- getOneBIC(Y, EM[[i]]$cp.final)$bic
  }
  if (is.null(type)==T) {
      if (bic.v[3] == min(bic.v)) {
        mins = 3
      } else {
        mins = which.min(bic.v)
      }
  }else{
    if (type=="dup"){
       if (bic.v[2] <= bic.v[3]) {mins=2}
       if (bic.v[2] >  bic.v[3]) {mins=3}
    }
    if (type=="del"){
      if (bic.v[1] <= bic.v[3]) {mins=1}
      if (bic.v[1] >  bic.v[3]) {mins=3}
    }
  }
  EM.f = EM[[mins]]
  newcp = EM.f$cp.final
  h = EM.f$index.final
  cnv.state <- getState(EM = EM.f, mins)
  return(list(newcp = newcp, h = h, cnv.state = cnv.state$cnv.state, 
              cnv.start = cnv.state$cnv.start, cnv.end = cnv.state$cnv.end))
} 

getState <-function (EM = EM, mins) {
  state      = EM$state.new
  cp.f       = EM$cp.final
  start.index  = which(state!=2)  
  cnv.start     = cp.f[start.index]
  cnv.start = cnv.start[!is.na(cnv.start)]
  cnv.end = cp.f[start.index+1]
  cnv.end = cnv.end[!is.na(cnv.end)]
  cnv.state = state[start.index]
  
  if (TRUE %in% (mins == 3)) {
    if(EM$mu.final[1] < 0 & EM$mu.final[3] > 0) {
      cnv.state[which(cnv.state == 1)] = "del"
      cnv.state[which(cnv.state == 3)] = "dup"
    }else if(EM$mu.final[1] < 0 & EM$mu.final[3] < 0){
      cnv.state[which(cnv.state != 2)] = "del"
    }else if(EM$mu.final[1] > 0 & EM$mu.final[3] > 0){
      cnv.state[which(cnv.state != 2)] = "dup"    
    }else if(EM$mu.final[1] > 0 & EM$mu.final[3] < 0){
      cnv.state[which(cnv.state == 1)] = "dup"
      cnv.state[which(cnv.state == 3)] = "del"
    }
  }
  if (TRUE %in% (mins == 1)) {
    cnv.state[which(cnv.state == 1)] = "del"
  }
  else if (TRUE %in%(mins == 2)) {
    cnv.state[which(cnv.state == 1)] = "dup"
  }
  return (list(cnv.state = cnv.state, cnv.start = cnv.start, cnv.end = cnv.end))
}
