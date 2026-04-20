CNVout <-function(lrr,  map, h1 = 5, h2 = 10, h3 = 15,  L = 100, chr, sigma = NULL, type=NULL, precise=10000, FINV=FINV, alpha = 0.01, thre = 10, dis.thre = 5, outname = outname){
  cp.f <- vector("list", dim(lrr)[2]) 
  #par(mfrow=c(4,4))
  for (i in 1:dim(lrr)[2]) {
      #print(paste("remove low coverage",i,sep=""))
      #removeindex=which(lrr[,i]!=Modes(lrr[,i]))
      #Y=lrr[removeindex,i]
      #map1=map[removeindex]
      #hist(Y,main = paste(chr,i,sep = "and"))
      cp <- tryCatch( modifiedSaRa(Y=lrr[,i], alpha = alpha, h1 = h1, h2 = h2, h3 = h3, L= L, FINV=FINV, type=type) ,error=function(e){})
      if(is.null(cp)==TRUE | length(cp$cnv.state)==0){
        next
      } else {
        cnv.state = cp$cnv.state
      cp <- modifiedSaRa(Y=lrr[,i], alpha = alpha, h1 = h1, h2 = h2, h3 = h3, L= L, FINV=FINV,sigma=sigma,type=type)
        cnv.start = cp$cnv.start
        cnv.end   = cp$cnv.end
        cnv.len = cnv.end - cnv.start + 1 	    
        RM.index <- which(cnv.len < thre)
        cnv.len <- cnv.len[-RM.index]
        if (TRUE %in% (length(RM.index) > 0)) {
          cnv.state = cnv.state[-RM.index]
          cnv.start = cnv.start[-RM.index]
          cnv.end   = cnv.end[-RM.index]
        }
        if (length(cnv.end)==0){
          next
        }
        mer <- data.frame(cnv.start = cnv.start, cnv.end = cnv.end, cnv.state = cnv.state, dis.thre =  dis.thre)    
        cnv.state = mer$cnv.state
        cnv.start = mer$cnv.start
        cnv.end   = mer$cnv.end
        cp.f[[i]] = sort(c(cnv.start, cnv.end)) 
        if (TRUE %in% (length(cnv.state)==0)){
          cp.f[[i]] = NULL 
        }else{
          cnv.n  <- length(cnv.state)
          output <- matrix(NA, cnv.n, 7)
          for (j in 1 : cnv.n) {
            output[j,1] <- i  #Save individual id
            index.start <- cnv.start[j]
            index.end   <- cnv.end[j]
            output[j,2] <- chr #chromosome number
            output[j,3] <- map[index.start,2]       #cnv start marker name
            output[j,4] <- map[index.end,3]
            output[j,5] <- (as.numeric(output[j,4]) - as.numeric(output[j, 3]))/1000 #length of CNV
            output[j,6] <- index.end - index.start + 1 #NumSNP
            output[j,7] <- as.character(cnv.state[j]) #copy number state, duplication or deletion
          }
        }}
        write.table(output, paste(outname, ".csv", sep=""), append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
      }
  return (list(cp = cp.f))
}


merge <- function(cnv.start, cnv.end, cnv.state, dis.thre) {
  cnv.len = cnv.end - cnv.start + 1 	
  if (TRUE %in% (length(cnv.start) < 2)) {
    cnv.start <- cnv.start
    cnv.end <- cnv.end
    cnv.state <- cnv.state
  }else {
    for (i in 2:length(cnv.start)) {
      distance <- cnv.start[i]-cnv.end[i-1]-1
      if(TRUE %in% (distance <= dis.thre) & (TRUE %in% (cnv.state[i]==cnv.state[i-1]))) {
        cnv.start <- cnv.start[-i]
        cnv.end <- cnv.end[-(i-1)]
        cnv.state <- cnv.state[-(i-1)]
      }
    }
  }
  return(list(cnv.start = cnv.start, cnv.end =  cnv.end, cnv.state =cnv.state))
}

