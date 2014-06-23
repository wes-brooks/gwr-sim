interpolate.bw = function(bw.trace, S) {
    #Reorder the points on the trace by bandwidth:
    b = bw.trace[order(bw.trace[,1]),c(1,2)]
    b = b[which(b[,2]!=Inf),]
    
    #Linearly interpolate the bandwidth likelihood trace
    xxx = seq(b[1,1], tail(b[,1],1), length.out=1001)
    smooth = approxfun(b)
    s = smooth(xxx) - mean(smooth(xxx))

    #Now restrict our attention to the region of the densest 99.9% of bandwidth probability mass
    p = cumsum(exp(-s / 2)) / sum(exp(-s / 2))
    maxi = min(which(p > 0.9995)[1], 1001, na.rm=TRUE)
    mini = max(tail(which(p < 0.0005),1), 1, na.rm=TRUE)
    xxx = seq(xxx[mini], xxx[maxi], length.out=1001)
    ss = smooth(xxx) - mean(smooth(xxx))

    #Get the CDF of bandwidth within the region of greatest density
    pp = cumsum(exp(-ss / 2)) / sum(exp(-ss / 2))

    #Draw some typical bandwidths from the CDF and produce a model with each.
    bws = xxx[sapply(runif(S), function(x) which(x<pp)[1])]

    return(bws)
}