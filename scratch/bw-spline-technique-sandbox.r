#Make the bw distribution so we can draw from it.
#First, get the AICc-vs-bw results in order
b = bw.lagr[['trace']][order(bw.lagr[['trace']][,1]),c(1,2)]

#Fit a spline through the AICc-vs-bw observations and then use it to smooth AICc across the entire range of the tested bandwidths.
spline = smooth.spline(b)
xxx = seq(b[1,1], tail(b[,1],1), length.out=1001)
smooth = predict(spline, xxx)

#Now restrict our attention to the region of the densest 99% of bandwidth probability mass
maxi = which(cumsum(exp(-0.5*smooth$y))/sum(exp(-0.5*smooth$y))>0.995)[1]
mini = tail(which(cumsum(exp(-0.5*smooth$y))/sum(exp(-0.5*smooth$y))<0.005),1)
xxx = seq(xxx[mini], xxx[maxi], length.out=1001)
smooth = predict(spline, xxx)

#Get the CDF of bandwidth within the region of greatest density
pp = cumsum(exp(-0.5*smooth$y))/sum(exp(-0.5*smooth$y))

#Draw some typical bandwidths from the CDF and produce a model with each.
bws = xxx[sapply(runif(19), function(x) which(x<pp)[1])]






#Make the bw distribution so we can draw from it.
#First, get the AICc-vs-bw results in order
b = bw.lagr2[['trace']][order(bw.lagr2[['trace']][,1]),c(1,2)]

#Fit a spline through the AICc-vs-bw observations and then use it to smooth AICc across the entire range of the tested bandwidths.
spline = smooth.spline(b)
#xxx = seq(b[1,1], tail(b[,1],1), length.out=1001)
xxx = seq(0, 1, length.out=1001)
smooth = predict(spline, xxx)
yy = smooth$y - mean(smooth$y)

#Now restrict our attention to the region of the densest 99% of bandwidth probability mass
maxi = which(cumsum(exp(-0.5*yy))/sum(exp(-0.5*yy))>0.995)[1]
mini = tail(which(cumsum(exp(-0.5*yy))/sum(exp(-0.5*yy))<0.005),1)
xxx = seq(xxx[mini], xxx[maxi], length.out=1001)
smooth = predict(spline, xxx)

#Get the CDF of bandwidth within the region of greatest density
pp = cumsum(exp(-0.5*smooth$y))/sum(exp(-0.5*smooth$y))

#Draw some typical bandwidths from the CDF and produce a model with each.
bws = xxx[sapply(runif(200), function(x) which(x<pp)[1])]









#Make the bw distribution so we can draw from it.
#First, get the AICc-vs-bw results in order
b2 = bw.lagr2[['trace']][order(bw.lagr2[['trace']][,1]),c(1,2)]

#Fit a spline through the AICc-vs-bw observations and then use it to smooth AICc across the entire range of the tested bandwidths.
spline2 = smooth.spline(b2)
#xxx2 = seq(b2[1,1], tail(b2[,1],1), length.out=1001)
xxx2 = seq(0, 1, length.out=1001)
smooth2 = predict(spline2, xxx2)
yy2 = smooth2$y

#Now restrict our attention to the region of the densest 99% of bandwidth probability mass
maxi2 = which(cumsum(exp(-0.5*yy2))/sum(exp(-0.5*yy2))>0.995)[1]
mini2 = tail(which(cumsum(exp(-0.5*yy2))/sum(exp(-0.5*yy2))<0.005),1)
xxx2 = seq(xxx2[mini2], xxx2[maxi2], length.out=1001)
smooth2 = predict(spline2, xxx2)

#Get the CDF of bandwidth within the region of greatest density
pp2 = cumsum(exp(-0.5*smooth2$y))/sum(exp(-0.5*smooth2$y))

#Draw some typical bandwidths from the CDF and produce a model with each.
bws2 = xxx2[sapply(runif(200), function(x) which(x<pp2)[1])]









#Make the bw distribution so we can draw from it.
#First, get the AICc-vs-bw results in order
b = bw.lagr[['trace']][order(bw.lagr[['trace']][,1]),c(1,2)]

#Fit a spline through the AICc-vs-bw observations and then use it to smooth AICc across the entire range of the tested bandwidths.
spline = smooth.spline(b)
#xxx = seq(b[1,1], tail(b[,1],1), length.out=1001)
xxx = seq(0, 1, length.out=1001)
smooth = predict(spline, xxx)
yy = smooth$y - mean(smooth$y)

#Now restrict our attention to the region of the densest 99% of bandwidth probability mass
maxi = which(cumsum(exp(-0.5*yy))/sum(exp(-0.5*yy))>0.995)[1]
mini = tail(which(cumsum(exp(-0.5*yy))/sum(exp(-0.5*yy))<0.005),1)
xxx = seq(xxx[mini], xxx[maxi], length.out=1001)
smooth = predict(spline, xxx)

#Get the CDF of bandwidth within the region of greatest density
pp = cumsum(exp(-0.5*yy))/sum(exp(-0.5*yy))

#Draw some typical bandwidths from the CDF and produce a model with each.
bws = xxx[sapply(runif(19), function(x) which(x<pp)[1])]
