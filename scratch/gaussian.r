write.log = function(message, file, append=TRUE) {
    sink(file, append=append)
    cat(paste(message, "\n", sep=""))
    sink()
}
write.log('entry', 'result.txt', append=FALSE)

require(sp)
require(splancs)
require(foreach)
require(iterators)
require(multicore)
require(geoR)
require(SGL)
require(lagr)
require(doMC)
registerDoMC(1)
write.log('installations complete', 'result.txt')

process = 1
cluster = NA
N = 30 #number of width and length divisions in the domain
settings = 18 #number of unique simulation settings
functions = 3 #number of function types
coord = seq(0, 1, length.out=N) #coordinates of the generated observations

write.log(paste('process:', process, sep=''), 'result.txt')

#Establish the simulation parameters
tau = rep(0.1, settings) #tau is the spatial autocorrelation range parameter for the covariates
rho = rep(c(rep(0,2), rep(0.5,2), rep(0.9,2)), functions) #rho is the correlation of the covariates
sigma.tau = rep(0, settings) #sigma.tau is the spatial autocorrelation range parameter of the noise term
sigma = rep(c(0.5,1), settings/2) #sigma is the variance of the noise term
params = data.frame(tau, rho, sigma.tau, sigma)

#Simulation parameters are based on the value of process
parameters = params[process,]

#Seed the RNG
set.seed(process)
write.log(paste('seed:', process, sep=''), 'result.txt')

#Generate the covariates:
if (parameters[['tau']] > 0) {
    write.log(paste('generating GRFs with tau of ', parameters[['tau']], sep=''), 'result.txt')
    d1 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))$data
    d2 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))$data
    d3 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))$data
    d4 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))$data
    d5 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))$data
} else {
    write.log('generating GRFs with tau of 0', 'result.txt')
    d1 = rnorm(N**2, mean=0, sd=1)
    d2 = rnorm(N**2, mean=0, sd=1)
    d3 = rnorm(N**2, mean=0, sd=1)
    d4 = rnorm(N**2, mean=0, sd=1)
    d5 = rnorm(N**2, mean=0, sd=1)
}
loc.x = rep(coord, times=N)
loc.y = rep(coord, each=N)
write.log('generated GRFs', 'result.txt')

#Use the Cholesky decomposition to correlate the random fields:
S = matrix(parameters[['rho']], 5, 5)
diag(S) = rep(1, 5)
L = chol(S)
write.log('generated correlation matrix', 'result.txt')

#Force correlation on the Gaussian random fields:
D = as.matrix(cbind(d1, d2, d3, d4, d5)) %*% L
X1 = matrix(D[,1], N, N)
X2 = matrix(D[,2], N, N)
X3 = matrix(D[,3], N, N)
X4 = matrix(D[,4], N, N)
X5 = matrix(D[,5], N, N)
write.log('generated Xs', 'result.txt')

#simulate the noise term, either with or without spatial correlation
if (parameters[['sigma.tau']] == 0) {epsilon = rnorm(N**2, mean=0, sd=parameters[['sigma']])}
if (parameters[['sigma.tau']] > 0) {epsilon = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(parameters[['sigma']]**2,parameters[['sigma.tau']]))$data}
write.log('generated epsilon', 'result.txt')

#calculate the B1 coefficient surface for the appropriate function type
if ((process-1) %/% 6 == 0) {
    B1 = matrix(rep(ifelse(coord<=0.4, 0, ifelse(coord<0.6,5*(coord-0.4),1)), N), N, N)
} else if ((process-1) %/% 6 == 1) {
    B1 = matrix(rep(coord, N), N, N)
} else if ((process-1) %/% 6 == 2) {
    Xmat = matrix(rep(rep(coord, times=N), times=N), N**2, N**2)
    Ymat = matrix(rep(rep(coord, each=N), times=N), N**2, N**2)
    D = (Xmat-0.5)**2 + (Ymat-0.5)**2
    d = D[,435]
    B1 = matrix(max(d)-d, N, N)
    B1 = B1 / max(B1)
}
write.log('generated B1 coefficient surface', 'result.txt')

#Generate the response variable and set up the data.frame:
mu = X1*B1
Y = mu + epsilon
sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), X5=as.vector(X5), loc.x, loc.y)
write.log('generated data', 'result.txt')

#Prep oracle:
vars = cbind(B1=as.vector(B1!=0))
oracle = list()
for (i in 1:N**2) { 
    oracle[[i]] = character(0)
    if (vars[i,'B1']) { oracle[[i]] = c(oracle[[i]] , "X1") }
}

#Prep GWR:
allvars = replicate(N**2, c('X1', 'X2', 'X3', 'X4', 'X5'), simplify=FALSE)



#MODELS:
S = 100 #number of draws from the bandwidth distribution for each replication

#LAGR:
write.log('making lagr model.', 'result.txt')
bw.lagr = lagr.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, varselect.method='AICc', kernel=epanechnikov, tol.bw=0.01, bw.type='knn', verbose=TRUE, bwselect.method='AICc', resid.type='pearson')

#Make the bw distribution so we can draw from it.
#First, get the AICc-vs-bw results in order
b = bw.lagr[['trace']][order(bw.lagr[['trace']][,1]),c(1,2)]

#Fit a spline through the AICc-vs-bw observations and then use it to smooth AICc across the entire range of the tested bandwidths.
spline = smooth.spline(b)
xxx = seq(b[1,1], tail(b[,1],1), length.out=1001)
smooth = predict(spline, xxx)
smooth = smooth$y - mean(smooth$y)

#Now restrict our attention to the region of the densest 99% of bandwidth probability mass
maxi = min(1001, which(cumsum(exp(-smooth / 2)) / sum(exp(-smooth / 2)) > 0.995)[1])
mini = max(1, tail(which(cumsum(exp(-smooth / 2))/sum(exp(-smooth / 2)) < 0.005),1))
xxx = seq(xxx[mini], xxx[maxi], length.out=1001)
smooth = predict(spline, xxx)$y

#Get the CDF of bandwidth within the region of greatest density
pp = cumsum(exp(-smooth / 2)) / sum(exp(-smooth / 2))

#Draw some typical bandwidths from the CDF and produce a model with each.
bws = xxx[sapply(runif(S), function(x) which(x<pp)[1])]
models.lagr = list()
models.lagr[[1]] = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, varselect.method='AICc', bw=0.15, kernel=epanechnikov, bw.type='knn', verbose=TRUE)
for (bw in bws) {
    indx = sample(1:nrow(sim), replace=TRUE)
    boot = sim[indx,]
    models.lagr[[length(models.lagr)+1]] = lagr(Y~X1+X2+X3+X4+X5-1, data=boot, family='gaussian', fit.loc=sim[,c('loc.x','loc.y')], coords=boot[,c('loc.x','loc.y')], longlat=FALSE, varselect.method='AICc', bw=bw, kernel=epanechnikov, bw.type='knn', verbose=TRUE)
}



#ORACLE:
write.log('making oracular model.', 'result.txt')
bw.oracle = lagr.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, oracle=oracle, kernel=epanechnikov, tol.bw=0.01, bw.type='knn', verbose=TRUE, bwselect.method='AICc', resid.type='pearson')

#Make the bw distribution so we can draw from it.
#First, get the AICc-vs-bw results in order
b = bw.oracle[['trace']][order(bw.oracle[['trace']][,1]),c(1,2)]

#Fit a spline through the AICc-vs-bw observations and then use it to smooth AICc across the entire range of the tested bandwidths.
spline = smooth.spline(b)
xxx = seq(b[1,1], tail(b[,1],1), length.out=1001)
smooth = predict(spline, xxx)
smooth = smooth$y - mean(smooth$y)

#Now restrict our attention to the region of the densest 99% of bandwidth probability mass
maxi = which(cumsum(exp(-smooth / 2)) / sum(exp(-smooth / 2)) > 0.995)[1]
mini = tail(which(cumsum(exp(-smooth / 2))/sum(exp(-smooth / 2)) < 0.005),1)
xxx = seq(xxx[mini], xxx[maxi], length.out=1001)
smooth = predict(spline, xxx)$y

#Get the CDF of bandwidth within the region of greatest density
pp = cumsum(exp(-smooth / 2)) / sum(exp(-smooth / 2))

#Draw some typical bandwidths from the CDF and produce a model with each.
bws = xxx[sapply(runif(S), function(x) which(x<pp)[1])]
models.oracle = list()
models.oracle[[1]] = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, oracle=oracle, bw=bw.oracle[['bw']], kernel=epanechnikov, bw.type='knn', verbose=TRUE)
for (bw in bws) {
    indx = sample(1:nrow(sim), replace=TRUE)
    boot = sim[indx,]
    models.oracle[[length(models.oracle)+1]] = lagr(Y~X1+X2+X3+X4+X5-1, data=boot, family='gaussian', fit.loc=sim[,c('loc.x','loc.y')], coords=boot[,c('loc.x','loc.y')], longlat=FALSE, oracle=oracle, bw=bw, kernel=epanechnikov, bw.type='knn', verbose=TRUE)
}



#GWR:
write.log('making gwr model.', 'result.txt')
bw.gwr = lagr.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, oracle=allvars, kernel=epanechnikov, tol.bw=0.01, bw.type='knn', verbose=TRUE, bwselect.method='AICc', resid.type='pearson')

#Make the bw distribution so we can draw from it.
#First, get the AICc-vs-bw results in order
b = bw.gwr[['trace']][order(bw.gwr[['trace']][,1]),c(1,2)]

#Fit a spline through the AICc-vs-bw observations and then use it to smooth AICc across the entire range of the tested bandwidths.
spline = smooth.spline(b)
xxx = seq(b[1,1], tail(b[,1],1), length.out=1001)
smooth = predict(spline, xxx)
smooth = smooth$y - mean(smooth$y)

#Now restrict our attention to the region of the densest 99% of bandwidth probability mass
maxi = which(cumsum(exp(-smooth / 2)) / sum(exp(-smooth / 2)) > 0.995)[1]
mini = tail(which(cumsum(exp(-smooth / 2))/sum(exp(-smooth / 2)) < 0.005),1)
xxx = seq(xxx[mini], xxx[maxi], length.out=1001)
smooth = predict(spline, xxx)$y

#Get the CDF of bandwidth within the region of greatest density
pp = cumsum(exp(-smooth / 2)) / sum(exp(-smooth / 2))

#Draw some typical bandwidths from the CDF and produce a model with each.
bws = xxx[sapply(runif(S), function(x) which(x<pp)[1])]
models.gwr = list()
models.gwr[[1]] = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, oracle=allvars, bw=bw.gwr[['bw']], kernel=epanechnikov, bw.type='knn', verbose=TRUE)
for (bw in bws) {
    indx = sample(1:nrow(sim), replace=TRUE)
    boot = sim[indx,]
    models.gwr[[length(models.gwr)+1]] = lagr(Y~X1+X2+X3+X4+X5-1, data=boot, family='gaussian', fit.loc=sim[,c('loc.x','loc.y')], coords=boot[,c('loc.x','loc.y')], longlat=FALSE, oracle=allvars, bw=bw, kernel=epanechnikov, bw.type='knn', verbose=TRUE)
}


#OUTPUT:

#LAGR:
write.log('summarize LAGR model.', 'result.txt')
save(bw.lagr, file=paste("bw.", cluster, ".", process, ".lagr.RData", sep=""))
save(models.lagr, file=paste("models.", cluster, ".", process, ".lagr.RData", sep=""))





#oracle:
write.log('summarize oracle model.', 'result.txt')
save(bw.oracle, file=paste("bw.", cluster, ".", process, ".oracle.RData", sep=""))
save(models.oracle, file=paste("models.", cluster, ".", process, ".oracle.RData", sep=""))





#GWR:
write.log('summarize gwr model.', 'result.txt')
save(bw.gwr, file=paste("bw.", cluster, ".", process, ".gwr.RData", sep=""))
save(models.gwr, file=paste("models.", cluster, ".", process, ".gwr.RData", sep=""))

write.log('done.', 'result.txt')
