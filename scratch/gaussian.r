library(lagr)
require(sp)
require(splancs)
require(foreach)
require(geoR)
require(SGL)
library(lattice)

require(doMC)
registerDoMC(3)

process = 1
cluster = NA
N = 20 #number of width and length divisions in the domain
settings = 18 #number of unique simulation settings
functions = 3 #number of function types
coord = seq(0, 1, length.out=N) #coordinates of the generated observations

#Establish the simulation parameters
tau = rep(0.1, settings) #tau is the spatial autocorrelation range parameter for the covariates
rho = rep(c(rep(0,2), rep(0.5,2), rep(0.9,2)), functions) #rho is the correlation of the covariates
sigma.tau = rep(0, settings) #sigma.tau is the spatial autocorrelation range parameter of the noise term
sigma = rep(c(0.5,1), settings/2) #sigma is the variance of the noise term
params = data.frame(tau, rho, sigma.tau, sigma)

#Simulation parameters are based on the value of process
parameters = params[process,]

#Seed the RNG
set.seed(8)

#Generate the covariates:
if (parameters[['tau']] > 0) {
    d1 = RFsimulate(RMexp(var=1, scale=parameters[['tau']]), x=coord, y=coord)@data[[1]]
    d2 = RFsimulate(RMexp(var=1, scale=parameters[['tau']]), x=coord, y=coord)@data[[1]]
    d3 = RFsimulate(RMexp(var=1, scale=parameters[['tau']]), x=coord, y=coord)@data[[1]]
    d4 = RFsimulate(RMexp(var=1, scale=parameters[['tau']]), x=coord, y=coord)@data[[1]]
    d5 = RFsimulate(RMexp(var=1, scale=parameters[['tau']]), x=coord, y=coord)@data[[1]]
} else {
    d1 = rnorm(N**2, mean=0, sd=1)
    d2 = rnorm(N**2, mean=0, sd=1)
    d3 = rnorm(N**2, mean=0, sd=1)
    d4 = rnorm(N**2, mean=0, sd=1)
    d5 = rnorm(N**2, mean=0, sd=1)
}
loc.x = rep(coord, times=N)
loc.y = rep(coord, each=N)

#Use the Cholesky decomposition to correlate the random fields:
S = matrix(parameters[['rho']], 5, 5)
diag(S) = rep(1, 5)
L = chol(S)

#Force correlation on the Gaussian random fields:
D = as.matrix(cbind(d1, d2, d3, d4, d5)) %*% L
X1 = matrix(D[,1], N, N)
X2 = matrix(D[,2], N, N)
X3 = matrix(D[,3], N, N)
X4 = matrix(D[,4], N, N)
X5 = matrix(D[,5], N, N)

#simulate the noise term, either with or without spatial correlation
if (parameters[['sigma.tau']] == 0) {epsilon = rnorm(N**2, mean=0, sd=parameters[['sigma']])}
if (parameters[['sigma.tau']] > 0) {epsilon = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(parameters[['sigma']]**2,parameters[['sigma.tau']]))$data}

#calculate the B1 coefficient surface for the appropriate function type
B1 = RFsimulate(RMexp(var=2.5, scale=1), x=coord, y=coord)@data[[1]]
B2 = RFsimulate(RMexp(var=0.5, scale=1), x=coord, y=coord)@data[[1]]
B3 = RFsimulate(RMexp(var=0.1, scale=1), x=coord, y=coord)@data[[1]]
B4 = RFsimulate(RMexp(var=0.02, scale=1), x=coord, y=coord)@data[[1]]

#if ((process-1) %/% 6 == 0) {
#    B1 = matrix(rep(ifelse(coord<=0.4, 0, ifelse(coord<0.6,5*(coord-0.4),1)), N), N, N)
#} else if ((process-1) %/% 6 == 1) {
#    B1 = matrix(rep(coord, N), N, N)
#} else if ((process-1) %/% 6 == 2) {
#    Xmat = matrix(rep(rep(coord, times=N), times=N), N**2, N**2)
#    Ymat = matrix(rep(rep(coord, each=N), times=N), N**2, N**2)
#    D = (Xmat-0.5)**2 + (Ymat-0.5)**2
#    d = D[,435]
#    B1 = matrix(max(d)-d, N, N)
#    B1 = B1 / max(B1)
#}

#Generate the response variable and set up the data.frame:
eta = X1*B1
p = exp(eta) / (1+exp(eta))
mu = exp(eta)

#binomial, gaussian, poisson response:
Y = rbinom(N**2, p=p, size=1)
Y = eta + epsilon
Y = rpois(N**2, mu)

sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), X5=as.vector(X5), loc.x, loc.y)
fl = cbind(loc.x, loc.y)[sample(1:N**2, 10),]

#MODELS:
S = 100 #number of draws from the bandwidth distribution for each replication
model = lagr(Y~X1+X2+X3+X4+X5, data=sim, family=gaussian, coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', bw=0.2, kernel=epanechnikov, bw.type='knn', verbose=TRUE)

stop()

cc = sapply(model[['model']], function(x) x[['coef']])
wireframe(matrix(cc[2,], 30, 30))


#LAGR:
write.log('making lagr model.', logfile)
#bw.lagr = lagr.tune(Y~X1+X2+X3+X4+X5, data=sim, family='binomial', coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', kernel=epanechnikov, tol.bw=0.01, bw.type='knn', verbose=TRUE, bwselect.method='AIC', resid.type='deviance')
#save(bw.lagr, file=paste("bw.", cluster, ".", process, ".lagr.RData", sep=""))

#Draw some typical bandwidths from the CDF and produce a model with each.
#bws = interpolate.bw(bw.lagr[['trace']], S=20)

#Now restrict our attention to the region of the densest 99% of bandwidth probability mass
maxi = min(1001, which(cumsum(exp(-smooth / 2)) / sum(exp(-smooth / 2)) > 0.995)[1])
mini = max(1, tail(which(cumsum(exp(-smooth / 2))/sum(exp(-smooth / 2)) < 0.005),1))
xxx = seq(xxx[mini], xxx[maxi], length.out=1001)
smooth = predict(spline, xxx)$y

#save(model, file=paste("lagr.model.", i, ".RData", sep=""))
rm(model)
rm(coefs)
gc()


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
write.log('making oracular model.', logfile)
vars = cbind(B1=as.vector(B1!=0))
oracle = list()
for (i in 1:N**2) { 
    oracle[[i]] = character(0)
    if (vars[i,'B1']) { oracle[[i]] = c(oracle[[i]] , "X1") }
}
bw.oracle = lagr.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, oracle=oracle, kernel=epanechnikov, tol.bw=0.01, bw.type='knn', verbose=TRUE, bwselect.method='AICc', resid.type='pearson')
#save(bw.oracle, file=paste("bw.", cluster, ".", process, ".oracle.RData", sep=""))

#Draw some typical bandwidths
#bws = interpolate.bw(bw.oracle[['trace']], S=20)

#Generate the models
#i=0
model = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, oracle=oracle, bw=bw.oracle[['bw']], kernel=epanechnikov, bw.type='knn', verbose=TRUE)
coefs = t(sapply(model[['model']][['models']], function(x) x[['coef']]))
write.table(coefs, file=paste("coefs", "oracle", process, "csv", sep="."))

#save(model, file=paste("oracle.model.", i, ".RData", sep=""))
rm(model)
rm(coefs)
gc()

#for (bw in bws) {
#    i = i+1
#    indx = sample(1:nrow(sim), replace=TRUE)
#    boot = sim[indx,]
#    model = lagr(Y~X1+X2+X3+X4+X5-1, data=boot, family='gaussian', fit.loc=sim[,c('loc.x','loc.y')], coords=boot[,c('loc.x','loc.y')], longlat=FALSE, oracle=oracle, bw=bw, kernel=epanechnikov, bw.type='knn', verbose=TRUE)
#
#    save(model, file=paste("oracle.model.", i, ".RData", sep=""))
#    rm(model)
#    gc()
#}



#GWR:
write.log('making gwr model.', logfile)
allvars = replicate(N**2, c('X1', 'X2', 'X3', 'X4', 'X5'), simplify=FALSE)
bw.gwr = lagr.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, oracle=allvars, kernel=epanechnikov, tol.bw=0.01, bw.type='knn', verbose=TRUE, bwselect.method='AICc', resid.type='pearson')
#save(bw.gwr, file=paste("bw.", cluster, ".", process, ".gwr.RData", sep=""))

#Draw some typical bandwidths from the CDF and produce a model with each.
#bws = interpolate.bw(bw.gwr[['trace']], S=20)

#i=0
model = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, oracle=allvars, bw=bw.gwr[['bw']], kernel=epanechnikov, bw.type='knn', verbose=TRUE)
coefs = t(sapply(model[['model']][['models']], function(x) x[['coef']]))
write.table(coefs, file=paste("coefs", "gwr", process, "csv", sep="."))

#save(model, file=paste("gwr.model.", i, ".RData", sep=""))
rm(model)
rm(coefs)
gc()

#for (bw in bws) {
#    i = i+1
#    indx = sample(1:nrow(sim), replace=TRUE)
#    boot = sim[indx,]
#    model = lagr(Y~X1+X2+X3+X4+X5-1, data=boot, family='gaussian', fit.loc=sim[,c('loc.x','loc.y')], coords=boot[,c('loc.x','loc.y')], longlat=FALSE, oracle=allvars, bw=bw, kernel=epanechnikov, bw.type='knn', verbose=TRUE)
#
#    save(model, file=paste("gwr.model.", i, ".RData", sep=""))
#    rm(model)
#    gc()
#}

#We're finished
write.log('done.', logfile)
