write.log = function(message, file, append=TRUE) {
    sink(file, append=append)
    cat(paste(message, "\n", sep=""))
    sink()
}

#Set the location from which to load packages
Sys.setenv(R_LIBS="rlibs")
.libPaths(new="rlibs")

require(sp)
require(splancs)
require(foreach)
#require(iterators)
#require(multicore)
require(geoR)
require(SGL)
require(lagr)
#require(doMC)
#registerDoMC(7)

#Read the process number from the command line
args <- commandArgs(trailingOnly = TRUE)
cluster = NA
process = as.numeric(args[1])

logfile = paste("result.", process, ".txt", sep="")
write.log('installations complete', logfile, append=FALSE)

B = 100 #number of replications for each setting
N = 30 #number of width and length divisions in the domain
settings = 18 #number of unique simulation settings
functions = 3 #number of function types
coord = seq(0, 1, length.out=N) #coordinates of the generated observations

#Establish the simulation parameters
tau = rep(0.1, settings) #tau is the spatial autocorrelation range parameter for the covariates
rho = rep(c(rep(0,2), rep(0.5,2), rep(0.9,2)), functions) #rho is the correlation of the covariates
sigma.tau = rep(0, settings) #sigma.tau is the spatial autocorrelation range parameter of the noise term
sigma = rep(c(0.5,1), settings/2) #sigma is the variance of the noise term
params = data.frame(tau, rho, sigma.tau, sigma)

#Read the cluster and process arguments
write.log(paste('process:', process, sep=''), logfile)

#Simulation parameters are based on the value of process
setting = process %/% B + 1
parameters = params[setting,]
set.seed(process)
write.log(paste('seed: ', process, sep=''), logfile)
write.log(paste("parameters: ", paste(parameters, collapse=","), sep=""), logfile)

#Generate the covariates:
if (parameters[['tau']] > 0) {
    write.log(paste('generating GRFs with tau of ', parameters[['tau']], sep=''), logfile)
    d1 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))$data
    d2 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))$data
    d3 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))$data
    d4 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))$data
    d5 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))$data
} else {
    write.log('generating GRFs with tau of 0', logfile)
    d1 = rnorm(N**2, mean=0, sd=1)
    d2 = rnorm(N**2, mean=0, sd=1)
    d3 = rnorm(N**2, mean=0, sd=1)
    d4 = rnorm(N**2, mean=0, sd=1)
    d5 = rnorm(N**2, mean=0, sd=1)
}
loc.x = rep(coord, times=N)
loc.y = rep(coord, each=N)
write.log('generated GRFs', logfile)

#Use the Cholesky decomposition to correlate the random fields:
S = matrix(parameters[['rho']], 5, 5)
diag(S) = rep(1, 5)
L = chol(S)
write.log('generated correlation matrix', logfile)

#Force correlation on the Gaussian random fields:
D = as.matrix(cbind(d1, d2, d3, d4, d5)) %*% L
X1 = matrix(D[,1], N, N)
X2 = matrix(D[,2], N, N)
X3 = matrix(D[,3], N, N)
X4 = matrix(D[,4], N, N)
X5 = matrix(D[,5], N, N)
write.log('generated Xs', logfile)

#simulate the noise term, either with or without spatial correlation
if (parameters[['sigma.tau']] == 0) {epsilon = rnorm(N**2, mean=0, sd=parameters[['sigma']])}
if (parameters[['sigma.tau']] > 0) {epsilon = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(parameters[['sigma']]**2,parameters[['sigma.tau']]))$data}
write.log('generated epsilon', logfile)

#calculate the B1 coefficient surface for the appropriate function type
if ((setting-1) %/% 6 == 0) {
    B1 = matrix(rep(ifelse(coord<=0.4, 0, ifelse(coord<0.6,5*(coord-0.4),1)), N), N, N)
} else if ((setting-1) %/% 6 == 1) {
    B1 = matrix(rep(coord, N), N, N)
} else if ((setting-1) %/% 6 == 2) {
    Xmat = matrix(rep(rep(coord, times=N), times=N), N**2, N**2)
    Ymat = matrix(rep(rep(coord, each=N), times=N), N**2, N**2)
    D = (Xmat-0.5)**2 + (Ymat-0.5)**2
    d = D[,435]
    B1 = matrix(max(d)-d, N, N)
    B1 = B1 / max(B1)
}
write.log('generated B1 coefficient surface', logfile)

#Generate the response variable and set up the data.frame:
mu = X1*B1
Y = mu + epsilon
sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), X5=as.vector(X5), loc.x, loc.y)
write.log('generated data', logfile)
write.table(sim, file=paste("simdata", process, "csv", sep="."))


#MODELS:
#LAGR:
write.log('making lagr model.', logfile)
#bw.lagr = lagr.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, varselect.method='AICc', kernel=epanechnikov, tol.bw=0.01, bw.type='knn', verbose=TRUE, bwselect.method='AICc', resid.type='pearson')
#save(bw.lagr, file=paste("bw", "lagr", process, "RData", sep="."))

#Draw some typical bandwidths from the CDF and produce a model with each.
#bws = interpolate.bw(bw.lagr[['trace']], S=20)

#i = 0
#model = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, varselect.method='AICc', bw=bw.lagr[['bw']], kernel=epanechnikov, bw.type='knn', verbose=TRUE)
model = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, varselect.method='AICc', bw=0.25, kernel=epanechnikov, bw.type='knn', verbose=TRUE)
coefs = t(sapply(model[['model']][['models']], function(x) x[['coef']]))
write.table(coefs, file=paste("coefs", "lagr", process, "csv", sep="."))

#save(model, file=paste("lagr.model.", i, ".RData", sep=""))
rm(model)
rm(coefs)
gc()

#for (bw in bws) {
#    i = i+1
#    indx = sample(1:nrow(sim), replace=TRUE)
#    boot = sim[indx,]
#    model = lagr(Y~X1+X2+X3+X4+X5-1, data=boot, family='gaussian', fit.loc=sim[,c('loc.x','loc.y')], coords=boot[,c('loc.x','loc.y')], longlat=FALSE, varselect.method='AICc', bw=bw, kernel=epanechnikov, bw.type='knn', verbose=TRUE)
#
#    save(model, file=paste("lagr.model.", i, ".RData", sep=""))
#    rm(model)
#    gc()
#}



#ORACLE:
write.log('making oracular model.', logfile)
vars = cbind(B1=as.vector(B1!=0))
oracle = list()
for (i in 1:N**2) { 
    oracle[[i]] = character(0)
    if (vars[i,'B1']) { oracle[[i]] = c(oracle[[i]] , "X1") }
}
#bw.oracle = lagr.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, oracle=oracle, kernel=epanechnikov, tol.bw=0.01, bw.type='knn', verbose=TRUE, bwselect.method='AICc', resid.type='pearson')
#save(bw.oracle, file=paste("bw", "oracle", process, "RData", sep="."))

#Draw some typical bandwidths
#bws = interpolate.bw(bw.oracle[['trace']], S=20)

#Generate the models
#i=0
#model = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, oracle=oracle, bw=bw.oracle[['bw']], kernel=epanechnikov, bw.type='knn', verbose=TRUE)
model = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, oracle=oracle, bw=0.25, kernel=epanechnikov, bw.type='knn', verbose=TRUE)
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
#bw.gwr = lagr.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, oracle=allvars, kernel=epanechnikov, tol.bw=0.01, bw.type='knn', verbose=TRUE, bwselect.method='AICc', resid.type='pearson')
#save(bw.gwr, file=paste("bw", "gwr", process, "RData", sep="."))

#Draw some typical bandwidths from the CDF and produce a model with each.
#bws = interpolate.bw(bw.gwr[['trace']], S=20)

#i=0
#model = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, oracle=allvars, bw=bw.gwr[['bw']], kernel=epanechnikov, bw.type='knn', verbose=TRUE)
model = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, oracle=allvars, bw=0.25, kernel=epanechnikov, bw.type='knn', verbose=TRUE)
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
