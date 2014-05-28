write.log = function(message, file, append=TRUE) {
    sink(file, append=append)
    cat(paste(message, "\n", sep=""))
    sink()
}

write.log('entry', 'result.txt', append=FALSE)

#Set ourselves up to import the packages:
r = getOption("repos")
r["CRAN"] = "http://cran.wustl.edu"
options(repos = r)
rm(r)
dir.create("rlibs")
Sys.setenv(R_LIBS="rlibs")
.libPaths(new="rlibs")

install.packages("sp")
install.packages("foreach")
install.packages("iterators")
install.packages("multicore")
install.packages("doMC")
install.packages("geoR")
install.packages("glmnet")

require(sp)
require(splancs)
require(foreach)
require(iterators)
require(multicore)
require(geoR)
require(glmnet)
require(SGL)
require(lagr)

write.log('installations complete', 'result.txt')

seeds = as.vector(read.csv("seeds.txt", header=FALSE)[,1])
B = 100
N = 30
settings = 18
functions = 3
coord = seq(0, 1, length.out=N)

#Establish the simulation parameters
tau = rep(0.1, settings)
rho = rep(c(rep(0,2), rep(0.5,2), rep(0.9,2)), functions)
sigma.tau = rep(0, settings)
sigma = rep(c(0.5,1), settings/2)
params = data.frame(tau, rho, sigma.tau, sigma)

#Read the cluster and process arguments
args = scan('jobid.txt', 'character')
args = strsplit(args, '\\n', fixed=TRUE)[[1]]
cluster = args[1]
process = as.integer(args[2]) - 1

write.log(paste('process:', process, sep=''), 'result.txt')

#Simulation parameters are based on the value of process
setting = process %/% B + 1
parameters = params[setting,]
set.seed(seeds[process+1])

write.log(paste('seed:', seeds[process+1], sep=''), 'result.txt')

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

write.log('generated GRFs', 'result.txt')

loc.x = rep(coord, times=N)
loc.y = rep(coord, each=N)

#Use the Cholesky decomposition to correlate the random fields:
S = matrix(parameters[['rho']], 5, 5)
diag(S) = rep(1, 5)
L = chol(S)

write.log('generated correlation matrix', 'result.txt')

#Force correlation on the Gaussian random fields:
D = as.matrix(cbind(d1, d2, d3, d4, d5)) %*% L
    
#
X1 = matrix(D[,1], N, N)
X2 = matrix(D[,2], N, N)
X3 = matrix(D[,3], N, N)
X4 = matrix(D[,4], N, N)
X5 = matrix(D[,5], N, N)

write.log('generated Xs', 'result.txt')

if (parameters[['sigma.tau']] == 0) {epsilon = rnorm(N**2, mean=0, sd=parameters[['sigma']])}
if (parameters[['sigma.tau']] > 0) {epsilon = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(parameters[['sigma']]**2,parameters[['sigma.tau']]))$data}

write.log('generated epsilon', 'result.txt')

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

write.log('generated B1 coefficient surface', 'result.txt')

mu = X1*B1
Y = mu + epsilon

sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), X5=as.vector(X5), loc.x, loc.y)

vars = cbind(B1=as.vector(B1!=0))#, B2=as.vector(B2!=0), B3=as.vector(B3!=0))
oracle = list()
for (i in 1:N**2) { 
    oracle[[i]] = character(0)
    if (vars[i,'B1']) { oracle[[i]] = c(oracle[[i]] , "X1") }
}

write.log('generated data', 'result.txt')


#MODELS:
write.log('making lagr model.', 'result.txt')
bw.lagr = lagr.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, varselect.method='AICc', kernel=epanechnikov, tol.bw=0.01, bw.type='knn', verbose=FALSE, bwselect.method='AICc', resid.type='pearson')

#Make the bw distribution so we can draw from it.
#First, get the AICc-vs-bw results in order
b = bw.lagr[['trace']][order(bw.lagr[['trace']][,1]),c(1,2)]

#Fit a spline through the AICc-vs-bw observations and then use it to smooth AICc across the entire range of the tested bandwidths.
spline = smooth.spline(b)
xxx = seq(b[1,1], tail(b[,1],1), length.out=1001)
smooth = predict(spline, xxx)

#Now restrict our attention to the region of the densest 99.98% of bandwidth probability mass
maxi = which(cumsum(exp(-0.5*smooth$y))/sum(exp(-0.5*smooth$y))>0.9999)[1]
mini = tail(which(cumsum(exp(-0.5*smooth$y))/sum(exp(-0.5*smooth$y))<0.0001),1)
xxx = seq(xxx[mini], xxx[maxi], length.out=1001)
smooth = predict(spline, xxx)

#Get the CDF of bandwidth within the region of greatest density
pp = cumsum(exp(-0.5*smooth$y))/sum(exp(-0.5*smooth$y))

#Draw some typical bandwidths from the CDF and produce a model with each.
bws = xxx[sapply(runif(19), function(x) which(x<pp)[1])]
models = list()
models[[1]] = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], N=100, longlat=FALSE, varselect.method='AICc', bw=bw.lagr[['bw']], kernel=epanechnikov, bw.type='knn', simulation=TRUE, verbose=FALSE)
for (bw in bws) {
    models[[length(models)+1]] = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], N=100, longlat=FALSE, varselect.method='AICc', bw=bw, kernel=epanechnikov, bw.type='knn', simulation=TRUE, verbose=FALSE)
}


write.log('making oracular model.', 'result.txt')
bw.oracular = gwglmnet.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', oracle=oracle, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, kernel=epanechnikov, tol.bw=0.01, bw.type='knn', verbose=FALSE, bwselect.method='AICc', resid.type='pearson')
model.oracular = gwglmnet(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', oracle=oracle, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, bw=bw.oracular[['bw']], kernel=epanechnikov, bw.type='knn', simulation=TRUE, verbose=FALSE)

#write.log('making gwr model.', 'result.txt')
#oracle2 = lapply(1:nrow(sim), function(x) {return(c("X1", "X2", "X3", "X4", "X5"))})
#bw.gwr = gwglmnet.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', oracle=oracle, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select='BIC', gweight=spherical, tol.bw=0.01, bw.method='knn', parallel=FALSE, interact=TRUE, verbose=FALSE, shrunk.fit=FALSE, bw.select='AICc', resid.type='pearson')
#model.gwr = gwglmnet(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', oracle=oracle2, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, mode.select='BIC', bw=bw.oracular[['bw']], gweight=spherical, bw.method='knn', simulation=TRUE, parallel=FALSE, interact=TRUE, verbose=FALSE, shrunk.fit=FALSE)



#OUTPUT:

#First, write the data
write.log('write the data.', 'result.txt')
write.table(sim, file=paste("Data.", cluster, ".", process, ".csv", sep=""), sep=',', row.names=FALSE)


#LAGR:
write.log('summarize LAGR model.', 'result.txt')
write.log(paste('LAGR bandwidth: ', bw.lagr[['bw']], sep=''), 'result.txt')
vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')

coefs = t(sapply(model.lagr[['model']][['models']], function(x) as.vector(x[['coef']])))
write.table(coefs, file=paste("CoefEstimates.", cluster, ".", process, ".lagr.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

output = matrix(NA,0,2)
output = t(sapply(model.lagr[['model']][['models']], function(x) c(x[['sigma2']][x[['s']]], x[['fitted']])))
colnames(output) = c("s2", "fitted")
write.table(output, file=paste("MiscParams.", cluster, ".", process, ".lagr.csv", sep=""), sep=',', row.names=FALSE)






#For oracle property:
write.log('summarize oracle model.', 'result.txt')
write.log(paste('Oracle bandwidth: ', bw.oracular[['bw']], sep=''), 'result.txt')
vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')

coefs = t(sapply(model.oracular[['model']][['models']], function(x) as.vector(x[['coef']])))
write.table(coefs, file=paste("CoefEstimates.", cluster, ".", process, ".oracle.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

#params = c('bw', 'sigma2', 'loss.local', 'fitted')
#target = params[1]
#output = sapply(1:N**2, function(y) {model.oracular[['model']][['models']][[y]][[target]]})

#for (i in 2:length(params)) {
#    target = params[i]
#    output = cbind(output, sapply(1:N**2, function(y) {model.oracular[['model']][['models']][[y]][[target]]}))
#}
#write.table(output, file=paste("MiscParams.", cluster, ".", process, ".oracle.csv", sep=""), col.names=params, sep=',', row.names=FALSE)






#For all vars:
#write.log('summarize GWR-LLE model.', 'result.txt')
#vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')

#coefs = t(sapply(1:N**2, function(y) {as.vector(model.gwr[['model']][['models']][[y]][['coef']])}))
#write.table(coefs, file=paste("CoefEstimates.", cluster, ".", process, ".gwr.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

#params = c('bw', 'sigma2', 'loss.local', 'fitted')
#target = params[1]
#output = sapply(1:N**2, function(y) {model.gwr[['model']][['models']][[y]][[target]]})

#for (i in 2:length(params)) {
#    target = params[i]
#    output = cbind(output, sapply(1:N**2, function(y) {model.gwr[['model']][['models']][[y]][[target]]}))
#}
#write.table(output, file=paste("MiscParams.", cluster, ".", process, ".gwr.csv", sep=""), #col.names=params, sep=',', row.names=FALSE)

write.log('done.', 'result.txt')
