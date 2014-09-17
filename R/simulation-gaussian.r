library(lagr)
library(RandomFields)

if (interactive()) {
    process = 6
} else {
    args <- commandArgs(trailingOnly = TRUE)
    process = as.numeric(args[1])
}
print(process)

#Establish the simulation parameters
tau = rep(0.1, 18) #tau is the spatial autocorrelation range parameter for the covariates
rho = rep(c(rep(0, 2), rep(0.5, 2), rep(0.9, 2)), 6) #rho is the correlation of the covariates
sigma.tau = rep(0, 18) #sigma.tau is the spatial autocorrelation range parameter of the noise term
sigma = rep(c(0.5, 1), 9) #sigma is the variance of the noise term
size = c(rep(100, 6), rep(200, 6), rep(400, 6))
params = data.frame(tau, rho, sigma.tau, sigma, size)

#Simulation parameters are based on the value of process
parameters = params[process,]

#Seed the RNG
set.seed(process %% 6)

#Generate the covariates:
N = 20 #number of width and length divisions in the domain
coord = seq(0, 1, length.out=N)
if (parameters[['tau']] > 0) {
    covariate.model = RMexp(var=1, scale=parameters[['tau']]) + RMnugget(var=0.2)
    d1 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
    d2 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
    d3 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
    d4 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
    d5 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
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

#Simulate the noise term, either with or without spatial correlation
if (parameters[['sigma.tau']] == 0) {epsilon = rnorm(N**2, mean=0, sd=parameters[['sigma']])}
if (parameters[['sigma.tau']] > 0) {epsilon = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(parameters[['sigma']]**2,parameters[['sigma.tau']]))$data}

#Calculate the coefficient surfaces
B1 = RFsimulate(RMexp(var=2.5, scale=1), x=coord, y=coord)@data[[1]]
B2 = RFsimulate(RMexp(var=0.5, scale=1), x=coord, y=coord)@data[[1]]
B3 = RFsimulate(RMexp(var=0.1, scale=1), x=coord, y=coord)@data[[1]]
B4 = RFsimulate(RMexp(var=0.02, scale=1), x=coord, y=coord)@data[[1]]

#Generate the response variable and set up the data.frame:
eta = X1*B1 + X2*B2 + X3+B3 + X4*B4
Y = eta + epsilon
sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), X5=as.vector(X5), loc.x, loc.y)

indx = sample(400, parameters[['size']])
data = sim[indx,]
h = 0.5*parameters[['size']]**(-1/6)

write.table(cbind(x=loc.x[indx], y=loc.y[indx], B1=B1[indx], B2=B2[indx], B3=B3[indx], B4=B4[indx], Y=Y[indx]), file=paste("output/truth", process, "csv", sep="."))

#LAGR:
model = lagr(Y~X1+X2+X3+X4+X5, data=data, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', bw=h, kernel=epanechnikov, bw.type='knn', verbose=TRUE)

#Write LAGR coefficients:
coefs = t(sapply(model[['model']], function(x) x[['coef']]))
write.table(cbind(x=loc.x[indx], y=loc.y[indx], coefs), file=paste("output/coefs", "lagr", process, "csv", sep="."))

#Write LAGR results:
fits = as.vector(sapply(model[['model']], function(x) x[['fitted']]))
resids = Y[indx] - fits
result = data.frame(x=loc.x[indx], y=loc.y[indx], fitted=fits, residual=resids)
write.table(result, file=paste("output/results", "lagr", process, "csv", sep="."))

#Clean up:
rm(model)
rm(coefs)
rm(result)
gc()


#GWR:
allvars = replicate(length(indx)**2, c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5'), simplify=FALSE)
gwr = lagr(Y~X1+X2+X3+X4+X5, data=data, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, oracle=allvars, bw=0.h, kernel=epanechnikov, bw.type='knn', verbose=TRUE)

#Write the GWR coefficients
cgwr = t(sapply(gwr[['model']], function(x) x[['coef']]))
write.table(cbind(x=loc.x[indx], y=loc.y[indx], cgwr), file=paste("output/coefs", "gwr", process, "csv", sep="."))

#Write the GWR results:
fgwr = as.vector(sapply(gwr[['model']], function(x) x[['fitted']]))
rgwr = Y[indx] - fgwr
result = data.frame(x=loc.x[indx], y=loc.y[indx], fitted=fgwr, residual=rgwr)
write.table(result, file=paste("output/results", "gwr", process, "csv", sep="."))

#Clean up:
rm(gwr)
rm(cgwr)
rm(result)
gc()
