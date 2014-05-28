covariates = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')

#For each output, make a list where each covariate has a vector
mu = sapply(covariates, function(x) return(vector()), USE.NAMES=TRUE)
med = sapply(covariates, function(x) return(vector()), USE.NAMES=TRUE)
sigma = sapply(covariates, function(x) return(vector()), USE.NAMES=TRUE)
qq = sapply(covariates, function(x) return(vector()), USE.NAMES=TRUE)
pzero = sapply(covariates, function(x) return(vector()), USE.NAMES=TRUE)

truth = list('(Intercept)'=rep(0,900), X1=B1, X2=rep(0,900), X3=rep(0,900), X4=rep(0,900), X5=rep(0,900))

for (l in 1:900) {
    cat(paste("location: ", l, "\n", sep=""))
    coef_loc = coefs[coefs$location==l,]
    
    #Compute the summary diagnostics
    for (v in covariates) {
        mu[[v]] = c(mu[[v]], mean(coef_loc[,v]))
        med[[v]] = c(med[[v]], median(coef_loc[,v]))
        sigma[[v]] = c(sigma[[v]], sd(coef_loc[,v]))
        qq[[v]] = c(qq[[v]], ecdf(coef_loc[,v])(truth[[v]][l]))
        pzero[[v]] = c(pzero[[v]], sum(coef_loc[,v]==0)/2000)
    }
}