coefs = matrix(NA, nrow=0, ncol=9)

for (i in 1:length(models)) {
    coefs_temp = t(sapply(models[[i]][['model']][['models']], function(x) x[['coef']]))
    coefs_temp = cbind(iter=i, bw=models[[i]][['bw']], loc=1:900, coefs_temp)
    coefs = rbind(coefs, coefs_temp)
}

coefs = data.frame(coefs)
colnames(coefs)[4] = "(Intercept)"