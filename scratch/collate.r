coefs = matrix(NA, nrow=0, ncol=8)

for (i in 1:length(models)) {
    coefs_temp = matrix(NA, nrow=0, ncol=7)
    
    for (j in 1:length(models[[i]][['model']][['models']])) {
        coefs_loc = matrix(NA, nrow=0, ncol=6)
        
        for (m in models[[i]][['model']][['models']][[j]][['coeflist']]) {
            coefs_loc = rbind(coefs_loc, as.vector(m))
        }
          
        #Add the location indicator to coefs_loc and then add coefs_loc to coefs.
        coefs_loc = cbind(j, coefs_loc)
        #colnames(coefs_loc) = NULL
        coefs_temp = rbind(coefs_temp, coefs_loc)
    }
    
    coefs_temp = cbind(i, coefs_temp)
    #colnames(coefs_temp) = NULL
    coefs = rbind(coefs, coefs_temp)
}

coefs = data.frame(coefs)
colnames(coefs) = c("iteration", "location", "(Intercept)", "X1", "X2", "X3", "X4", "X5")