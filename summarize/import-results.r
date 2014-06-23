vars = c("X.Intercept.", "X1", "X2", "X3", "X4", "X5")
models = c("lagr", "oracle", "gwr")
settings = 1:18

output = sapply(settings, function(y)
    sapply(models, function(t) 
        sapply(vars, function(x)
            return(matrix(ncol=0,nrow=900)), simplify=FALSE), simplify=FALSE), simplify=FALSE)

for (k in settings) {
    for (m in models) {
        cat(paste('setting: ', k, '; model: ', m, "\n", sep=""))
        for (j in 1:100) {
            cc = read.table(paste("coefs", m, k, j, "csv", sep='.'))
        
            for (v in vars) {
                output[[k]][[m]][[v]] = cbind(output[[k]][[m]][[v]], cc[[v]])
            }
        }
    }
}


#compute the B1 coefficient surfaces
N = 30 #number of width and length divisions in the domain
coord = seq(0, 1, length.out=N) #coordinates of the generated observations

if ((setting-1) %/% 6 == 0) {
    step = matrix(rep(ifelse(coord<=0.4, 0, ifelse(coord<0.6,5*(coord-0.4),1)), N), N, N)
} else if ((setting-1) %/% 6 == 1) {
    gradient = matrix(rep(coord, N), N, N)
} else if ((setting-1) %/% 6 == 2) {
    Xmat = matrix(rep(rep(coord, times=N), times=N), N**2, N**2)
    Ymat = matrix(rep(rep(coord, each=N), times=N), N**2, N**2)
    D = (Xmat-0.5)**2 + (Ymat-0.5)**2
    d = D[,435]
    parabola = matrix(max(d)-d, N, N)
    parabola = parabola / max(parabola)
}


wireframe(matrix(apply(output[[1]][['lagr']][['X1']], 1, function(x) mean(x==0)), 30, 30))
wireframe(matrix(apply(output[[1]][['lagr']][['X1']], 1, function(x) mean(x)), 30, 30))
wireframe(matrix(apply(output[[7]][['lagr']][['X1']], 1, function(x) mean(x)), 30, 30))


#Calculate the MISE for the first coefficient:
mise = matrix(0, nrow=0, ncol=3)
for (setting in 1:18) {
    #Get the proper coefficient surface
    if (setting < 7) {B1 = step}
    else if (setting < 13) {B1 = gradient}
    else {B1 = parabola}
    
    #compute the MISE
    row = sapply(c('lagr', 'gwr', 'oracle'), function(selection.method) mean((matrix(apply(output[[setting]][[selection.method]][['X1']],1,mean), 30, 30) - B1)**2))
    mise = rbind(mise, row)
}
rownames(mise) = NULL
round(mise, 4)



#Calculate the MISE for the other coefficients:
vars = c('X2', 'X3', 'X4', 'X5')
for (covariate in vars) {
    mise = matrix(0, nrow=0, ncol=3)
    for (setting in 1:18) {    
        #compute the MISE
        row = sapply(c('lagr', 'gwr', 'oracle'), function(selection.method) mean((matrix(apply(output[[setting]][[selection.method]][[covariate]],1,mean), 30, 30))**2))
        mise = rbind(mise, row)
    }
    rownames(mise) = NULL
    round(mise, 4)
}


#Calculate the frequency of an exact zero for the coefficients X2,...,X5:
vars = c('X2', 'X3', 'X4', 'X5')
for (covariate in vars) {
    pzero = matrix(0, nrow=0, ncol=3)
    for (setting in 1:18) {    
        #compute the proportion of zeros
        row = sapply(c('lagr', 'gwr', 'oracle'), function(selection.method) mean((matrix(apply(output[[setting]][[selection.method]][[covariate]],1,mean), 30, 30))**2))
        mise = rbind(pzero, row)
    }
    rownames(pzero) = NULL
    round(pzero, 4)
}


#Calculate the frequency of an exact zero for the coefficients X2,...,X5:
vars = c('X2', 'X3', 'X4', 'X5')
pz = matrix(0, nrow=0, ncol=3)
for (setting in 1:18) { 
    pzero = matrix(0, nrow=0, ncol=3)
    for (covariate in vars) {   
        #compute the proportion of zeros
        row = sapply(c('lagr', 'gwr', 'oracle'), function(selection.method) mean((matrix(apply(output[[setting]][[selection.method]][[covariate]],1,mean), 30, 30))**2))
        pzero = rbind(pzero, row)
    }
    pzero = colMeans(pzero)
    pz = rbind(pz, pzero)
}
rownames(pzero) = NULL
round(pzero, 4)