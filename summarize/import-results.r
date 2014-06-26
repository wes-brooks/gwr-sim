vars = c("X.Intercept.", "X1", "X2", "X3", "X4", "X5")
methods = c("lagr", "oracle", "gwr")
settings = 1:18

setwd("~/Dropbox/gwr-sim")

output = sapply(settings, function(y)
    sapply(methods, function(t) 
        sapply(vars, function(x)
            return(matrix(ncol=0,nrow=900)), simplify=FALSE), simplify=FALSE), simplify=FALSE)

for (k in settings) {
    for (m in methods) {
        cat(paste('setting: ', k, '; model: ', m, "\n", sep=""))
        for (j in 1:100) {
            cc = read.table(paste("coefs", m, k, j, "csv", sep='.'))
        
            for (v in vars) {
                output[[k]][[m]][[v]] = cbind(output[[k]][[m]][[v]], cc[[v]])
            }
        }
    }
}


vars2 = c("Y", "X1", "X2", "X3", "X4", "X5")
simdata = sapply(settings, function(y)
    sapply(vars2, function(x)
        return(matrix(ncol=0,nrow=900)), simplify=FALSE), simplify=FALSE)

for (k in settings) {
    cat(paste('setting: ', k, "\n", sep=""))
    for (j in 1:100) {
        cc = read.table(paste("simdata", k, j, "csv", sep='.'))
        
        for (v in vars2) {
            simdata[[k]][[v]] = cbind(simdata[[k]][[v]], cc[[v]])
        }
    }
}



#compute the B1 coefficient surfaces
N = 30 #number of width and length divisions in the domain
coord = seq(0, 1, length.out=N) #coordinates of the generated observations

step = matrix(rep(ifelse(coord<=0.4, 0, ifelse(coord<0.6,5*(coord-0.4),1)), N), N, N)
gradient = matrix(rep(coord, N), N, N)
parabola = 1 - ((rep(coord, times=N)-0.5)**2 + (rep(coord, each=N)-0.5)**2)/0.5


#Calculate the MISE for the first coefficient:
mise = matrix(0, nrow=0, ncol=3)
for (setting in 1:18) {
    #Get the proper coefficient surface
    if (setting < 7) {B1 = step}
    else if (setting < 13) {B1 = gradient}
    else {B1 = parabola}
    
    #compute the MISE
    row = sapply(c('lagr', 'gwr', 'oracle'), 
                 function(selection.method) {
                     mean(sweep(output[[setting]][[selection.method]][['X1']], 1, as.vector(B1))**2)
                 })
    mise = rbind(mise, row)
}
rownames(mise) = NULL
round(mise, 4)
xtable(mise)



#Calculate the MISE for the response:
misey = matrix(0, nrow=0, ncol=3)
for (setting in 1:18) {    
    #compute the MISE
    row = sapply(c('lagr', 'gwr', 'oracle'), 
         function(selection.method) {
             B1 = output[[setting]][[selection.method]][['X1']]
             X1 = simdata[[setting]][['X1']]
             intercept = output[[setting]][[selection.method]][['X.Intercept.']]
             truth = simdata[[setting]][['Y']]
             err = matrix(NA, nrow=900, ncol=0)
             
             for (j in 1:100) {
                 err = cbind(err, truth[,j] - X1[,j]*B1[,j] - intercept[,j])
             }
             
             return(mean(err**2))
         })
    misey = rbind(misey, row)
}
rownames(misey) = NULL
round(misey, 4)
xtable(misey)


bold.misey = matrix(FALSE, nrow=18, ncol=3)
ital.misey = matrix(FALSE, nrow=18, ncol=3)
std = matrix(rbind(c(0.25, 0.25, 0.25), c(1, 1, 1)), nrow=18, ncol=3)

bold.misey = matrix(FALSE, nrow=nrow(mise), ncol=ncol(mise))
for (i in 1:nrow(mise)){
    bold.misey[i,which.min(abs(misey[i,] - std[i,]))] = TRUE
}

#Which entries should be italicised (for having the second-lowest MISE)?
ital.misey = matrix(FALSE, nrow=nrow(mise), ncol=ncol(mise))
for (i in 1:nrow(mise)){
    ital.misey[i,order(abs(misey[i,] - std[i,]))[2]] = TRUE
}

xtable.printbold(xtable(misey), which.bold=bold.misey, which.ital=ital.misey)


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
pzero = sapply(vars, function(t) return(vector()), simplify=FALSE)
for (covariate in vars) {
    for (setting in 1:18) {    
        #compute the proportion of zeros
        z = mean(apply(output[[setting]][['lagr']][[covariate]],1,'==',0)**2)
        pzero[[covariate]] = c(pzero[[covariate]], z)
    }
}

zz = vector()
for (i in 1:18) {
    zz = c(zz, mean(sapply(pzero, function(x) x[i])))
}

