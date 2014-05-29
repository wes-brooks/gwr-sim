m = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, varselect.method='AICc', bw=0.18, kernel=epanechnikov, bw.type='knn', verbose=FALSE)

#Use a profile likelihood to determine model-averaged coefficients for the LAGR model and to find how likely is each local coefficient to be zero:
zz = sapply(m[['model']][['models']], function(y) apply(y[['model']][['beta']], 1, function(x) sum((x==0)*w)/sum(w)))[1:5,]
cc = sapply(m[['model']][['models']], function(y) apply(y[['model']][['beta']], 1, function(x) sum(x*w)/sum(w)))[1:5,]