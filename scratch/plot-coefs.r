library(RandomFields)
library(brooks)
library(dplyr)

N = 20 #number of width and length divisions in the domain
coord = seq(0, 1, length.out=N)

#Calculate the coefficient surfaces
B1 = RFsimulate(RMexp(var=2.5, scale=1), x=coord, y=coord)@data[[1]] %>% matrix(20,20)
B2 = RFsimulate(RMexp(var=0.5, scale=1), x=coord, y=coord)@data[[1]] %>% matrix(20,20)
B3 = RFsimulate(RMexp(var=0.1, scale=1), x=coord, y=coord)@data[[1]] %>% matrix(20,20)
B4 = RFsimulate(RMexp(var=0.02, scale=1), x=coord, y=coord)@data[[1]] %>% matrix(20,20)

pdf("~/Desktop/coefs.pdf", 9, 2.88)
par(oma=c(1,1,1,1))
par('mar'=c(4,1,1,1))
layout(matrix(1:4,1,4))

matplot(B1 %>% matrix(20,20), border=NA, show.legend=TRUE, axes=FALSE, xlab=expression(beta[1]), ylab=NA, cex.lab=2)
matplot(B2 %>% matrix(20,20), border=NA, show.legend=TRUE, axes=FALSE, xlab=expression(beta[2]), ylab=NA, cex.lab=2)
matplot(B3 %>% matrix(20,20), border=NA, show.legend=TRUE, axes=FALSE, xlab=expression(beta[3]), ylab=NA, cex.lab=2)
matplot(B4 %>% matrix(20,20), border=NA, show.legend=TRUE, axes=FALSE, xlab=expression(beta[4]), ylab=NA, cex.lab=2)

dev.off()
