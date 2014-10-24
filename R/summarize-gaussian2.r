library(dplyr)

MISEX = list('lagr'=list(), 'gwr'=list())
MISEY = list('lagr'=list(), 'gwr'=list())
freq.zero = list()

for (b in 1:5) {
for (i in 1:3) {
for (process in 1:6) {
#for (j in 1:18) {
    j = (b-1)*3*6 + (i-1)*6 + process
    truth = read.table(paste("output/truth", j, "csv", sep="."))
    coef.gwr = read.table(paste("output/coefs.gwr", j, "csv", sep="."))
    coef.lagr = read.table(paste("output/coefs.lagr", j, "csv", sep="."))

    results.gwr = read.table(paste("output/results.gwr", j, "csv", sep="."))
    results.lagr = read.table(paste("output/results.lagr", j, "csv", sep="."))

    MISEX[['lagr']][[j]] = list()
    MISEX[['gwr']][[j]] = list()
    freq.zero[[j]] = list()

    for (k in 1:3) {
        B = paste("B", k, sep="")
        X = paste("X", k, sep="")

        MISEX[['lagr']][[j]][[k]] = (truth[[B]] - coef.lagr[[X]])**2 %>% mean
        MISEX[['gwr']][[j]][[k]] = (truth[[B]] - coef.gwr[[X]])**2 %>% mean

        freq.zero[[j]][[k]] = (coef.lagr[[X]] == 0) %>% mean
    }
    MISEX[['lagr']][[j]][[4]] = (coef.lagr$X4)**2 %>% mean
    MISEX[['gwr']][[j]][[4]] = (coef.gwr$X4)**2 %>% mean
    freq.zero[[j]][[4]] = (coef.lagr$X4 == 0) %>% mean

    MISEY[['lagr']][[j]] = (results.lagr$residual)**2 %>% mean
    MISEY[['gwr']][[j]] = (results.gwr$residual)**2 %>% mean
}
}
}
