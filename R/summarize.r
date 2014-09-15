library(dplyr)

MISEX = list('lagr'=list(), 'gwr'=list())
MISEY = list('lagr'=list(), 'gwr'=list())
freq.zero = list()

for (j in 1:18) {
    truth = read.table(paste("output/truth", j, "csv", sep="."))
    coef.gwr = read.table(paste("output/coefs.gwr", j, "csv", sep="."))
    coef.lagr = read.table(paste("output/coefs.lagr", j, "csv", sep="."))

    results.gwr = read.table(paste("output/results.gwr", j, "csv", sep="."))
    results.lagr = read.table(paste("output/results.lagr", j, "csv", sep="."))

    MISEX[['lagr']][[j]] = list()
    MISEX[['gwr']][[j]] = list()
    freq.zero[[j]] = list()

    for (k in 1:4) {
        B = paste("B", k, sep="")
        X = paste("X", k, sep="")

        MISEX[['lagr']][[j]][[k]] = (truth[[B]] - coef.lagr[[X]])**2 %>% mean
        MISEX[['gwr']][[j]][[k]] = (truth[[B]] - coef.gwr[[X]])**2 %>% mean

        freq.zero[[j]][[k]] = (coef.lagr[[X]] == 0) %>% mean
    }
    MISEX[['lagr']][[j]][[5]] = (coef.lagr$X5)**2 %>% mean
    MISEX[['gwr']][[j]][[5]] = (coef.gwr$X5)**2 %>% mean
    freq.zero[[j]][[5]] = (coef.lagr$X5 == 0) %>% mean

    MISEY[['lagr']][[j]] = (results.lagr$residual)**2 %>% mean
    MISEY[['gwr']][[j]] = (results.gwr$residual)**2 %>% mean
}