require(dplyr)
require(mvtnorm)
require(tidyverse)
require(reshape2)
require(ggplot2)
require(doParallel)
require(microbenchmark)

library(tictoc)

theme_set(theme_light())

# Generieren eines Zufallsprozesses X(s) = mu(s) + sigma(s)*epsilon(s) --------

mu <- function(S) {
    sin(S * 5 * pi) * 5 * S
} # Erwartungswertfunktion

sigma <- function(S){
    rep(1, length(S))
} # Varianzfunktion

eps <- function(S) {
    outer(S, S, FUN = function(x, y) exp( -( x - y )^2 / 0.0005 ))
} # Fehlerprozess

linearModell <- function(n, S, mu, eps) {
    mu(S) + sigma(S) * t(mvtnorm::rmvnorm(n, sigma = eps(S)))
} # Lineares Modell


# HAUPTTEIL: Beispielhafter Durchlauf aller Funktionen --------------------

#### Initialisierung globaler Variablen -----------------------------------

randomseed <- 1
samplesize.max <- 100
iterations <- 100
gridpoints <- 200               # Anzahl equidistanter Gitterpunkte

alpha <- 0.05
level <- 0                      # Grenzwert aus der Null-Hypothese
S <- seq(0, 1, length.out = gridpoints)  # Grundmenge S = [0,1]
partitions <- partition_seq(S, pieces = 3)  # Folge von Partitionen (degenerierend, feiner werdend)


#### Testfaelle generieren --------------------------------------------------

set.seed(randomseed)

# Realisierungen von X
cluster <- makeCluster(detectCores() - 1)
registerDoParallel(cluster)
data.all <- foreach(i = 1:iterations) %dopar% {
    data.frame(linearModell(samplesize.max, S, mu, eps))
}
stopCluster(cluster)

data.mu <- data.frame(S = S, Erwartungswert = mu(S))
S0 = S[mu(S) > level]

# Darstellung von Realisierungen
data.plot <- reshape2::melt(
    data.frame(S = S, data.all[[1]])[, 1:10],
    id.vars = "S",
    variable.name = 'Sample',
    value.name = "Wert"
)
p.ex <- ggplot() +
    labs(x = "Grundmenge S", y = "", colour = "Realisierung") +
    theme(legend.position = "none")
p.ex + geom_line(data = data.plot, aes(S, Wert, colour = Sample), alpha = 0.5) +
    geom_line(data = data.mu, aes(S, Erwartungswert), colour = "red")

#### Konstruktion --------------------------------------------------

# U <- ConfSet(data.all[[3]], S, partitions, alpha, level, pmethod = "mboot")
# cat("S0 nicht in U: ", length( S0[ !(S0 %in% U)] ) / length(S0),
#     "\nU nicht in S0: ", length( U[ !(U %in% S0)] ) / length(U))
# 
# U <- ConfSet2(data.all[[3]], S, partitions, alpha, level, pmethod = "mboot")
# cat("S0 nicht in U: ", length( S0[ !(S0 %in% U)] ) / length(S0),
#     "\nU nicht in S0: ", length( U[ !(U %in% S0)] ) / length(U))

mbm <- microbenchmark(
    V1 = ConfSet(data.all[[1]], S, partitions, alpha, level, pmethod = "mboot"),
    V2 = ConfSet2(data.all[[1]], S, partitions, alpha, level, pmethod = "mboot"),
    times = 5
)
autoplot(mbm)

samplesize.list <- c(20,50,100,200,500)

tic("GKF")
results.gkf <- sapply(samplesize.list, function(N) {
    cluster <- makeCluster(detectCores() - 1)
    registerDoParallel(cluster)
    results.U <-
        foreach(data = data.all, .export = ls(globalenv()), .packages = c("dplyr", "stats")) %dopar% {
            ConfSet(data[, 1:N], S, partitions, alpha, level, pmethod = "tgkf")
        }
    stopCluster(cluster)

    fullcovered <-
        sapply(results.U, function(U)
            all(S[data.mu$Erwartungswert > level] %in% U))
    covering.rate <-
        length(fullcovered[fullcovered == T]) / length(fullcovered)

})
toc()

tic("MBoot")
results.mb <- sapply(samplesize.list, function(N) {
    cluster <- makeCluster(detectCores() - 1)
    registerDoParallel(cluster)
    results.U <-
        foreach(data = data.all, .export = ls(globalenv()), .packages = c("dplyr", "stats")) %dopar% {
            ConfSet(data[, 1:N], S, partitions, alpha, level, pmethod = "mboot")
        }
    stopCluster(cluster)
    
    fullcovered <-
        sapply(results.U, function(U)
            all(S[data.mu$Erwartungswert > level] %in% U))
    covering.rate <-
        length(fullcovered[fullcovered == T]) / length(fullcovered)
    
})
toc()

#### Auswertung --------------------------------------------------

df.gkf <- data.frame(samplesize = samplesize.list,
                     res = results.gkf, method = "t-GKF")
# df.mb <- data.frame(samplesize = samplesize.list,
#                     res = results.mb, method = "Multiplier Bootstrap")
# data.plot <- rbind(df.gkf, df.mb)

ggplot(data = df.gkf, aes(x = samplesize, y = res, color = method)) +
    labs(x = "Anzahl an Samples", y = "Überdeckungsrate", color = "Methode") +
    theme(legend.position = "top") +
    geom_line() + geom_point() +
    geom_hline(yintercept = 1 - alpha, linetype = 'dashed', col = 'red')



