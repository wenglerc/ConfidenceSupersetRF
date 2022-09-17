library(dplyr)
library(mvtnorm)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(microbenchmark)

theme_set(theme_light())

# Generieren eines Zufallsprozesses X(s) = mu(s) + sigma(s)*epsilon(s) --------

mu <- function(S) {
    sin(S * 2 * pi) * 3
} # Erwartungswertfunktion

sigma <- function(S){
    rep(1, length(S))
} # Varianzfunktion

eps <- function(S) {
    outer(S, S, FUN = function(x, y) exp( -( x - y )^2 / 0.0005 ))
} # Fehlerprozess

linearModell <- function(n, S, mu, eps) {
    mu(S) + sigma(S) * t(rmvnorm(n, sigma = eps(S)))
} # Lineares Modell


# HAUPTTEIL: Beispielhafter Durchlauf aller Funktionen --------------------

#### Initialisierung globaler Variablen -----------------------------------

randomseed <- 1
samplesize.max <- 200
iterations <- 1000
gridpoints <- 200               # Anzahl equidistanter Gitterpunkte

alpha <- 0.05
level <- 0                      # Grenzwert aus der Null-Hypothese
S <- seq(0, 1, length.out = gridpoints)  # Grundmenge S = [0,1]
partitions <- partition_seq(S, pieces = 2)  # Folge von Partitionen (degenerierend, feiner werdend)


#### Testfaelle generieren --------------------------------------------------

set.seed(randomseed)

# Realisierungen von X
data.all <- lapply(1:iterations, function(i) {
    data.frame(linearModell(samplesize.max, S, mu, eps))
})
data.mu <- data.frame(S = S, Erwartungswert = mu(S))

# Darstellung von Realisierungen
data.plot <- melt(
    data.frame(S = S, data.all[[1]])[, 1:10],
    id.vars = "S",
    variable.name = 'Sample',
    value.name = "Wert"
)
p <- ggplot() +
    labs(x = "Grundmenge S", y = "", colour = "Realisierung") +
    theme(legend.position = "none")
p + geom_line(data = data.plot, aes(S, Wert, colour = Sample), alpha = 0.5) +
    geom_line(data = data.mu, aes(S, Erwartungswert), colour = "red")

#### Konstruktion --------------------------------------------------

# U <- ConfSet(data, S, partitions, alpha, level)
# cat("Genauigkeit: ", length(U[!(U %in% S[mu(S) > level])]) / length(U))

# mbm <- microbenchmark(
#     V1 = ConfSet(data, S, partitions, alpha, level),
#     V2 = ConfSet2(data, S, partitions, alpha, level),
#     times = 10
# )
# 
# autoplot(mbm)

samplesize.list <- c(10,20,50,100,150,200)

results <- sapply(samplesize.list, function(N){
    results.U <- lapply(data.all, function(data) ConfSet(data[,1:N], S, partitions, alpha, level))
    fullcovered <- sapply(results.U, function(U) all(S[data.mu$Erwartungswert > level] %in% U))
    
    covering.rate <- length(fullcovered[fullcovered == T]) / length(fullcovered)
})



