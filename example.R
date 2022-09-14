library(dplyr)
library(mvtnorm)
library(tidyverse)
library(reshape2)
library(ggplot2)

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
samplesize <- 50
gridwidth <- 1e-3               # Schrittweite der Diskretisierung der Grundmenge

alpha <- 0.05
level <- 0                      # Grenzwert aus der Null-Hypothese
S <- seq(0, 1, by = gridwidth)  # Grundmenge S = [0,1]
partitions <- partition_seq(S)  # Folge von Partitionen (degenerierend, feiner werdend)


#### Testfaelle generieren --------------------------------------------------

set.seed(randomseed)

# Realisierungen von X
data <- data.frame(linearModell(samplesize, S, mu, eps))
data.mu <- data.frame(S = S, Erwartungswert = mu(S))
data.standardized <- data.frame(
    S = S,
    Mittelwert.normiert = rowSums(data[,-1]) / sqrt(samplesize) / apply(data[,-1], 1, sd)
    )

# Darstellung einiger Realisierungen
data.plot <- melt(data.frame(S=S, data)[, 1:10], id.vars = "S", variable.name = 'Sample', value.name = "Wert")
p <- ggplot() +
    labs(x ="Grundmenge S", y = "", colour = "Realisierung") +
    theme(legend.position="none")
p + geom_line(data = data.plot, aes(S, Wert, colour = Sample), alpha = 0.5) +
    geom_line(data = data.mu, aes(S, Erwartungswert), colour = "red")

