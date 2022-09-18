require(dplyr)
require(stats)
require(pracma)
require(doParallel)

library(tictoc)

# Erstellung einer Partitionsfolge (degenerierend, verfeinernd) -----------

#' Erzeugt eine degenerative Folge von Partitionen, die immer feiner werden.
#'
#' @param S Grundmenge, zu der eine Folge von Partitionen erstellt werden soll
#' @param max_iterations Maximale Anzahl an Folgegliedern
#'
#' @return Eine Liste mit Listen, die immer feiner werdende Aufspaltungen
#' der Menge S enthalten 
#' @export
#'
#' @examples
partition_seq <- function(S, pieces = 2, max_iterations = 2000) {
    iterations <- min(length(S), max_iterations)
    
    result <- lapply(1:(log(iterations, base = pieces) + 1), function(iter) {
        split(S, cut(seq_along(S),
                     # Intervalle in pieces Teile teilen
                     breaks = pieces ^ iter,
                     labels = FALSE))
    })
    
    return(result)
}


# Expected Euler Characteristic via Gaussian Kinematic Formular -----------

tGKF <- function(S, data, level, basisEC = 1) {
    
    # EC Dichten p_0 und p_1
    degree = ncol(data) - 1
    p0 <- 1 - pt(level, degree)
    p1 <- (2 *pi)^(-1) * (1 + level ^ 2 / degree)^(-(degree - 1) / 2)
    
    # LKC_0
    L0 <- basisEC
    
    # LKC_1 gemaess Schaetzung
    residuals <-
        (data - rowSums(data) / ncol(data)) / apply(data, 1, sd)
    
    if (length(S) > 1) {
        df.residuals <- apply(residuals, 2, function(values) {
            func <- stats::splinefun(S, values, method = "natural")
            pracma::fderiv(func, S, method = "central")
        }) %>% apply(1, stats::sd)
        L1 <- # Trapezregel
            sum(diff(S) / 2 * (df.residuals[-length(df.residuals)] + df.residuals[-1]))
    } else {
        L1 <- 0
    }
    
    # Bestimme die EEC
    p0 * L0 + p1 * L1
}

# Multiplier Bootstrap -----------
mBoot <- function(S, data, level, iter = 100) {
    # Residuen
    samplesize = ncol(data)
    residuals <-
        sqrt(samplesize / (samplesize - 1)) * (data - rowSums(data) / samplesize)
    
    # Standardnormalverteilte Multiplier
    g <- rnorm(samplesize * iter) %>% matrix(samplesize, iter)
    
    # Approximation der bedingten Wahrscheinlichkeit der Summe gewichteter Residuen
    cluster <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cluster)
    
    mboot.distr <- foreach::foreach(i = 1:iter, .combine = "c") %dopar% {
        # Standardabweichung nach Anwendung der Multiplier
        sd.gX <- apply(t(t(data) * g[,i]), 1, sd)
        # Maximum der Summen gewichteter Residuen über der Grundmenge
        max(sqrt(samplesize) ^ (-1) * colSums(t(residuals / sd.gX) * g[,i]))
    }
    parallel::stopCluster(cluster)
    
    # TODO: Richtigkeit checken
    length(mboot.distr[mboot.distr >= level]) / length(mboot.distr)
    
}


# Berechnung der Konfidenz-Obermenge auf Basis einer Partitionsfolge ------

#' Berechnung einer Konfidenz-Obermenge des Excursion Sets
#' S_0 = \eqn{\{\mu(s) > level\}} bzgl. eines Zufallsfelds
#' \eqn{X_s(omega) = mu(s) + eps_s(omega)} zu einer kompakten Grundmenge S mit
#' Eulercharakteristik 1.
#' Die Berechnung basiert auf der Metode aus
#' "False Discovery Control for Random Fields" (2004) - M. Perone Pacifico,
#' C. Genovese, I. Verdinelli and L. Wasserman
#'
#' @param df Realisierungen des Zufallsfeldes
#' @param S Grundmenge S (1D und EC(S) = 1)
#' @param partitions Partitionierung der Grundmenge (degenerierend und verfeinernd)
#' @param alpha Signifikanzniveau für die Konfidenzmenge
#' @param level Grenzwert des Excursion Sets
#'
#' @return Konfidenzmenge U mit \eqn{P(S_0 \subset U) \geq 1 - \alpha}
#' @export
#'
#' @examples
ConfSet <- function(data,
                    S,
                    partitions,
                    alpha,
                    level,
                    pmethod = "tgkf"){
    
    # Initialisierungen ----
    samplesize <- ncol(data)
    t.statistic <- - sqrt(samplesize) * (rowSums(data)/samplesize - level) / apply(data, 1, sd)
    
    U <- c()

    # Aussere Schleife: Loop über die Folge an Partitionen ----
    for (partS in partitions) {
        # Reduktion der aktuell betrachteten Partition auf jene Teilmengen,
        # die zuvor noch nicht zu U gezaehlt wurden: partS^*_n = partS_n \ U_{n-1}
        partS <- lapply(partS, function(vec) if (!all(vec %in% U)) vec) %>%
            Filter(Negate(is.null), .)
        if (all(sapply(partS, is.null))) break
        
        N <- length(partS)
        # Innere Schleife: Berechnung von U_n nach Partition partS^*_n ----
        
        # Berechnung und Sortierung der realisierten Werte der Teststatistik
        teststatistics <- lapply(partS, function(vec) t.statistic[S %in% vec] %>% max()) %>%
            unlist()
        sortorder <- order(teststatistics, decreasing = T)
        teststatistics <- sort(teststatistics, decreasing = T)
        
        # Vereinigungen V_k ab dem k-ten Partitionselement (geordnet wie Teststatisik)
        sortedUnions <-
            lapply(1:N, function(i) {
                partS[sortorder[i:N]] %>% unlist() %>% unname() %>% sort()
            })

        # Berechnung der p-Werte
        accept.index = 1
        for (i in N:1) {
            if (pmethod == "tgkf") {
                pvalue <- tGKF(sortedUnions[[i]],
                               data[S %in% sortedUnions[[i]], ],
                               teststatistics[i],
                               basisEC = 1) %>% unname()
            } else if (pmethod == "mboot") {
                pvalue <- mBoot(sortedUnions[[i]],
                                data[S %in% sortedUnions[[i]], ],
                                teststatistics[i]) %>% unname()
            }
            
            # Akzeptiere V_k, wenn p(x_k, V_k) >= alpha
            if (pvalue < alpha) {
                accept.index <- ifelse(i == N, 0, i + 1)
                break
            }
        }

        # Erweitere U_n = U_{n-1} + V_k
        if (accept.index > 0) U <- sort(union(U, sortedUnions[[accept.index]]))

    }

    return(U)
}
