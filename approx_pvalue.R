library(dplyr)

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
partition_seq <- function(S, max_iterations = 2000) {
    iterations <- min(length(S), max_iterations)
    
    result <- lapply(1:(log2(iterations) + 1), function(iter) {
        split(S, cut(seq_along(S),
                     # Intervalle halbieren
                     breaks = 2 ^ iter,
                     labels = FALSE))
    })
    
    return(result)
}


# Expected Euler Characteristic via Gaussian Kinematic Formular -----------

tGKF <- function(S, data, basisEC, level) {
    # Euler Charakteristik der Grundmenge
    chi <- basisEC
    
    # EC Dichten p_0 und p_1
    degree = ncol(data) - 1
    p0 <- 1 - pt(level, degree)
    p1 <- (1 + level ^ 2 / degree) ^ (-(degree - 1) / 2) / (2 * pi)
    
    # LKC_0
    L0 <- basisEC
    
    # LKC_1 gemaess Schaetzung
    residuals <- (data - rowSums(data) / ncol(data) ) / apply(data, 1, sd)
    
    df.residuals <- # Differenzenquotient
        ( (rbind(residuals[-1, ], residuals[length(residuals), ]) - residuals ) / diff(S)) %>%
        t() %>% apply(2, var) %>% sqrt()
    
    L1 <- # Trapezregel
        sum(diff(S) / 2 * (df.residuals[-length(df.residuals)] + df.residuals[-1]))
    
    # Bestimme das Ergebnis
    return(p0 * L0 + p1 * L1)
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
                    level){
    
    # Initialisierungen ----
    samplesize <- ncol(data)
    t.statistic <- - (rowSums(data) - level) / sqrt(samplesize) / apply(data, 1, sd)
    
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
        
        # Berechnung der p-Werte per Expected Euler Charakteristic
        basisEC <- 1    # EC von S = [0,1]
        pvalues <- sapply(1:N, function(i) {
            tGKF(sortedUnions[[i]],
                 data[S %in% sortedUnions[[i]],],
                 basisEC,
                 teststatistics[i])
        })
        
        # Akzeptiere V_k, wenn p(x_k, V_k) >= alpha
        tempU <- c()
        k <- min(which(pvalues >= alpha))
        if (k <= N) tempU <- sortedUnions[[k]]
        
        # tempU zu U hinzufuegen
        U <- union(U, tempU)
    }
    
    return(U)
}