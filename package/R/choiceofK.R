choiceofK <- function(AIC, BIC, RSS, K, filename) {
    Kmax <- length(AIC)
    pdf(filename, width = 13, height = 4)
    par(mfrow = c(1, 3))
    plot(K, RSS, type = "b", xlab = "Number of latent variables")
    plot(K, AIC, type = "b", xlab = "Number of latent variables")
    plot(K, BIC, type = "b", xlab = "Number of latent variables")
    dev.off()
}