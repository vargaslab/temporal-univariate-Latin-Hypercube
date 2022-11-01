EyeModel_sim <- function (CoorX, CoorY, Variable, Variable_Sim, Direccion, Tol, NIntervalos, 
          Lags, Npares, Modelo, Nugget, SillYNugget, Alcance, MainTitle) 
{
  Sill <- SillYNugget - Nugget
  DatosGeo <- as.geodata(cbind(CoorX, CoorY, Variable), header = TRUE, 
                         coords.col = 1:2, data.col = 3)
  Nuevo <- variog(DatosGeo, breaks = c(seq(0, Lags * NIntervalos, 
                                           Lags)), trend = "cte", lambda = 1, estimator.type = "classical", 
                  nugget.tolerance = 0, direction = Direccion, tolerance = Tol, 
                  unit.angle = "degrees", pairs.min = Npares)
  DatosGeo_sim <- as.geodata(cbind(CoorX, CoorY, Variable_Sim), header = TRUE, 
                         coords.col = 1:2, data.col = 3)
  Nuevo_sim <- variog(DatosGeo_sim, breaks = c(seq(0, Lags * NIntervalos, 
                                           Lags)), trend = "cte", lambda = 1, estimator.type = "classical", 
                  nugget.tolerance = 0, direction = Direccion, tolerance = Tol, 
                  unit.angle = "degrees", pairs.min = Npares)
  Modelo1 <- Modelo
  if (Modelo == 1) {
    Modelo <- "exponential"
  }
  if (Modelo == 2) {
    Modelo <- "spherical"
  }
  if (Modelo == 3) {
    Modelo <- "gaussian"
  }
  if (Modelo1 == 1) {
    ModeloA <- "exponential"
  }
  if (Modelo1 == 2) {
    ModeloA <- "spherical"
  }
  if (Modelo1 == 3) {
    ModeloA <- "gaussian"
  }
  VMod <- (SillYNugget) - cov.spatial(Nuevo$u, cov.model = ModeloA, 
                                      cov.pars = c(Sill, Alcance))
  Error2 <- (Nuevo_sim$v - VMod)^2
  RMSE <- sqrt(sum(Error2)/length(VMod))
  
  MaxX <- max(Nuevo$u, Alcance) + 0.22 * (max(Nuevo$u, Alcance))
  MaxY <- max(Nuevo$v, Nuevo$var.mark, Sill + Nugget) + 0.22 * 
    (max(Nuevo$v, Nuevo$var.mark, Sill + Nugget))
  orden <- matrix(c(1, 2), ncol = 1, nrow = 2, byrow = T)
  div <- layout(orden, widths = 9, heights = c(3.5, 0.8))
  layout.show(div)
  par(mar = c(5, 5, 2, 2))
  plot(Nuevo$u, Nuevo$v, pch = 21, col = "black", bg = "black", 
       xlim = c(0, MaxX), ylim = c(0, MaxY), xlab = "Distance", 
       ylab = "Semivariance", cex.lab = 1.5, cex.axis = 1.2, 
       main = MainTitle)
  points(Nuevo_sim$u, Nuevo_sim$v, pch = 21, col = "green", bg = "green", 
       xlim = c(0, MaxX), ylim = c(0, MaxY), xlab = "Distance", 
       ylab = "Semivariance", cex.lab = 1.5, cex.axis = 1.2, 
       main = MainTitle)
  lines.variomodel(cov.model = Modelo, cov.pars = c(Sill, Alcance), 
                   nug = Nugget, col = "Black", lwd = 2, max.dist = Lags * 
                     NIntervalos)
  abline(h = Nuevo$var.mark, lty = 5, lwd = 2, col = "Red")
  abline(h = Sill + Nugget, lty = 4, lwd = 2, col = "darkviolet")
  abline(v = Alcance, lty = 3, lwd = 2, col = "Green")
  legend(0, MaxY, legend=c("Reference", "Simulated", Modelo, "Variance", "Sill", "Range"), 
         col = c("black", "green","black", "red", "darkviolet", "green"), lty = c(1,1,1, 
                                                                                 5, 4, 3), lwd = c(NA,NA,2,2,2,2), pch = c(19,19,NA,NA,NA,NA), bty = "n")
  ModEle <- cbind(Nugget, Sill + Nugget, Alcance, RMSE)
  rownames(ModEle) <- c(Modelo)
  colnames(ModEle) <- c("Nugget", "Sill", "Range", "RMSE")
  par(mar = c(0.5, 0.5, 0, 0.5))
  plot(0, 0, type = "n", xlim = c(0, 60), ylim = c(0, 7), xaxt = "n", 
       yaxt = "n", xlab = "", ylab = "")
  text(4, 5, labels = "Model", cex = 1.1)
  text(17, 5, labels = "Nugget", cex = 1.1)
  text(32, 5, labels = "Sill", cex = 1.1)
  text(47, 5, labels = "Range", cex = 1.1)
  text(58, 5, labels = "RMSE", cex = 1.1)
  text(4, 3, labels = Modelo, cex = 1.1)
  #  text(17, 3, labels = sprintf("%.4f", Nugget), cex = 1.1)
  #  text(32, 3, labels = sprintf("%.4f", Sill + Nugget), cex = 1.1)
  #  text(47, 3, labels = sprintf("%.4f", Alcance), cex = 1.1)
  #  text(58, 3, labels = sprintf("%.4f", RMSE), cex = 1.1)
  text(17, 3, labels = sprintf("%.9g", Nugget), cex = 1.1)
  text(32, 3, labels = sprintf("%.9g", Sill + Nugget), cex = 1.1)
  text(47, 3, labels = sprintf("%.9g", Alcance), cex = 1.1)
  text(58, 3, labels = sprintf("%f", RMSE), cex = 1.1)
  return(ModEle)
}
