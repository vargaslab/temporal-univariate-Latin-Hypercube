#--------------------------------------------------------------#
#                   Variography Analysis                       #
#--------------------------------------------------------------#

root_dir<-getwd()

# Creates a folder to store results for Variographic Analysis (VA)
dir.create(paste(getwd(),"/Results/Time_Series/VA", sep=""))

VA_dir<-paste(root_dir,"/Results/Time_Series/VA",sep="")

#--------------------------------------------------------------#
#                    Temporal Distribution                      #
#--------------------------------------------------------------#

# input
log_list=list(CO2, CH4, N2O)
log_xlabel_list=list(CO2label, CH4label, N2Olabel)
n_logs= 3
plot_mean=TRUE
plot_median=TRUE
fontproportion=1.0
png(paste(VA_dir,"/CO2_CH4_N2O_Temporal_Distr.png",sep=""), 
    bg = "white", width = 1000, height = 1000, res = 100)
par(mfrow = c(n_logs, 1), mar = c(4, 6, 2, 2), cex.lab= 1.5)
p=unlist(log_list[1])
log_label=unlist(log_xlabel_list[1])
s=summary(p)
plot(Time, p,  type = "p", 
     ylab = log_label , xlab = "", bty="o")
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Time, p, type = "p", 
     ylab =log_label , xlab = "", bty="o" )
if (plot_mean) {
  abline(h = s[4], col = "Red")
}  
if (plot_median) {
  abline(h = s[3], col = "Blue", lty = 2)
}
# legend(50, 20, legend=c("Mean", "Median"),
#        col=c("red", "blue"), lty=1:2, box.lty=0)
# par(adj=0)  # Para justificar a izquierda
# text(x=10, y=23, '(a)', cex=4)

for (i in 2:(n_logs-1)) {
  p=unlist(log_list[i])
  log_label=unlist(log_xlabel_list[i])
  s=summary(p)
  plot(Time, p,  type = "p", main = "",
       ylab = log_label , xlab = "", bty="o")
  grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
  par(new=TRUE)
  plot(Time, p, type = "p", 
       ylab =log_label , xlab = "", bty="o")
  if (plot_mean) {
    abline(h = s[4], col = "Red")
  }
  if (plot_median) {
    abline(h = s[3], col = "Blue", lty = 2)
  }
}
# par(adj=0)  # Para justificar a izquierda
# text(x=10, y=0.0005, '(b)', cex=4)

par(adj=0.5)
for (i in n_logs) {
  p=unlist(log_list[i])
  log_label=unlist(log_xlabel_list[i])
  s=summary(p)
  plot(Time, p,  type = "p", main = "",
       ylab = log_label , xlab = "Time (days)", bty="o" )
  grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
  par(new=TRUE)
  plot(Time, p, type = "p", 
       ylab =log_label , xlab = "Time (days)", bty="o")
  if (plot_mean) {
    abline(h = s[4], col = "Red")
  }
  if (plot_median) {
    abline(h = s[3], col = "Blue", lty = 2)
  }
}
# legend(200, 15, legend=c("Mean", "Median"),
#        col=c("red", "blue"), lty=1:2, box.lty=0)

# par(adj=0)  # Para justificar a izquierda
# text(x=10, y=0.032, '(c)', cex=4)
dev.off()


#--------------------------------------------------------------#
#               Trend (Stationarity) Analysis                  #
#--------------------------------------------------------------#

t = Time
# Estimation of the experimental variogram
N_lags<- 10 # 250 #10 #length(t)/2.5
lag_value <- 16 #25 #12 # min(dist(t)) # or delta t
ti = numeric(length(t))
ti_Stat <- Estadisticas(ti)
png(paste(VA_dir,"/CO2_VarioEstimation_10lags.png",sep=""), bg = "white", width = 1500, height = 1000, res = 150)
CO2_VarioEstimation<-Variograma(t, ti, 
                                         CO2, 0, 90, N_lags, lag_value, 1, MainTitle="Empirical variogram", xlab= "Time (days)", ylab = "Semivariogram") # CO2 temporal variogram
dev.off()

N_lags<- 10 # 250 #10 #length(t)/2.5
lag_value <-  16 #25 #12 # min(dist(t)) # or delta t
png(paste(VA_dir,"/CH4_VarioEstimation_10lags.png",sep=""), bg = "white", width = 1000, height = 1000, res = 150)
CH4_VarioEstimation<-Variograma(t, ti, 
                                CH4, 0, 90, N_lags, lag_value, 1, "", "Time (days)") # CO2 temporal variogram
dev.off()

N_lags<- 10 # 250 #10 #length(t)/2.5
lag_value <-  16 #25 #12 # min(dist(t)) # or delta t
png(paste(VA_dir,"/N2O_VarioEstimation_10lags.png",sep=""), bg = "white", width = 1000, height = 1000, res = 150)
N2O_VarioEstimation<-Variograma(t, ti, 
                                N2O, 0, 90, N_lags, lag_value, 1, "", "Time (days)") # CO2 temporal variogram
dev.off()


# Manual Variogram Model Fitting
#variogram model (1- exponential, 2- spherical, 3- gaussian)
CO2_vario_model<- 3
CO2_nugget<-4.0
CO2_sill_and_nugget<- 28
CO2_rank <- 80
png(paste(VA_dir,"/CO2_VarioEyeEstimation.png",sep=""),   bg = "white", width = 1000, height = 1000, res = 150)
CO2_EyeModelVarioFit<-EyeModel(t, ti, 
                                           CO2, 0, 90, N_lags, lag_value, 1, 
                                    CO2_vario_model, CO2_nugget, CO2_sill_and_nugget, CO2_rank,
                                    "") # A fitted temporal variogram model
dev.off()


CH4_vario_model<- 2
CH4_nugget<-0.07
CH4_sill_and_nugget<- 0.15
CH4_rank <- 110
png(paste(VA_dir,"/CH4_VarioEyeEstimation.png",sep=""),   bg = "white", width = 1000, height = 1000, res = 150)
CH4_EyeModelVarioFit<-EyeModel(t, ti, 
                               CH4, 0, 90, N_lags, lag_value, 1, 
                               CH4_vario_model, CH4_nugget, CH4_sill_and_nugget, CH4_rank,
                               "") # A fitted temporal variogram model
dev.off()

N2O_vario_model<- 2
N2O_nugget<-2.6
N2O_sill_and_nugget<- 2.6
N2O_rank <- 0
png(paste(VA_dir,"/N2O_VarioEyeEstimation.png",sep=""),   bg = "white", width = 1000, height = 1000, res = 150)
N2O_EyeModelVarioFit<-EyeModel(t, ti, 
                               N2O, 0, 90, N_lags, lag_value, 1, 
                               N2O_vario_model, N2O_nugget, N2O_sill_and_nugget, N2O_rank,
                               "") # A fitted temporal variogram model
dev.off()
