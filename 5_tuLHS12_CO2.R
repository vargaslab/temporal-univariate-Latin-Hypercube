#----------------------Optimization--------------------------------#
dir.create(paste(getwd(),"/Results/Time_Series/Optimization", sep=""))
ee_dir<-paste(root_dir,"/Results/Time_Series/Optimization",sep="")


#-----------------------------------------Sampling_12_CO2--------------------------#
n_sub <- 12
Fn_optimal_subsampling <- seq(1/n_sub,1,length.out = n_sub)
# subsampling by inverse function
CO2_optimal_subsampling <- sapply(Fn_optimal_subsampling, Fn.inv.Bernshtein, CO2)
CO2_optimal_sub <- quantile(CO2,Fn_optimal_subsampling)
plot(CO2_optimal_subsampling,CO2_optimal_sub)
#Temp_optimal_subsampling <- sapply(Fn_optimal_subsampling, Fn.inv.Bernshtein, Temp)
optimal_time_position_CO2 <- NULL
#optimal_time_position_Temp <- NULL
for (i in 1:n_sub) {
  optimal_time_position_CO2 <- rbind(optimal_time_position_CO2, which.min(abs(CO2-CO2_optimal_subsampling[i])))
  #optimal_time_position_Temp <- rbind(optimal_time_position_Temp, which.min(abs(Temp-Temp_optimal_subsampling[i])))
}
(CO2_optimal_subsample <- CO2[optimal_time_position_CO2])
CO2_optimal_subsampling
CO2_optimal_subsample_Stat <- Estadisticas(CO2_optimal_subsample)
Fn_sub_opt <- rank(CO2_optimal_subsample)/length(CO2_optimal_subsample)
(sum(abs(Fn_sub_opt-Fn_opt_sub))) # 0.5


#--------------------temporal sub-sampling defined by user---------------------#
# Temporal subsampling
CO2_temporal_subsample <- CO2[seq(from= 1,to =length(CO2), by= (length(CO2)/n_sub))]
Fn_temporal_subsampling <- rank(CO2_temporal_subsample)/length(CO2_temporal_subsample)
CO2_temporal_subsample_Stat <- Estadisticas(CO2_temporal_subsample)

Fn_CO2_temporal_subsample <- rank(CO2_temporal_subsample)/length(CO2_temporal_subsample)
(sum(abs(Fn_CO2_temporal_subsample-Fn_optimal_subsampling))) # 0.5
t_temporal_sub <- Time[seq(from= 1,to =length(CO2), by= (length(CO2)/n_sub))]
ti_temporal_sub <- rep(0,length(t_temporal_sub))
Fn_t_temporal_sub <- rank(t_temporal_sub)/length(t_temporal_sub)


#-------------------------------------------Objective functions---------------------------------#
#var_sub <- CO2_optimal_subsample
var_sub <- CO2_temporal_subsample
Fn_opt_sub <- seq(1/n_sub,1,length.out = n_sub)
Fn_limit <- seq(0,1,length.out = n_sub+1)
Time_limit <- quantile(Time,Fn_limit)
CO2_limit <- quantile(CO2,Fn_limit)
CO2_optimal_sub <- unname(quantile(CO2,Fn_opt_sub))


FuncO2 <- function(var_sub) {
  # var_sub <- CO2_optimal_subsample
  N_lags<- 10 #length(t)/2.5
  lag_value <- 16 # min(dist(t)) # or delta t
  pos_sub <- NULL
  for (i in 1:n_sub) {
    pos_sub <- rbind(pos_sub,which.min(abs(CO2-var_sub[i])))
  }
  t_sub <- t[pos_sub]
  ti_sub <- numeric(length(t_sub))
  Vario_var_sub <- variog(as.geodata(cbind(t_sub, ti_sub, var_sub), 
                                     coords.col = 1:2, data.col = 3), breaks = c(seq(0, lag_value * N_lags, lag_value)),
                          trend = "cte", lambda = 1, estimator.type = "classical", nugget.tolerance = 0, direction = 0, 
                          tolerance = 90, unit.angle = "degrees", pairs.min = 1)
  return (sum((Vario_var_sub$v-CO2_VarioEstimation$Semivarianzas)^2)) # # initial = 613.2208
  # Vario_var_sub <- Variograma(t_sub, ti_sub,
  #                             var_sub, 0, 90, N_lags, lag_value, 1, "", "") # CO2 temporal variogram
  #return (sum(abs(Vario_var_sub$Semivarianzas-CO2_VarioEstimation$Semivarianzas)))
}

#----------------Función conjuntas--------------------------------------------#
FuncConj  <- function(var_sub) {
  total =   FuncO2(var_sub)  
  return (total)
}

#-------------------------------optimization--------------------------------#
#n_sub = 12
Fn_optimal_t_sub <- seq(1/n_sub,1,length.out = n_sub)
Fn_opt_sub <- seq(1/n_sub,1,length.out = n_sub)
CO2_VarioEstimation$Semivarianzas
# stratified probability space
# CO2_optimal_sub1 <- res$sampled_data[,2]
# CO2_optimal_sub_order <- CO2_optimal_sub1[order(CO2_optimal_sub1)]
CO2_optimal_sub <- unname(quantile(CO2,Fn_opt_sub)) #CO2_optimal_sub_order    #
CO2_order <- CO2[order(CO2)]
pos_otp_order <- NULL
for (i in 1:n_sub) {
  pos_otp_order <- rbind(pos_otp_order,which.min(abs(CO2_order-CO2_optimal_sub[i])))
}
#pos_otp_order <- which(CO2_order %in% CO2_optimal_sub)
#group_sub <- t(matrix(CO2[order(CO2)], n_sub, byrow=TRUE)) #t(matrix(CO2, n_sub, byrow=TRUE)) #  #
#group_sub <- t(matrix(CO2, n_sub, byrow=TRUE))
tol_sub <- 200 #round(length(CO2)/n_sub)  # round(length(CO2)/(2*n_sub))
Pop_Inicial <- matrix(data=NA,nrow=tol_sub,ncol=n_sub)
for (i in 1:(n_sub)) {
  for (j in 1:tol_sub) {
    pos <- pos_otp_order[i]
    Pop_Inicial[j,i] <- CO2_order[pos-j]
  }
}
subvector_min = apply(Pop_Inicial, 2, min)
subvector_max = apply(Pop_Inicial, 2, max)
system.time(outDEoptim_CO2 <- DEoptim(fn=FuncConj, lower= subvector_min,
                                      upper= subvector_max,
                                      control = DEoptim.control(VTR = 0.0000001,strategy = 3, itermax =1000, reltol = 1e-8, NP= length(Pop_Inicial[,1]), CR = 0.5, F = 0.8, initialpop = Pop_Inicial)))

# user  system elapsed
# 44.813  13.145  57.946   # itermax =1000, n_sub = 48, tol = 5, outDEoptim_CO2$optim$bestval = 7.800542
var_sub <- as.numeric(outDEoptim_CO2$optim$bestmem)
Estadisticas(var_sub)

# var_sub <- as.numeric(outDEoptim_CO2$optim$bestmem)
# Estadisticas(var_sub)
minPoints<-4
maxPoints<-10
CO2_temp_opt_subsample_vario_model<- 3
CO2_temp_opt_subsample_nugget<- 4.0
CO2_temp_opt_subsample_sill_and_nugget<- 28
CO2_temp_opt_subsample_rank <- 80
pos_sub_CO2 <- NULL
for (i in 1:n_sub) {
  pos_sub_CO2 <- rbind(pos_sub_CO2,which.min(abs(CO2-var_sub[i])))
}
t_sub_CO2 <- t[pos_sub_CO2]
ti_sub_CO2 <- numeric(length(t_sub_CO2))
# pay attention
var_sub_CO2 <- CO2[pos_sub_CO2]
var_sub_Temp <- Temp[pos_sub_CO2]
CO2_temp_opt_subsample <- var_sub_CO2
CO2_temp_opt_subsample_Stat <- Estadisticas(CO2_temp_opt_subsample)
Fn_var_sub_CO2 <- rank(var_sub_CO2)/(length(var_sub_CO2))
Fn_var_sub <- rank(var_sub_CO2)/(length(var_sub_CO2))
cbind(var_sub_CO2, CO2_optimal_subsample)
#plot(var_sub_CO2, CO2_optimal_subsample)
cor(var_sub_Temp,var_sub_CO2)
plot(var_sub_Temp,var_sub_CO2)

#png(paste(ee_dir,"/CO2_temp_opt_subsample_VarioEyeEstimation.png",sep=""),   bg = "white", width = 1000, height = 1000, res = 200)
CO2_temp_opt_subsample_EyeModelVarioFit<-EyeModel(t_sub_CO2, ti_sub_CO2, 
                                                  CO2_temp_opt_subsample, 0, 90, N_lags, lag_value, 1, 
                                                  CO2_temp_opt_subsample_vario_model, CO2_temp_opt_subsample_nugget, CO2_temp_opt_subsample_sill_and_nugget, 
                                                  CO2_temp_opt_subsample_rank,
                                                  "") # A fitted temp_opt variogram model
#dev.off()
main_CO2 <-  expression( (a) ~  CO[2]: ~ 12 ~ samples)

#-------------------------------------------------------------------------------------#
#----------------------------------Comparing the results------------------------------#
#-------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
log_list=list(var_sub_CO2, var_sub_CH4, var_sub_N2O)
log_xlabel_list=list(CO2label, CH4label, N2Olabel)
n_logs= 3
plot_mean=TRUE
plot_median=TRUE
fontproportion=1.0
png(paste(ee_dir,"/CO2_CH4_N2O_Temp_Temporal_Distr_n_sub_12_opt.png",sep=""), bg = "white", width = 1300, height = 1500, res = 200)
par(mfrow = c(n_logs, 1), mar = c(4, 6, 4, 2), cex.lab= 1.3, cex.axis = 1.3)
p=CO2
s=summary(p)
plot(Time, p,  type = "p", xaxt = "n", main = "",
     ylab = CO2label , xlab = "", pch = 1, cex = 0.8)
axis(1, at= seq(0,850,50) )
#grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
#par(new=TRUE)
# plot(Time, p, type = "p", xaxt = "n",
#      ylab = CO2label , xlab = "", pch = 1)
# axis(1, at= seq(0,850,50) )
if (plot_mean) {
  abline(h = s[4], col = "Red", lwd = 2)
}  
if (plot_median) {
  abline(h = s[3], col = "Blue", lty = 2,  lwd = 2)
}
# legend(300, 150, legend=c("Mean", "Median"),
#        col=c("red", "blue"), lty=1:2, box.lty=0)
points(t_sub_CO2, var_sub_CO2, col = "red", lty=1, lwd = 1, pch = 19)
#points(Time[seq(from= 1,to =length(CO2), by= (length(CO2)/n_sub))], CO2_temporal_subsample, col = "blue", lty=1, lwd = 3, pch = 19)
abline(v=t_sub_CO2, col= 'red', lty =2, lwd=1)
# par(adj=0)  # Para justificar a izquierda
# text(x=10, y=23, '(a)', cex=2)
#-------
p=CH4
s=summary(p)
plot(Time, p,  type = "p", xaxt = "n", main = "",
     ylab = CH4label , xlab = "", pch = 1)
axis(1, at= seq(0,850,50) )
# grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
# par(new=TRUE)
# plot(Time, p, type = "p", xaxt = "n",
#      ylab = CH4label , xlab = "", pch = 1)
# axis(1, at= seq(0,850,50) )
if (plot_mean) {
  abline(h = s[4], col = "Red", lwd = 2)
}  
if (plot_median) {
  abline(h = s[3], col = "Blue", lty = 2, lwd = 2)
}
# legend(300, 150, legend=c("Mean", "Median"),
#        col=c("red", "blue"), lty=1:2, box.lty=0)
points(t_sub_CH4, var_sub_CH4, col = "red", lty=1, lwd = 1, pch = 19)
#points(Time[seq(from= 1,to =length(CH4), length.out =n_sub)], CH4_temporal_subsample, col = "blue", lty=1, lwd = 3, pch = 19)
abline(v=t_sub_CH4, col= 'red', lty =2, lwd=1)
# par(adj=0)  # Para justificar a izquierda
# text(x=10, y=0.0005, '(b)', cex=2)
# par(adj=0.5)
#------------
p=N2O
s=summary(p)
plot(Time, p,  type = "p", xaxt = "n", main = "",
     ylab = N2Olabel , xlab = Timelabel, pch = 1)
axis(1, at= seq(0,850,50) )
# grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
# par(new=TRUE)
# plot(Time, p, type = "p", xaxt = "n",
#      ylab = N2Olabel , xlab = Timelabel, pch = 1)
# axis(1, at= seq(0,850,50) )
if (plot_mean) {
  abline(h = s[4], col = "Red", lwd = 2)
}  
if (plot_median) {
  abline(h = s[3], col = "Blue", lty = 2, lwd = 2)
}
# legend(300, 150, legend=c("Mean", "Median"),
#        col=c("red", "blue"), lty=1:2, box.lty=0)
points(t_sub_N2O, var_sub_N2O, col = "red", lty=1, lwd = 1, pch = 19)
#points(Time[seq(from= 1,to =length(N2O), length.out =n_sub)], N2O_temporal_subsample, col = "blue", lty=1, lwd = 3, pch = 19)
abline(v=t_sub_N2O, col= 'red', lty =2, lwd=1)
# par(adj=0)  # Para justificar a izquierda
# text(x=10, y=0.032, '(c)', cex=2)
dev.off()
#-------------------------Univariate Statistics----------------------------#
cbind(CO2_Stat,CO2_temp_opt_subsample_Stat[,2], CO2_temporal_subsample_Stat[,2])

n_bins = histf_CO2$breaks
png(paste(ee_dir,"/CO2_HistBoxPlotFreq_temp_opt_12.png",sep=""), bg = "white",  width = 1000, height = 1000, res = 200)
histf_CO2_temp_opt_subsample <- HistBoxplot(x=CO2_temp_opt_subsample,
                                           mean = CO2_temp_opt_subsample_Stat[5,2], median = CO2_temp_opt_subsample_Stat[4,2], main ="", 
                                           xlab = CO2label, ylab = expression(bold("Relative frequency (%)")), AbsFreq = FALSE, PercentFreq = TRUE
                                           ,nbins = n_bins, cex.lab= 1.5, cex.axis = 1.3)
dev.off()

png(paste(ee_dir,"/CO2_HistBoxPlotFreq_temporal_12.png",sep=""), bg = "white",  width = 1000, height = 1000, res = 200)
histf_CO2_temporal_subsample <- HistBoxplot(x=CO2_temporal_subsample, 
                                            mean = CO2_temporal_subsample_Stat[5,2], median = CO2_temporal_subsample_Stat[4,2], main ="", 
                                            xlab = CO2label, ylab = expression(bold("Relative frequency (%)")), AbsFreq = FALSE, PercentFreq = TRUE
                                            ,nbins = n_bins, cex.lab= 1.5, cex.axis = 1.3)
dev.off()


n_bins = histf_CH4$breaks
#n_bins = nclass.Sturges(CH4_temp_opt_subsample)
png(paste(ee_dir,"/CH4_HistBoxPlotFreq_temp_opt_24.png",sep=""), bg = "white",  width = 1000, height = 1000, res = 200)
histf_CH4_temp_opt_subsample <- HistBoxplot(x=CH4_temp_opt_subsample,
                                            mean = CH4_temp_opt_subsample_Stat[5,2], median = CH4_temp_opt_subsample_Stat[4,2], main ="", 
                                            xlab = CH4label, ylab = expression(bold("Relative frequency (%)")), AbsFreq = FALSE, PercentFreq = TRUE
                                            ,nbins = n_bins, cex.lab= 1.5, cex.axis = 1.3)
dev.off()
#n_bins = nclass.Sturges(CH4_temporal_subsample)
png(paste(ee_dir,"/CH4_HistBoxPlotFreq_temporal_24.png",sep=""), bg = "white",  width = 1000, height = 1000, res = 200)
histf_CH4_temporal_subsample <- HistBoxplot(x=CH4_temporal_subsample, 
                                            mean = CH4_temporal_subsample_Stat[5,2], median = CH4_temporal_subsample_Stat[4,2], main ="", 
                                            xlab = CH4label, ylab = expression(bold("Relative frequency (%)")), AbsFreq = FALSE, PercentFreq = TRUE
                                            ,nbins = n_bins, cex.lab= 1.5, cex.axis = 1.3)
dev.off()

n_bins = histf_N2O$breaks 
png(paste(ee_dir,"/N2O_HistBoxPlotFreq_temp_opt_24.png",sep=""), bg = "white",  width = 1000, height = 1000, res = 200)
histf_N2O_temp_opt_subsample <- HistBoxplot(x=N2O_temp_opt_subsample,
                                            mean = N2O_temp_opt_subsample_Stat[5,2], median = N2O_temp_opt_subsample_Stat[4,2], main ="", 
                                            xlab = N2Olabel, ylab = expression(bold("Relative frequency (%)")), AbsFreq = FALSE, PercentFreq = TRUE
                                            ,nbins = n_bins, cex.lab= 1.5, cex.axis = 1.3)
dev.off()

#n_bins = nclass.Sturges(N2O_temporal_subsample)
png(paste(ee_dir,"/N2O_HistBoxPlotFreq_temporal_24.png",sep=""), bg = "white",  width = 1000, height = 1000, res = 200)
histf_N2O_temporal_subsample <- HistBoxplot(x=N2O_temporal_subsample, 
                                            mean = N2O_temporal_subsample_Stat[5,2], median = N2O_temporal_subsample_Stat[4,2], main ="", 
                                            xlab = N2Olabel, ylab = expression(bold("Relative frequency (%)")), AbsFreq = FALSE, PercentFreq = TRUE
                                            ,nbins = n_bins, cex.lab= 1.5, cex.axis = 1.3)
dev.off()




CO2_temp_opt_subsample_Stat <- Estadisticas(CO2_temp_opt_subsample)
Comparison_CO2_sps12_sts12_Stat <- cbind(CO2_Stat, CO2_temp_opt_subsample_Stat[,2], CO2_temporal_subsample_Stat[,2])
colnames(Comparison_CO2_sps12_sts12_Stat)[colnames(Comparison_CO2_sps12_sts12_Stat)=="Values"] <- "CO2_data"
colnames(Comparison_CO2_sps12_sts12_Stat)[colnames(Comparison_CO2_sps12_sts12_Stat)=="CO2_temp_opt_subsample_Stat[, 2]"] <- "CO2_sps_12"
colnames(Comparison_CO2_sps12_sts12_Stat)[colnames(Comparison_CO2_sps12_sts12_Stat)=="CO2_temporal_subsample_Stat[, 2]"] <- "CO2_sts_12"
write.csv(Comparison_CO2_sps12_sts12_Stat , file = "Results/Time_Series/AED/Comparison_CO2_sps12_sts12_Stat.csv")

CH4_temp_opt_subsample_Stat <- Estadisticas(CH4_temp_opt_subsample)
Comparison_CH4_sps12_sts12_Stat <- cbind(CH4_Stat, CH4_temp_opt_subsample_Stat[,2], CH4_temporal_subsample_Stat[,2])
colnames(Comparison_CH4_sps12_sts12_Stat)[colnames(Comparison_CH4_sps12_sts12_Stat)=="Values"] <- "CH4_data"
colnames(Comparison_CH4_sps12_sts12_Stat)[colnames(Comparison_CH4_sps12_sts12_Stat)=="CH4_temp_opt_subsample_Stat[, 2]"] <- "CH4_sps_12"
colnames(Comparison_CH4_sps12_sts12_Stat)[colnames(Comparison_CH4_sps12_sts12_Stat)=="CH4_temporal_subsample_Stat[, 2]"] <- "CH4_sts_12"
write.csv(Comparison_CH4_sps12_sts12_Stat , file = "Results/Time_Series/AED/Comparison_CH4_sps12_sts12_Stat.csv")

N2O_temp_opt_subsample_Stat <- Estadisticas(N2O_temp_opt_subsample)
Comparison_N2O_sps12_sts12_Stat <- cbind(N2O_Stat, N2O_temp_opt_subsample_Stat[,2], N2O_temporal_subsample_Stat[,2])
colnames(Comparison_N2O_sps12_sts12_Stat)[colnames(Comparison_N2O_sps12_sts12_Stat)=="Values"] <- "N2O_data"
colnames(Comparison_N2O_sps12_sts12_Stat)[colnames(Comparison_N2O_sps12_sts12_Stat)=="N2O_temp_opt_subsample_Stat[, 2]"] <- "N2O_sps_12"
colnames(Comparison_N2O_sps12_sts12_Stat)[colnames(Comparison_N2O_sps12_sts12_Stat)=="N2O_temporal_subsample_Stat[, 2]"] <- "N2O_sts_12"
write.csv(Comparison_N2O_sps12_sts12_Stat , file = "Results/Time_Series/AED/Comparison_N2O_sps12_sts12_Stat.csv")

t_sub_CO2_12 <- t_sub_CO2
t_sub_CH4_12 <- t_sub_CH4
t_sub_N2O_12 <- t_sub_N2O
temporal_sub_12 <- Time_temporal_subsample
Time_opt_temp_12 <- data.frame(cbind(t_sub_CO2_12,t_sub_CH4_12,t_sub_N2O_12, temporal_sub_12))
write.csv(Time_opt_temp_12 , file = "Results/Time_Series/AED/Time_opt_temp_12.csv")


# One-sample t-test
res_fixed_CO2 <- t.test(CO2_temporal_subsample, mu = mean(CO2), conf.level = 0.95)
res_opt_CO2 <- t.test(CO2_temp_opt_subsample, mu = mean(CO2), conf.level = 0.95)
# Printing the results
res_fixed_CO2
res_opt_CO2

ks_test_fixed_CO2 <- ks.test(CO2_temporal_subsample, CO2, conf.level = 0.95)
ks_test_opt_CO2 <- ks.test(CO2_temp_opt_subsample, CO2, conf.level = 0.95)
ks_test_fixed_CO2
ks_test_opt_CO2

# One-sample t-test
res_fixed_CH4 <- t.test(CH4_temporal_subsample, mu = mean(CH4), conf.level = 0.95)
res_opt_CH4 <- t.test(CH4_temp_opt_subsample, mu = mean(CH4), conf.level = 0.95)
# Printing the results
res_fixed_CH4
res_opt_CH4

ks_test_fixed_CH4 <- ks.test(CH4_temporal_subsample, CH4, conf.level = 0.95)
ks_test_opt_CH4 <- ks.test(CH4_temp_opt_subsample, CH4, conf.level = 0.95)
ks_test_fixed_CH4
ks_test_opt_CH4


# One-sample t-test
res_fixed_N2O <- t.test(N2O_temporal_subsample, mu = mean(N2O), conf.level = 0.95)
res_opt_N2O <- t.test(N2O_temp_opt_subsample, mu = mean(N2O), conf.level = 0.95)
# Printing the results
res_fixed_N2O
res_opt_N2O

ks_test_fixed_N2O <- ks.test(N2O_temporal_subsample, N2O, conf.level = 0.95)
ks_test_opt_N2O <- ks.test(N2O_temp_opt_subsample, N2O, conf.level = 0.95)
ks_test_fixed_N2O
ks_test_opt_N2O

#-----------------------------Variogram---------------------------------------#

png(paste(ee_dir,"/CO2_VarioEyeEstimation_overlap.png",sep=""),   bg = "white", width = 1000, height = 1000, res = 150)
CO2_EyeModelVarioFit<-EyeModel_overlap(t, ti, 
                                       CO2,  cbind(t_sub_CO2,t_temporal_sub), 
                                       cbind(ti_sub_CO2,ti_temporal_sub), cbind(CO2_temp_opt_subsample,CO2_temporal_subsample), 0, 90, N_lags, lag_value, 1, 
                                       CO2_vario_model, CO2_nugget, CO2_sill_and_nugget, CO2_rank,
                                       MainTitle="", xlab= Timelabel, ylab = expression(bold("Semivariance")) ) # A fitted temporal variogram model
dev.off()

png(paste(ee_dir,"/CH4_VarioEyeEstimation_overlap.png",sep=""),   bg = "white", width = 1000, height = 1000, res = 150)
CH4_EyeModelVarioFit<-EyeModel_overlap(t, ti, 
                                       CH4,  cbind(t_sub_CH4,t_temporal_sub), 
                                       cbind(ti_sub_CH4,ti_temporal_sub), cbind(CH4_temp_opt_subsample,CH4_temporal_subsample), 0, 90, N_lags, lag_value, 1, 
                                       CH4_vario_model, CH4_nugget, CH4_sill_and_nugget, CH4_rank,
                                       MainTitle="", xlab= Timelabel, ylab = expression(bold("Semivariance")) ) # A fitted temporal variogram model
dev.off()

png(paste(ee_dir,"/N2O_VarioEyeEstimation_overlap.png",sep=""),   bg = "white", width = 1000, height = 1000, res = 150)
N2O_EyeModelVarioFit<-EyeModel_overlap_N2O(t, ti, 
                                           N2O,  cbind(t_sub_N2O,t_temporal_sub), 
                                           cbind(ti_sub_N2O,ti_temporal_sub), cbind(N2O_temp_opt_subsample,N2O_temporal_subsample), 0, 90, N_lags, lag_value, 1, 
                                           N2O_vario_model, N2O_nugget, N2O_sill_and_nugget, N2O_rank,
                                           MainTitle="", xlab= Timelabel, ylab = expression(bold("Semivariance")) ) # A fitted temporal variogram model
dev.off()


png(paste(ee_dir,"/CO2_temp_opt_subsample_VarioEstimation_10lags.png",sep=""), bg = "white", width = 1000, height = 1000, res = 150)
CO2_temp_opt_subsample_VarioEstimation<-Variograma(t_sub_CO2, ti_sub_CO2, 
                                                   CO2_temp_opt_subsample, 0, 90, N_lags, lag_value, 1, MainTitle="Empirical variogram", xlab= "Time (days)", ylab = "Semivariogram") # CO2_temp_opt_subsample temporal variogram
dev.off()

png(paste(ee_dir,"/CO2_temporal_subsample_VarioEstimation_10lags.png",sep=""), bg = "white", width = 1000, height = 1000, res = 150)
CO2_temporal_subsample_VarioEstimation<-Variograma(t_temporal_sub, ti_temporal_sub, 
                                                   CO2_temporal_subsample, 0, 90, N_lags, lag_value, 1, MainTitle="Empirical variogram", xlab= "Time (days)", ylab = "Semivariogram") # CO2_temporal_subsample temporal variogram
dev.off()
# sum((CO2_temp_opt_subsample_VarioEstimation$Semivarianzas- CO2_VarioEstimation$Semivarianzas)^2)
sum(abs(CO2_temp_opt_subsample_VarioEstimation$Semivarianzas-CO2_VarioEstimation$Semivarianzas))
sum(abs(CO2_temporal_subsample_VarioEstimation$Semivarianzas-CO2_VarioEstimation$Semivarianzas))

png(paste(ee_dir,"/CH4_temp_opt_subsample_VarioEstimation_10lags.png",sep=""), bg = "white", width = 1000, height = 1000, res = 150)
CH4_temp_opt_subsample_VarioEstimation<-Variograma(t_sub_CH4, ti_sub_CH4, 
                                                   CH4_temp_opt_subsample, 0, 90, N_lags, lag_value, 1, MainTitle="Empirical variogram", xlab= "Time (days)", ylab = "Semivariogram") # CH4_temp_opt_subsample temporal variogram
dev.off()

png(paste(ee_dir,"/CH4_temporal_subsample_VarioEstimation_10lags.png",sep=""), bg = "white", width = 1000, height = 1000, res = 150)
CH4_temporal_subsample_VarioEstimation<-Variograma(t_temporal_sub, ti_temporal_sub, 
                                                   CH4_temporal_subsample, 0, 90, N_lags, lag_value, 1, MainTitle="Empirical variogram", xlab= "Time (days)", ylab = "Semivariogram") # CH4_temporal_subsample temporal variogram
dev.off()
# sum((CH4_temp_opt_subsample_VarioEstimation$Semivarianzas- CH4_VarioEstimation$Semivarianzas)^2)
sum(abs(CH4_temp_opt_subsample_VarioEstimation$Semivarianzas-CH4_VarioEstimation$Semivarianzas))
sum(abs(CH4_temporal_subsample_VarioEstimation$Semivarianzas-CH4_VarioEstimation$Semivarianzas))

png(paste(ee_dir,"/N2O_temp_opt_subsample_VarioEstimation_10lags.png",sep=""), bg = "white", width = 1000, height = 1000, res = 150)
N2O_temp_opt_subsample_VarioEstimation<-Variograma(t_sub_N2O, ti_sub_N2O, 
                                                   N2O_temp_opt_subsample, 0, 90, N_lags, lag_value, 1, MainTitle="Empirical variogram", xlab= "Time (days)", ylab = "Semivariogram") # N2O_temp_opt_subsample temporal variogram
dev.off()

png(paste(ee_dir,"/N2O_temporal_subsample_VarioEstimation_10lags.png",sep=""), bg = "white", width = 1000, height = 1000, res = 150)
N2O_temporal_subsample_VarioEstimation<-Variograma(t_temporal_sub, ti_temporal_sub, 
                                                   N2O_temporal_subsample, 0, 90, N_lags, lag_value, 1, MainTitle="Empirical variogram", xlab= "Time (days)", ylab = "Semivariogram") # N2O_temporal_subsample temporal variogram
dev.off()
# sum((N2O_temp_opt_subsample_VarioEstimation$Semivarianzas- N2O_VarioEstimation$Semivarianzas)^2)
sum(abs(N2O_temp_opt_subsample_VarioEstimation$Semivarianzas-N2O_VarioEstimation$Semivarianzas))
sum(abs(N2O_temporal_subsample_VarioEstimation$Semivarianzas-N2O_VarioEstimation$Semivarianzas))

#-----------------------------Calculation of cumulative GHG fluxes------------------------#
#---------------------------------Conversion----------------------------------------------#
# (μmol CO2 m^-2 s^-1) -----> (g C m^-2 day^-1)
#  μmol CO2 / 10^6     ----->  1 mol CO2 -----> 1 mol C -----> 12.01 g C 
#  μmol CO2 / 10^6     ----->  1 mol CO2 -----> 44.01 g CO2 
#  s^-1 *60*60*24      ----->  day^-1
CO2_GWP <- (CO2*60*60*24*44.01)/1000000
CO2_temp_opt_subsample_GWP <- (CO2_temp_opt_subsample*60*60*24*44.01)/1000000
CO2_temporal_subsample_GWP <- (CO2_temporal_subsample*60*60*24*44.01)/1000000
#
Mean_CO2_GWP <- mean(CO2_GWP)
sd_CO2_GWP <- sd(CO2_GWP)
#Cumsum_CO2_GWP <- cumsum(Mean_CO2_GWP*Time_GWP)
Mean_CO2_GWP*365
sd_CO2_GWP*365


Mean_CO2_temp_opt_subsample_GWP <- mean(CO2_temp_opt_subsample_GWP)
sd_CO2_temp_opt_subsample_GWP <- sd(CO2_temp_opt_subsample_GWP)
Mean_CO2_temp_opt_subsample_GWP*365
sd_CO2_temp_opt_subsample_GWP*365
Mean_CO2_temp_opt_subsample_GWP*length(CO2_GWP)/32
sd_CO2_temp_opt_subsample_GWP*length(CO2_GWP)/32
CI95_CO2_temp_opt_subsample_GWP <- unname(quantile(CO2_temp_opt_subsample_GWP, probs = c(0.025,0.975)))
CI95_CO2_temp_opt_subsample_GWP[1]*365
CI95_CO2_temp_opt_subsample_GWP[2]*365
CI95_CO2_temp_opt_subsample_GWP[1]*length(CO2_GWP)/32
CI95_CO2_temp_opt_subsample_GWP[2]*length(CO2_GWP)/32

Mean_CO2_temporal_subsample_GWP <- mean(CO2_temporal_subsample_GWP)
sd_CO2_temporal_subsample_GWP <- sd(CO2_temporal_subsample_GWP)
Mean_CO2_temporal_subsample_GWP*365
sd_CO2_temporal_subsample_GWP*365 
CI95_CO2_temporal_subsample_GWP <- unname(quantile(CO2_temporal_subsample_GWP, probs = c(0.025,0.975)))
CI95_CO2_temporal_subsample_GWP[1]*365
CI95_CO2_temporal_subsample_GWP[2]*365
CI95_CO2_temporal_subsample_GWP[1]*length(CO2_GWP)/32
CI95_CO2_temporal_subsample_GWP[2]*length(CO2_GWP)/32


Mean_CO2_temp_opt_subsample_GWP*258
sd_CO2_temp_opt_subsample_GWP*sqrt(258)
Mean_CO2_temporal_subsample_GWP*length(CO2_GWP)/32
sd_CO2_temporal_subsample_GWP*length(CO2_GWP)/32

CH4_GWP <- (CH4*60*60*24*16.04)/1000000
CH4_temp_opt_subsample_GWP <- (CH4_temp_opt_subsample*60*60*24*16.04)/1000000
CH4_temporal_subsample_GWP <- (CH4_temporal_subsample*60*60*24*16.04)/1000000
Mean_CH4_GWP <- mean(CH4_GWP)
sd_CH4_GWP <- sd(CH4_GWP)
#Cumsum_CH4_GWP <- cumsum(Mean_CH4_GWP*Time_GWP)
Mean_CH4_GWP*365
sd_CH4_GWP*365
Mean_CH4_GWP*length(CH4_GWP)/32
sd_CH4_GWP*length(CH4_GWP)/32

Mean_CH4_temp_opt_subsample_GWP <- mean(CH4_temp_opt_subsample_GWP)
sd_CH4_temp_opt_subsample_GWP <- sd(CH4_temp_opt_subsample_GWP)
Mean_CH4_temp_opt_subsample_GWP*365
sd_CH4_temp_opt_subsample_GWP*365
Mean_CH4_temp_opt_subsample_GWP*length(CH4_GWP)/32
sd_CH4_temp_opt_subsample_GWP*length(CH4_GWP)/32
CI95_CH4_temp_opt_subsample_GWP <- unname(quantile(CH4_temp_opt_subsample_GWP, probs = c(0.025,0.975)))
CI95_CH4_temp_opt_subsample_GWP[1]*365
CI95_CH4_temp_opt_subsample_GWP[2]*365
CI95_CH4_temp_opt_subsample_GWP[1]*length(CH4_GWP)/32
CI95_CH4_temp_opt_subsample_GWP[2]*length(CH4_GWP)/32

Mean_CH4_temporal_subsample_GWP <- mean(CH4_temporal_subsample_GWP)
sd_CH4_temporal_subsample_GWP <- sd(CH4_temporal_subsample_GWP)
Mean_CH4_temporal_subsample_GWP*365
sd_CH4_temporal_subsample_GWP*365
Mean_CH4_temporal_subsample_GWP*length(CH4_GWP)/32
sd_CH4_temporal_subsample_GWP*length(CH4_GWP)/32
CI95_CH4_temporal_subsample_GWP <- unname(quantile(CH4_temporal_subsample_GWP, probs = c(0.025,0.975)))
CI95_CH4_temporal_subsample_GWP[1]*365
CI95_CH4_temporal_subsample_GWP[2]*365
CI95_CH4_temporal_subsample_GWP[1]*length(CH4_GWP)/32
CI95_CH4_temporal_subsample_GWP[2]*length(CH4_GWP)/32

N2O_GWP <- (N2O*60*60*24*44.013)/1000000
N2O_temp_opt_subsample_GWP <- (N2O_temp_opt_subsample*60*60*24*44.013)/1000000
N2O_temporal_subsample_GWP <- (N2O_temporal_subsample*60*60*24*44.013)/1000000
Mean_N2O_GWP <- mean(N2O_GWP)
sd_N2O_GWP <- sd(N2O_GWP)
#Cumsum_N2O_GWP <- cumsum(Mean_N2O_GWP*Time_GWP)
Mean_N2O_GWP*365
sd_N2O_GWP*365
Mean_N2O_GWP*length(N2O_GWP)/32
sd_N2O_GWP*length(N2O_GWP)/32

Mean_N2O_temp_opt_subsample_GWP <- mean(N2O_temp_opt_subsample_GWP)
sd_N2O_temp_opt_subsample_GWP <- sd(N2O_temp_opt_subsample_GWP)
Mean_N2O_temp_opt_subsample_GWP*365
sd_N2O_temp_opt_subsample_GWP*365
Mean_N2O_temp_opt_subsample_GWP*length(N2O_GWP)/32
sd_N2O_temp_opt_subsample_GWP*length(N2O_GWP)/32
CI95_N2O_temp_opt_subsample_GWP <- unname(quantile(N2O_temp_opt_subsample_GWP, probs = c(0.025,0.975)))
CI95_N2O_temp_opt_subsample_GWP[1]*365
CI95_N2O_temp_opt_subsample_GWP[2]*365
CI95_N2O_temp_opt_subsample_GWP[1]*length(N2O_GWP)/32
CI95_N2O_temp_opt_subsample_GWP[2]*length(N2O_GWP)/32

Mean_N2O_temporal_subsample_GWP <- mean(N2O_temporal_subsample_GWP)
sd_N2O_temporal_subsample_GWP <- sd(N2O_temporal_subsample_GWP)
Mean_N2O_temporal_subsample_GWP*365
sd_N2O_temporal_subsample_GWP*365
Mean_N2O_temporal_subsample_GWP*length(N2O_GWP)/32
sd_N2O_temporal_subsample_GWP*length(N2O_GWP)/32
CI95_N2O_temporal_subsample_GWP <- unname(quantile(N2O_temporal_subsample_GWP, probs = c(0.025,0.975)))
CI95_N2O_temporal_subsample_GWP[1]*365
CI95_N2O_temporal_subsample_GWP[2]*365
CI95_N2O_temporal_subsample_GWP[1]*length(N2O_GWP)/32
CI95_N2O_temporal_subsample_GWP[2]*length(N2O_GWP)/32


