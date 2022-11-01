
#----------------------Optimization--------------------------------#
dir.create(paste(getwd(),"/Results/Time_Series/Optimization", sep=""))
ee_dir<-paste(root_dir,"/Results/Time_Series/Optimization",sep="")

#-----------------------------------------Sampling_12_CH4-----------------------#
n_sub <- 12
Fn_optimal_subsampling <- seq(1/n_sub,1,length.out = n_sub)
# subsampling by invevrse function
CH4_optimal_subsampling <- sapply(Fn_optimal_subsampling, Fn.inv.Bernshtein, CH4)
CH4_optimal_sub <- quantile(CH4,Fn_optimal_subsampling)
plot(CH4_optimal_subsampling,CH4_optimal_sub)
#Temp_optimal_subsampling <- sapply(Fn_optimal_subsampling, Fn.inv.Bernshtein, Temp)
optimal_time_position_CH4 <- NULL
#optimal_time_position_Temp <- NULL
for (i in 1:n_sub) {
  optimal_time_position_CH4 <- rbind(optimal_time_position_CH4, which.min(abs(CH4-CH4_optimal_subsampling[i])))
  #optimal_time_position_Temp <- rbind(optimal_time_position_Temp, which.min(abs(Temp-Temp_optimal_subsampling[i])))
}
(CH4_optimal_subsample <- CH4[optimal_time_position_CH4])
CH4_optimal_subsampling
CH4_optimal_subsample_Stat <- Estadisticas(CH4_optimal_subsample)
Fn_sub_opt <- rank(CH4_optimal_subsample)/length(CH4_optimal_subsample)
(sum(abs(Fn_sub_opt-Fn_opt_sub))) # 0.5

#------------------------------Random subsampling------------------------------#
# Random subsampling
Fn_random_subsampling <- runif(n_sub,0,1)
# subsampling by invevrse function
CH4_random_subsampling <- sapply(Fn_random_subsampling, Fn.inv.Bernshtein, CH4)
random_time_position <- NULL
for (i in 1:n_sub) {
  random_time_position <- rbind(random_time_position, which.min(abs(CH4-CH4_random_subsampling[i])))
}
(CH4_random_subsample <- CH4[random_time_position])
CH4_random_subsampling
CH4_random_subsample_Stat <- Estadisticas(CH4_random_subsample)
Fn_CH4_random_subsampling <- rank(CH4_random_subsampling)/length(CH4_random_subsampling)
(sum(abs(Fn_CH4_random_subsampling-Fn_opt_sub))) # 0.5

#--------------------temporal sub-sampling defined by user---------------------#
# Temporal subsampling
CH4_temporal_subsample <- CH4[seq(from= 1,to =length(CH4), by= (length(CH4)/n_sub))]
Fn_temporal_subsampling <- rank(CH4_temporal_subsample)/length(CH4_temporal_subsample)
CH4_temporal_subsample_Stat <- Estadisticas(CH4_temporal_subsample)
t_temporal_sub <- Time[seq(from= 1,to =length(CH4), by= (length(CH4)/n_sub))]
ti_temporal_sub <- rep(0,length(t_temporal_sub))
Fn_CH4_temporal_subsample <- rank(CH4_temporal_subsample)/length(CH4_temporal_subsample)
(sum(abs(Fn_CH4_temporal_subsample-Fn_optimal_subsampling))) # 0.5

#-------------------------------------------Objective functions---------------------------------#
# FuncO1  <- function(var_sub) {
#   Fn_var_sub <- rank(var_sub)/(length(var_sub))
#   return (sum((Fn_var_sub-Fn_opt_sub)^2)) # initial = 32.8, opt = 0.007018322
# }

var_sub <- CH4_optimal_subsample
FuncO2  <- function(var_sub) {
  # var_sub <- CH4_optimal_subsample
  N_lags<- 10 #length(t)/2.5
  lag_value <- 16 # min(dist(t)) # or delta t
  pos_sub <- NULL
  for (i in 1:n_sub) {
    pos_sub <- rbind(pos_sub,which.min(abs(CH4-var_sub[i])))
  }
  t_sub <- t[pos_sub]
  ti_sub <- numeric(length(t_sub))
  Vario_var_sub <- variog(as.geodata(cbind(t_sub, ti_sub, var_sub), 
                                     coords.col = 1:2, data.col = 3), breaks = c(seq(0, lag_value * N_lags, lag_value)),
                          trend = "cte", lambda = 1, estimator.type = "classical", nugget.tolerance = 0, direction = 0, 
                          tolerance = 90, unit.angle = "degrees", pairs.min = 1)
  return (sum((Vario_var_sub$v-CH4_VarioEstimation$Semivarianzas)^2)) # # initial = 9.635893e-14
  # Vario_var_sub <- Variograma(t_sub, ti_sub,
  #                             var_sub, 0, 90, N_lags, lag_value, 1, "", "") # CH4 temporal variogram
  #return (sum(abs(Vario_var_sub$Semivarianzas-CH4_VarioEstimation$Semivarianzas)))
}

# FuncO3  <- function(var_sub) {
#   var_dife <- var(CH4)- var(var_sub)
#   return (abs(var_dife)) # initial = 1.518811
# }

#----------------Función conjuntas--------------------------------------------#
FuncConj  <- function(var_sub) {
  total = FuncO2(var_sub)  
  return (total)
}

#-------------------------------optimization--------------------------------#
#n_sub = 12
Fn_optimal_t_sub <- seq(1/n_sub,1,length.out = n_sub)
Fn_opt_sub <- seq(1/n_sub,1,length.out = n_sub)
CH4_VarioEstimation$Semivarianzas
# stratified probability space
CH4_optimal_sub <- unname(quantile(CH4,Fn_opt_sub))
CH4_order <- CH4[order(CH4)]
pos_otp_order <- NULL
for (i in 1:n_sub) {
  pos_otp_order <- rbind(pos_otp_order,which.min(abs(CH4_order-CH4_optimal_sub[i])))
}
#pos_otp_order <- which(CH4_order %in% CH4_optimal_sub)
#group_sub <- t(matrix(CH4[order(CH4)], n_sub, byrow=TRUE)) #t(matrix(CH4, n_sub, byrow=TRUE)) #  #
#group_sub <- t(matrix(CH4, n_sub, byrow=TRUE))
tol_sub <- 200 #round(length(CH4)/n_sub)  # round(length(CH4)/(2*n_sub))
Pop_Inicial <- matrix(data=NA,nrow=tol_sub,ncol=n_sub)
for (i in 1:(n_sub)) {
  for (j in 1:tol_sub) {
    pos <- pos_otp_order[i]
    Pop_Inicial[j,i] <- CH4_order[pos-j]
  }
}

subvector_min = apply(Pop_Inicial, 2, min)
subvector_max = apply(Pop_Inicial, 2, max)

system.time(outDEoptim_CH4 <- DEoptim(fn=FuncO2, lower= subvector_min,
                                      upper= subvector_max,
                                      control = DEoptim.control(VTR = 1e-20, strategy = 2, itermax =1000, reltol = 1e-32, NP= length(Pop_Inicial[,1]), CR = 0.5, F = 0.8, initialpop = Pop_Inicial)))

# user  system elapsed CH4
# 162.748   8.323 171.348  # itermax =1000, n_sub = 48, tol = 20,  outDEoptim_CH4$optim$bestval = 2.649979e-17
var_sub <- as.numeric(outDEoptim_CH4$optim$bestmem)*1000
Estadisticas(var_sub)

minPoints<-4
maxPoints<-10
CH4_temp_opt_subsample_vario_model<- 2
CH4_temp_opt_subsample_nugget<- 0.07
CH4_temp_opt_subsample_sill_and_nugget<- 0.15
CH4_temp_opt_subsample_rank <- 110
pos_sub <- NULL
for (i in 1:n_sub) {
  pos_sub <- rbind(pos_sub,which.min(abs(CH4-var_sub[i])))
}
t_sub_CH4 <- Time[pos_sub]
ti_sub_CH4 <- numeric(length(t_sub_CH4))
# pay attention
var_sub_CH4 <- CH4[pos_sub]
CH4_temp_opt_subsample <- var_sub_CH4
CH4_temp_opt_subsample_Stat <- Estadisticas(CH4_temp_opt_subsample)
Fn_var_sub <- rank(var_sub)/(length(var_sub))

main_CH4 <- expression( (d) ~  CH[4]: ~ 12 ~ samples)
