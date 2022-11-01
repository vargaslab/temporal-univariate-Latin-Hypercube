
#----------------------Optimization--------------------------------#
dir.create(paste(getwd(),"/Results/Time_Series/Optimization", sep=""))
ee_dir<-paste(root_dir,"/Results/Time_Series/Optimization",sep="")

#-----------------------------------------Sampling_12_N2O-----------------------#
n_sub <- 12
Fn_optimal_subsampling <- seq(1/n_sub,1,length.out = n_sub)
# subsampling by invevrse function
N2O_optimal_subsampling <- sapply(Fn_optimal_subsampling, Fn.inv.Bernshtein, N2O)
N2O_optimal_sub <- quantile(N2O,Fn_optimal_subsampling)
plot(N2O_optimal_subsampling,N2O_optimal_sub)
#Temp_optimal_subsampling <- sapply(Fn_optimal_subsampling, Fn.inv.Bernshtein, Temp)
optimal_time_position_N2O <- NULL
#optimal_time_position_Temp <- NULL
for (i in 1:n_sub) {
  optimal_time_position_N2O <- rbind(optimal_time_position_N2O, which.min(abs(N2O-N2O_optimal_subsampling[i])))
  #optimal_time_position_Temp <- rbind(optimal_time_position_Temp, which.min(abs(Temp-Temp_optimal_subsampling[i])))
}
(N2O_optimal_subsample <- N2O[optimal_time_position_N2O])
N2O_optimal_subsampling
N2O_optimal_subsample_Stat <- Estadisticas(N2O_optimal_subsample)
Fn_sub_opt <- rank(N2O_optimal_subsample)/length(N2O_optimal_subsample)
(sum(abs(Fn_sub_opt-Fn_opt_sub))) # 0.5

#------------------------------Random subsampling------------------------------#
# Random subsampling
Fn_random_subsampling <- runif(n_sub,0,1)
# subsampling by invevrse function
N2O_random_subsampling <- sapply(Fn_random_subsampling, Fn.inv.Bernshtein, N2O)
random_time_position <- NULL
for (i in 1:n_sub) {
  random_time_position <- rbind(random_time_position, which.min(abs(N2O-N2O_random_subsampling[i])))
}
(N2O_random_subsample <- N2O[random_time_position])
N2O_random_subsampling
N2O_random_subsample_Stat <- Estadisticas(N2O_random_subsample)
Fn_N2O_random_subsampling <- rank(N2O_random_subsampling)/length(N2O_random_subsampling)
(sum(abs(Fn_N2O_random_subsampling-Fn_opt_sub))) # 0.5

#--------------------temporal sub-sampling defined by user---------------------#
# Temporal subsampling
N2O_temporal_subsample <- N2O[seq(from= 1,to =length(N2O), by= (length(N2O)/n_sub))]
Fn_temporal_subsampling <- rank(N2O_temporal_subsample)/length(N2O_temporal_subsample)
N2O_temporal_subsample_Stat <- Estadisticas(N2O_temporal_subsample)
t_temporal_sub <- Time[seq(from= 1,to =length(N2O), by= (length(N2O)/n_sub))]
ti_temporal_sub <- rep(0,length(t_temporal_sub))
Fn_N2O_temporal_subsample <- rank(N2O_temporal_subsample)/length(N2O_temporal_subsample)
(sum(abs(Fn_N2O_temporal_subsample-Fn_optimal_subsampling))) # 0.5


#-------------------------------------------Objective functions---------------------------------#
# FuncO1  <- function(var_sub) {
#   Fn_var_sub <- rank(var_sub)/(length(var_sub))
#   return (sum((Fn_var_sub-Fn_opt_sub)^2)) # initial = 32.8, opt = 0.007018322
# }

var_sub <- N2O_optimal_subsample
FuncO2  <- function(var_sub) {
  # var_sub <- N2O_optimal_subsample
  N_lags<- 10 #length(t)/2.5
  lag_value <- 16 # min(dist(t)) # or delta t
  pos_sub <- NULL
  for (i in 1:n_sub) {
    pos_sub <- rbind(pos_sub,which.min(abs(N2O-var_sub[i])))
  }
  t_sub <- t[pos_sub]
  ti_sub <- numeric(length(t_sub))
  Vario_var_sub <- variog(as.geodata(cbind(t_sub, ti_sub, var_sub), 
                                     coords.col = 1:2, data.col = 3), breaks = c(seq(0, lag_value * N_lags, lag_value)),
                          trend = "cte", lambda = 1, estimator.type = "classical", nugget.tolerance = 0, direction = 0, 
                          tolerance = 90, unit.angle = "degrees", pairs.min = 1)
  return (sum((Vario_var_sub$v-N2O_VarioEstimation$Semivarianzas)^2)) # # initial = 1.140355e-08
  # Vario_var_sub <- Variograma(t_sub, ti_sub,
  #                             var_sub, 0, 90, N_lags, lag_value, 1, "", "") # N2O temporal variogram
  #return (sum(abs(Vario_var_sub$Semivarianzas-N2O_VarioEstimation$Semivarianzas)))
}

# FuncO3  <- function(var_sub) {
#   var_dife <- var(N2O)- var(var_sub)
#   return (abs(var_dife)) # initial = 1.518811
# }

#----------------Función conjuntas--------------------------------------------#
FuncConj  <- function(var_sub) {
  total =   FuncO2(var_sub)  
  return (total)
}

#-------------------------------optimization--------------------------------#
#n_sub = 48
Fn_opt_sub <- seq(1/n_sub,1,length.out = n_sub)
N2O_VarioEstimation$Semivarianzas
# stratified probability space
N2O_optimal_sub <- unname(quantile(N2O,Fn_opt_sub))
N2O_order <- N2O[order(N2O)]
pos_otp_order <- NULL
for (i in 1:n_sub) {
  pos_otp_order <- rbind(pos_otp_order,which.min(abs(N2O_order-N2O_optimal_sub[i])))
}
#pos_otp_order <- which(N2O_order %in% N2O_optimal_sub)
#group_sub <- t(matrix(N2O[order(N2O)], n_sub, byrow=TRUE)) #t(matrix(N2O, n_sub, byrow=TRUE)) #  # 
#group_sub <- t(matrix(N2O, n_sub, byrow=TRUE))
tol_sub <- round(length(N2O)/(n_sub)) 
Pop_Inicial <- matrix(data=NA,nrow=tol_sub,ncol=n_sub)
for (i in 1:(n_sub)) {
  for (j in 1:tol_sub) {
    pos <- pos_otp_order[i]
    Pop_Inicial[j,i] <- N2O_order[pos-j]
  }
}

subvector_min = apply(Pop_Inicial, 2, min)
subvector_max = apply(Pop_Inicial, 2, max)

system.time(outDEoptim_N2O <- DEoptim(fn=FuncConj, lower= subvector_min,
                                      upper= subvector_max,
                                      control = DEoptim.control(VTR = 1e-16,strategy = 2, itermax =200, reltol = 1e-32, NP= length(Pop_Inicial[,1]), CR = 0.5, F = 0.8, initialpop = Pop_Inicial)))

# user  system elapsed
var_sub <- as.numeric(outDEoptim_N2O$optim$bestmem)*1000
Estadisticas(var_sub)

minPoints<-4
maxPoints<-10
N2O_temp_opt_subsample_vario_model<- 2
N2O_temp_opt_subsample_nugget<- 2.6
N2O_temp_opt_subsample_sill_and_nugget<- 2.6
N2O_temp_opt_subsample_rank <- 0
pos_sub_N2O <- NULL
for (i in 1:n_sub) {
  pos_sub_N2O <- rbind(pos_sub_N2O,which.min(abs(N2O-var_sub[i])))
}
t_sub_N2O <- Time[pos_sub_N2O]
ti_sub_N2O <- numeric(length(t_sub_N2O))
# pay attention
var_sub_N2O <- N2O[pos_sub_N2O]
N2O_temp_opt_subsample <- var_sub_N2O
N2O_temp_opt_subsample_Stat <- Estadisticas(N2O_temp_opt_subsample)
Fn_var_sub_N2O <- rank(var_sub_N2O)/(length(var_sub_N2O))

main_N2O <- expression( (g) ~  N[2]*O: ~ 12 ~ samples)
