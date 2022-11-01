#### ------- RGeoestad script in 2D ------- ####
#--------------------------------------------------------------#
#                   Data Manipulation                          #
#--------------------------------------------------------------#

root_dir<-getwd()

# Reading Data File in ASCII format space separated (.txt)
# -999.25 for non available values
#Data_File <- read.table(file=file.choose(),header=TRUE,na.strings="-999.25")

# Reading Data File in ASCII format comma separated (.csv)
# -999.25 for non available values
Data_Frame <- read.csv(file=file.choose(),header=T,na.strings="-999.25")
Data_Frame <- na.omit(Data_Frame[,c(3,9,10,11)])
#Data_Frame <- na.omit(Data_Frame[,c(3,9,10,14,15,16)])
## Data_Frame_No_NA <- na.omit(Data_Frame)
#Data_Frame <- Data_Frame[seq(from=1, to = length(Data_Frame[,1]), by= 15),]

# Creates a folder to store results for a case
dir.create(paste(getwd(),"/Results/Time_Series", sep=""))
result_dir<-paste(root_dir,"/Results/Time_Series",sep="")

# Creates a folder to store results for AED
dir.create(paste(getwd(),"/Results/Time_Series/AED", sep=""))

aed_dir<-paste(result_dir,"/AED",sep="")

#Time<-Data_Frame$DOEc
Time<-Data_Frame$DOYc
#Temp<-Data_Frame$S_Temp
#VWC<-Data_Frame$VWC
CO2<-Data_Frame$CO2_flux
CH4<-Data_Frame$CH4_flux*1000  
N2O<-Data_Frame$N2O_flux*1000

  
#--------------------------------------------------------------#
#               Exploratory Data Analysis                      #
#--------------------------------------------------------------#

###------------ Univariate Data Analysis---------------------###

# Basic Statistics Estimation 

# Basic Statistics
Time_Stat<-Estadisticas(Time)
#Temp_Stat<-Estadisticas(Temp)
#VWC_Stat<-Estadisticas(VWC)
CO2_Stat<-Estadisticas(CO2)
CH4_Stat<-Estadisticas(CH4)
N2O_Stat<-Estadisticas(N2O)
(n_bins=ceiling(sqrt(CO2_Stat[1,2])))
n_bins = 10
Timelabel <- expression(bold("Time (days)"))
#Templabel <- expression(bold("Temperature (°C)"))
#VWClabel <- expression(bold(paste(Volumetric ~ Water ~ Content ~ (m^3 / m^3))))
CO2label <- expression(bold(paste(CO[2] ~ efflux ~ (μmol ~ m^-2 ~ s^-1))))
CH4label <- expression(bold(CH[4] ~ efflux ~ (nmol ~ m^-2 ~ s^-1)))
N2Olabel <- expression(bold(N[2]*O ~ efflux ~ (nmol ~ m^-2 ~ s^-1)))


#-----------------------------CO2 data----------------------------------------#
# Histogram and Boxplot
# Save plot as image in png format
png(paste(aed_dir,"/CO2_HistBoxPlotCounts.png",sep=""), bg = "white", width = 1000, height = 1000, res = 150)
histc_CO2 <- HistBoxplot(x=CO2, mean = CO2_Stat[5,2], median = CO2_Stat[4,2], main = "",  
            xlab = CO2label, ylab = "Absolute frequency (count)", AbsFreq = TRUE, PercentFreq = FALSE,
            nbin = n_bins)
dev.off()
png(paste(aed_dir,"/CO2_HistBoxPlotFreq.png",sep=""), bg = "white",  width = 1000, height = 1000, res = 150)
histf_CO2 <- HistBoxplot(x=CO2, mean = CO2_Stat[5,2], median = CO2_Stat[4,2], main ="", 
            xlab = CO2label, ylab = "Relative frequency (%)", AbsFreq = FALSE, PercentFreq = TRUE
            ,nbins = n_bins) # 
dev.off()


n_bins = 20
#-----------------------------CH4 data----------------------------------------#
# Histogram and Boxplot
# Save plot as image in png format
png(paste(aed_dir,"/CH4_HistBoxPlotCounts.png",sep=""), bg = "white", width = 1000, height = 1000, res = 150)
histc_CH4 <- HistBoxplot(x=CH4, mean = CH4_Stat[5,2], median = CH4_Stat[4,2], main = "",  
                         xlab = CH4label, ylab = "Absolute frequency (count)", AbsFreq = TRUE, PercentFreq = FALSE,
                         nbin = n_bins)
dev.off()
png(paste(aed_dir,"/CH4_HistBoxPlotFreq.png",sep=""), bg = "white",  width = 1000, height = 1000, res = 150)
histf_CH4 <- HistBoxplot(x=CH4, mean = CH4_Stat[5,2], median = CH4_Stat[4,2], main ="", 
                         xlab = CH4label, ylab = "Relative frequency (%)", AbsFreq = FALSE, PercentFreq = TRUE,
                         nbin =n_bins)
dev.off()


n_bins = 30
#-----------------------------N2O data----------------------------------------#
# Histogram and Boxplot
# Save plot as image in png format
png(paste(aed_dir,"/N2O_HistBoxPlotCounts.png",sep=""), bg = "white", width = 1000, height = 1000, res = 150)
histc_N2O <- HistBoxplot(x=N2O, mean = N2O_Stat[5,2], median = N2O_Stat[4,2], main = "",  
                         xlab = N2Olabel, ylab = "Absolute frequency (count)", AbsFreq = TRUE, PercentFreq = FALSE,
                         nbin = n_bins)
dev.off()
png(paste(aed_dir,"/N2O_HistBoxPlotFreq.png",sep=""), bg = "white",  width = 1000, height = 1000, res = 150)
histf_N2O <- HistBoxplot(x=N2O, mean = N2O_Stat[5,2], median = N2O_Stat[4,2], main ="", 
                         xlab = N2Olabel, ylab = "Relative frequency (%)", AbsFreq = FALSE, PercentFreq = TRUE
                         , nbin =n_bins)
dev.off()


###------------ Bivariate Data Analysis---------------------###


png(paste(aed_dir,"/CO2-CH4_ScatterPlot.png",sep=""), bg = "white",  width = 1500, height = 1500, res = 230)
ScatterPlot(CH4, CO2, 15, 
            X = CH4_Stat[2,2], Xmax = CH4_Stat[7,2], 
            Ymin = CO2_Stat[2,2], Ymax = CO2_Stat[7,2], 
            XLAB = CH4label, YLAB = CO2label)
dev.off()

png(paste(aed_dir,"/CO2-N2O_ScatterPlot.png",sep=""), bg = "white",  width = 1500, height = 1500, res = 230)
ScatterPlot(N2O, CO2, 15, 
            X = N2O_Stat[2,2], Xmax = N2O_Stat[7,2], 
            Ymin = CO2_Stat[2,2], Ymax = CO2_Stat[7,2], 
            XLAB = N2Olabel, YLAB = CO2label)
dev.off()

png(paste(aed_dir,"/CH4-N2O_ScatterPlot.png",sep=""), bg = "white",  width = 1500, height = 1500, res = 230)
ScatterPlot(N2O, CH4, 15, 
            X = N2O_Stat[2,2], Xmax = N2O_Stat[7,2], 
            Ymin = CH4_Stat[2,2], Ymax = CH4_Stat[7,2], 
            XLAB = N2Olabel, YLAB = CH4label)
dev.off()

cor1 <- cor.test(CH4, CO2, conf.level = 0.95, method = "pearson")
cor2 <- cor.test(N2O, CO2, conf.level = 0.95, method = "pearson")
cor3 <- cor.test(N2O, CH4, conf.level = 0.95, method = "pearson")
cor1s <- cor.test(CH4, CO2, conf.level = 0.95, method = "spearman")
cor2s <- cor.test(N2O, CO2, conf.level = 0.95, method = "spearman")
cor3s <- cor.test(N2O, CH4, conf.level = 0.95, method = "spearman")
cor1k <- cor.test(CH4, CO2, conf.level = 0.95, method = "kendall")
cor2k <- cor.test(N2O, CO2, conf.level = 0.95, method = "kendall")
cor3k <- cor.test(N2O, CH4, conf.level = 0.95, method = "kendall")
cor1; cor2; cor3; cor1s; cor2s; cor3s; cor1k; cor2k; cor3k

