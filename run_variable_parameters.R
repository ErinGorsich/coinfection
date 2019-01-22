#####################################################
#####################################################
#####################################################
# Packages and source code
#####################################################
#####################################################
#####################################################
rm(list = ls())
library("deSolve")
library("doParallel")
library("foreach")

setwd("~/git/coinfection")
data.dir <- "~/Documents/collaborations/Bree-Carrie-BTBcoinfection_paper/code/24-May-2017 update/R figures/"
# Parameter files (must read fixed first)
source('fixedparams.R', chdir = TRUE)
source('chronicmildparams.R', chdir = TRUE)

# RHS & Summary stats calculations
source('rhs.R', chdir = TRUE)
source('get_summary_stats.R', chdir = TRUE)

bins <- 50  # currently at 12

#####################################################
#####################################################
#####################################################
# run model analyses
#####################################################
#####################################################
#####################################################

#####################################################
# Chronic, density dependence model
#####################################################
params <- c(fixedparams, chronicmildparams)

# make matrix for levelplots
dfDD <- data.frame(
    beta_tc = rep(seq(params$beta_tcL, params$beta_tcU, 
        length.out = bins), times = bins),
    alpha_tc = rep(seq(params$alpha_tcL, params$alpha_tcU, 
        length.out = bins), each = bins),
    type = rep("variable", bins*bins),
    RoCT = NA, EE_C = NA, EE_TB = NA, 
    EE_CinS = NA, EE_CinTB = NA, FinalN = NA)

# make matrix for lines manipulating one parameter at a time
dfDD2 <- data.frame(
	beta_tc = c(rep(params$beta_tc, bins), 
		seq(params$beta_tcL, params$beta_tcU, 
		length.out = bins)), 
	alpha_tc = c(seq(params$alpha_tcL, params$alpha_tcU, 
		length.out = bins),
		rep(params$alpha_tc, bins)),
	type = rep("observed", bins*2),
	RoCT = NA, EE_C = NA, EE_TB = NA, 
	EE_CinS = NA, EE_CinTB = NA, FinalN = NA)	
dfDD <- rbind(dfDD, dfDD2)
dfDD$rowid <- seq(1, length(dfDD[,1]), 1)
dfDD <- dfDD[, c(1:3, 10)]

# start cluster
cl <- makeCluster(3)
registerDoParallel(cl)

DDchronic <- foreach (d = iter(dfDD, by = "row"), .combine = rbind,
    .packages = "deSolve") %dopar%{
    params <- c(fixedparams, chronicmildparams)
    params$beta_tc <- d$beta_tc
    params$alpha_tc <- d$alpha_tc
    rowid <- d$rowid
    
    val <- get_chronic_summary_stats_DD(params)
    
    df <- data.frame(rowid = rowid,
        RoCT = val[1], EE_C = val[2], EE_TB = val[3],
        EE_CinS = val[4], EE_CinTB = val[5], FinalN = val[6])
    }
stopCluster(cl)
dfDDchronic <- as.data.frame(cbind(dfDD, DDchronic))

# save results
dfDDchronic$change_beta_tc <- dfDDchronic$beta_tc/params$beta_c
dfDDchronic$change_alpha_tc <- dfDDchronic$alpha_tc/ params$alpha_c
write.csv(dfDDchronic, 
    paste(data.dir, "chronic_densitydependent_2019.csv", sep = ""))
rm(params)

#####################################################
# Chronic frequency dependent model
#####################################################
params <- c(fixedparams, chronicmildparamsFD)

dfFD <- data.frame(
	beta_tc = rep(seq(params$beta_tcL, params$beta_tcU, 
		length.out = bins), times = bins),
	alpha_tc = rep(seq(params$alpha_tcL, params$alpha_tcU, 
	   length.out = bins), each = bins),
	type = rep("variable", bins*bins),
	RoCT = NA, EE_C = NA, EE_TB = NA, 
	EE_CinS = NA, EE_CinTB = NA, FinalN = NA)

dfFD2 <- data.frame(
	beta_tc = c(rep(params$beta_tc, bins), 
		seq(params$beta_tcL, params$beta_tcU, 
		length.out = bins)), 
	alpha_tc = c(seq(params$alpha_tcL, params$alpha_tcU, 
	   length.out = bins), 
	rep(params$alpha_tc, bins)),
	type = rep("observed", bins*2),
	RoCT = NA, EE_C = NA, EE_TB = NA, 
	EE_CinS = NA, EE_CinTB = NA, FinalN = NA )	
dfFD <- rbind(dfFD, dfFD2)	
dfFD$rowid <- seq(1, length(dfFD[,1]), 1)
dfFD <- dfFD[, c(1:3, length(dfFD))]
	
# start cluster
cl <- makeCluster(3)
registerDoParallel(cl)

FDchronic <- foreach (d = iter(dfFD, by = "row"), .combine = rbind,
    .packages = "deSolve") %dopar%{
    params <- c(fixedparams, chronicmildparamsFD)
    params$beta_tc <- d$beta_tc
    params$alpha_tc <- d$alpha_tc
    rowid <- d$rowid
    
    val <- get_chronic_summary_stats_FD(params)

    df <- data.frame(rowid = rowid,
        RoCT = val[1], EE_C = val[2], EE_TB = val[3],
        EE_CinS = val[4], EE_CinTB = val[5], FinalN = val[6])
}
stopCluster(cl)
dfFDchronic <- as.data.frame(cbind(dfFD, FDchronic))

dfFDchronic$change_beta_tc <- dfFDchronic$beta_tc/params$beta_c
dfFDchronic$change_alpha_tc <- dfFDchronic$alpha_tc/ params$alpha_c
write.csv(dfFDchronic, 
    paste(data.dir, "chronic_frequencydependent_2019.csv", sep = ""))
rm(params)
print("Chronic models finished")

#####################################################
# Acute Density dependent model, gamma
#####################################################
source('~/Documents/collaborations/Bree-Carrie-BTBcoinfection_paper/code/24-May-2017 update/R figures/acutemildparams.R', chdir = TRUE)

params <- c(fixedparams, acutemildhightrans)

dfDD <- data.frame(
    beta_tu = rep(seq(params$beta_tuL, params$beta_tuU, 
        length.out = bins), times = bins),
    gamma_tu = rep(seq(params$gamma_tuL, params$gamma_tuU, 
        length.out = bins), each = bins),
    type = rep("variable", bins*bins),
    RoAT = NA, MaxA = NA, TimeMax = NA, MaxR = NA,
    EE_FinalN = NA, EE_R = NA, EE_TB = NA, EE_I = NA,
    EE_RinnoTB = NA, EE_RinTB = NA)

dfDD2 <- data.frame(
    beta_tu = c(rep(params$beta_tu, bins), 
        seq(params$beta_tuL, params$beta_tuU, 
            length.out = bins)), 
    gamma_tu = c(seq(params$gamma_tuL, params$gamma_tuU, 
        length.out = bins), rep(params$gamma_tu, bins)),
    type = rep("observed", bins*2),
    RoAT = NA, MaxA = NA, TimeMax = NA, MaxR = NA,
    EE_FinalN = NA, EE_R = NA, EE_TB = NA, EE_I = NA,
    EE_RinnoTB = NA, EE_RinTB = NA)	
dfDD <- rbind(dfDD, dfDD)	
dfDD$rowid <- seq(1, length(dfDD[,1]), 1)
dfDD <- dfDD[, c(1:3, length(dfDD))]

# start cluster
cl <- makeCluster(3)
registerDoParallel(cl)

DDacute <- foreach (d = iter(dfDD, by = "row"), .combine = rbind,
    .packages = "deSolve") %dopar%{
    params <- c(fixedparams, chronicmildparams)
    params$beta_tu <- d$beta_tu
    params$gamma_tu <- d$gamma_tu
    rowid <- d$rowid
                          
    val <- get_acute_summary_stats_DD(params)
                          
    df <- data.frame(rowid = rowid,
        RoAT = val[1], MaxA = val[2], TimeMax = val[3], MaxR = val[4],
        EE_FinalN = val[5], EE_R = val[6], EE_TB = val[7], EE_I = val[8],
        EE_RinnoTB = val[9], EE_RinTB = val[10])
}
stopCluster(cl)
dfDDacute <- as.data.frame(cbind(dfDD, DDacute))

dfDDacute$change_beta_tu <- dfDDacute$beta_tu/params$beta_u
dfDDacute$change_gamma_tu <- dfDDacute$gamma_tu/ params$gamma_u
dfDDacute_recov<- dfDDacute
write.csv(dfDDacute_recov, 
    paste(data.dir, "acute_densitydependent_recov_2019.csv", sep = ""))
rm(params)

#####################################################
# Acute Density dependent model, mortality
#####################################################
params <- c(fixedparams, acutemildhightrans)

dfDD <- data.frame(
    beta_tu = rep(seq(params$beta_tuL, params$beta_tuU, 
        length.out = bins), times = bins),
    alpha_tu = rep(seq(params$alpha_u, 4*params$alpha_u, 
        length.out = bins), each = bins),
    type = rep("variable", bins*bins),
    RoAT = NA, MaxA = NA, TimeMax = NA, MaxR = NA,
    EE_FinalN = NA, EE_R = NA, EE_TB = NA, EE_I = NA,
    EE_RinnoTB = NA, EE_RinTB = NA)

dfDD2 <- data.frame(
    beta_tu = c(rep(params$beta_tu, bins), 
                seq(params$beta_tuL, params$beta_tuU, 
                    length.out = bins)), 
    alpha_tu = c(seq(params$alpha_tuL, params$alpha_tuU, 
                     length.out = bins), rep(params$alpha_tu, bins)),
    type = rep("observed", bins*2),
    RoAT = NA, MaxA = NA, TimeMax = NA, MaxR = NA,
    EE_FinalN = NA, EE_R = NA, EE_TB = NA, EE_I = NA,
    EE_RinnoTB = NA, EE_RinTB = NA)	
dfDD <- rbind(dfDD, dfDD)	
dfDD$rowid <- seq(1, length(dfDD[,1]), 1)
dfDD <- dfDD[, c(1:3, length(dfDD))]

# start cluster
cl <- makeCluster(3)
registerDoParallel(cl)

DDacute <- foreach (d = iter(dfDD, by = "row"), .combine = rbind,
    .packages = "deSolve") %dopar%{
    params <- c(fixedparams, chronicmildparams)
    params$beta_tu <- d$beta_tu
    params$alpha_tu <- d$alpha_tu
    rowid <- d$rowid
                        
    val <- get_acute_summary_stats_DD(params)
                        
    df <- data.frame(rowid = rowid,
        RoAT = val[1], MaxA = val[2], TimeMax = val[3], MaxR = val[4],
        EE_FinalN = val[5], EE_R = val[6], EE_TB = val[7], EE_I = val[8],
        EE_RinnoTB = val[9], EE_RinTB = val[10])
                    }
stopCluster(cl)
dfDDacute <- as.data.frame(cbind(dfDD, DDacute))

# DeSolve errors
dfDDacute$change_beta_tu <- dfDDacute$beta_tu/params$beta_u
dfDDacute$change_alpha_tu <- dfDDacute$alpha_tu/ params$alpha_u
write.csv(dfDDacute, 
    paste(data.dir, "acute_densitydependent_2019.csv", sep = ""))
rm(params)

#####################################################
# Acute Frequency dependent model 
#####################################################
params <- c(fixedparams, acutemildhightransFD)

dfFD <- data.frame(
	beta_tu = rep(seq(params$beta_tuL, params$beta_tuU, 
		length.out = bins), times = bins),
	gamma_tu = rep(seq(params$gamma_tuL, params$gamma_tuU, 
        length.out = bins), each = bins),
	type = rep("variable", bins*bins),
	RoAT = NA, MaxA = NA, TimeMax = NA, MaxR = NA,
	EE_FinalN = NA, 
	EE_R = NA, EE_TB = NA, EE_I = NA,
	EE_RinnoTB = NA, EE_RinTB = NA)

dfFD2 <- data.frame(
	beta_tu = c(rep(params$beta_tu, bins), 
		seq(params$beta_tuL, params$beta_tuU, 
		length.out = bins)), 
	gamma_tu = c(seq(params$gamma_tuL, params$gamma_tuU, 
	    length.out = bins), 
		rep(params$gamma_tu, bins)),
    	type = rep("observed", bins*2),
	RoAT = NA, MaxA = NA, TimeMax = NA, MaxR = NA,
    EE_FinalN = NA, EE_R = NA, EE_TB = NA, EE_I = NA,
    EE_RinnoTB = NA, EE_RinTB = NA)
dfFD <- rbind(dfFD, dfFD2)	
dfFD$rowid <- seq(1, length(dfFD[,1]), 1)
dfFD <- dfFD[, c(1:3, length(dfFD))]

# start cluster
cl <- makeCluster(3)
registerDoParallel(cl)

FDacute <- foreach (d = iter(dfDD, by = "row"), .combine = rbind,
    .packages = "deSolve") %dopar%{
    params <- c(fixedparams, chronicmildparams)
    params$beta_tu <- d$beta_tu
    params$gamma_tu <- d$gamma_tu
    rowid <- d$rowid

    val <- get_acute_summary_stats_FD(params)
    df <- data.frame(rowid = rowid,
        RoAT = val[1], MaxA = val[2], TimeMax = val[3], MaxR = val[4],
        EE_FinalN = val[5], EE_R = val[6], EE_TB = val[7], EE_I = val[8],
        EE_RinnoTB = val[9], EE_RinTB = val[10])
}
stopCluster(cl)
dfFDacute <- as.data.frame(cbind(dfFD, FDacute))

dfFDacute$change_beta_tu <- dfFDacute$beta_tu/ params$beta_u
dfFDacute$change_gamma_tu <- dfFDacute$gamma_tu / params$gamma_u
dfFDacute_recov<- dfFDacute
write.csv(dfFDacute_recov, 
    paste(data.dir, "acute_frequencydependent_recov_2019.csv", sep = ""))
rm(params)

#####################################################
# Acute Frequency dependent model, mortality 
#####################################################
params <- c(fixedparams, acutemildhightransFD)

dfFD <- data.frame(
	beta_tu = rep(seq(params$beta_tuL, params$beta_tuU, 
	   length.out = bins), times = bins),
	alpha_tu = rep(seq(params$alpha_u, 4*params$alpha_u, 
	   length.out = bins), each = bins),
	type = rep("variable", bins*bins),
	RoAT = NA, MaxA = NA, TimeMax = NA, MaxR = NA,
	EE_FinalN = NA, EE_R = NA, EE_TB = NA, EE_I = NA,
	EE_RinnoTB = NA, EE_RinTB = NA)

dfFD2 <- data.frame(
	beta_tu = c(rep(params$beta_tu, bins), 
		seq(params$beta_tuL, params$beta_tuU, 
		length.out = bins)), 
	alpha_tu = c(seq(params$alpha_tuL, params$alpha_tuU, 
	   length.out = bins), rep(params$alpha_tu, bins)),
	type = rep("observed", bins*2),
	RoAT = NA, MaxA = NA, TimeMax = NA, MaxR = NA,
    EE_FinalN = NA, EE_R = NA, EE_TB = NA, EE_I = NA,
    EE_RinnoTB = NA, EE_RinTB = NA)
dfFD <- rbind(dfFD, dfFD2)	
dfFD$rowid <- seq(1, length(dfFD[,1]), 1)
dfFD <- dfFD[, c(1:3, length(dfFD))]

# start cluster
cl <- makeCluster(3)
registerDoParallel(cl)

FDacute <- foreach (d = iter(dfDD, by = "row"), .combine = rbind,
    .packages = "deSolve") %dopar%{
    params <- c(fixedparams, chronicmildparams)
    params$beta_tu <- d$beta_tu
    params$alpha_tu <- d$alpha_tu
    rowid <- d$rowid

    val <- get_acute_summary_stats_FD(params)
    df <- data.frame(rowid = rowid,
        RoAT = val[1], MaxA = val[2], TimeMax = val[3], MaxR = val[4],
        EE_FinalN = val[5], EE_R = val[6], EE_TB = val[7], EE_I = val[8],
        EE_RinnoTB = val[9], EE_RinTB = val[10])
}

stopCluster(cl)
dfFDacute <- as.data.frame(cbind(dfFD, FDacute))

dfFDacute$change_beta_tu <- dfFDacute$beta_tu/ params$beta_u
dfFDacute$change_alpha_tu <- dfFDacute$alpha_tu/ params$alpha_u

write.csv(dfFDacute, 
    paste(data.dir, "acute_frequencydependent_2019.csv", sep = ""))
rm(params)