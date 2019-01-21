# Functions to calculate summary statistics for each model

get_chronic_summary_stats_DD = function(params){
	with(as.list(params),{
		
		# sent chronic free conditions at equilibrium
		R0_TB = (beta_t*K)/(mu+alpha_t)
		EE_TB = K/(b1*2)*( (b*b1-mu-alpha_t)/r - (1+b1)/R0_TB + 
			sqrt(((b*b1-mu-alpha_t)/r - (1+b1)/R0_TB)^2 + 
			(4*b1)/R0_TB*(1-1/R0_TB)))
		EE_suscept = K*1/R0_TB
		
		# Ro 
		Ro = (beta_c * EE_suscept) / (mu + alpha_c + beta_ct * EE_TB ) + 
    		(beta_c * EE_suscept) / (mu + alpha_c + beta_ct * EE_TB) *
    		(beta_ct * EE_TB) / (mu + alpha_ct) + 
    		(beta_tc * EE_TB) / (mu + alpha_tc)
		
		# get EE co-infection conditions
		times <- seq(0, 90000, 1)
		S0 <- EE_suscept-1
		It0 <- EE_TB
		Itc0 <- 0
		Ict0 <- 0
		Ic0 <- 1
		x0 <- c(S = S0, It = It0, Ic = Ic0, Itc = Itc0, Ict = Ict0)
		sol <- as.data.frame(
		    ode(x0, times, rhs_chronicDD, params, method = "ode45"))
		EE_N <- sol$S[length(sol$S)] + sol$Ic[length(sol$S)] + 
			sol$It[length(sol$Itc)] + sol$Itc[length(sol$Itc)] + 
			sol$Ict[length(sol$Ict)]
		EE_TB <- (sol$It[length(sol$Itc)] + sol$Itc[length(sol$Itc)] + 
			sol$Ict[length(sol$Ict)])/ EE_N
		EE_C <- (sol$Ic[length(sol$S)] + sol$Itc[length(sol$Itc)] + 
			sol$Ict[length(sol$Ict)])/ EE_N
		EE_CinS <- (sol$Ic[length(sol$Ic)]) / (sol$Ic[length(sol$Ic)] + 
			sol$S[length(sol$S)])
		EE_CinTB <- (sol$Itc[length(sol$Itc)] + 
			sol$Ict[length(sol$Ict)]) / (sol$It[length(sol$It)] + 
			sol$Itc[length(sol$Itc)] + sol$Ict[length(sol$Ict)])
		return(c(Ro, EE_C, EE_TB, EE_CinS, EE_CinTB, EE_N))	
		}
	)
}

get_chronic_summary_stats_FD = function(params){
	with(as.list(params),{
		# Set chronic free conditions at equilibrium
		R0_TB = (beta_t*K)/(mu+alpha_t)
		EE_TB = K/(b1*2)*((b*b1-mu-alpha_t)/r - (1+b1)/R0_TB +
			sqrt(((b*b1-mu-alpha_t)/r - (1+b1)/R0_TB)^2 +
			(4*b1)/R0_TB*(1-1/R0_TB)))
		EE_suscept = K*1/R0_TB;  
		Total_TB = EE_TB + EE_suscept
		
		# Ro 
		Ro = (EE_suscept / Total_TB) * (beta_c) / 
			(mu + alpha_c + beta_ct * EE_TB ) + 
    		(EE_suscept / Total_TB) * (beta_c) / 
    		(mu + alpha_c + beta_ct * EE_TB) * 
    		(EE_TB) * (beta_ct) / (mu + alpha_ct) +
    		(EE_TB / Total_TB) * (beta_tc) / (mu + alpha_tc)
		
		# get EE co-infection conditions
		times <- seq(0, 400000, 1)
		S0 <- EE_suscept-1; 
		It0 <- EE_TB
		Itc0 <- 0
		Ict0 <- 0
		Ic0 <- 1
		x0 <- c(S = S0, It = It0, Ic = Ic0, Itc = Itc0, Ict = Ict0)
		sol <- as.data.frame(
		    ode(x0, times, rhs_chronicFD, params, method = "ode45"))
		EE_N <- sol$S[length(sol$S)] + sol$Ic[length(sol$S)] + 
			sol$It[length(sol$Itc)] + sol$Itc[length(sol$Itc)] + 
			sol$Ict[length(sol$Ict)]
		EE_TB <- (sol$It[length(sol$Itc)] + sol$Itc[length(sol$Itc)] + 
			sol$Ict[length(sol$Ict)])/ EE_N
		EE_C <- (sol$Ic[length(sol$S)] + sol$Itc[length(sol$Itc)] + 
			sol$Ict[length(sol$Ict)])/ EE_N
		EE_CinS <- (sol$Ic[length(sol$Ic)]) / (sol$Ic[length(sol$Ic)] + 
			sol$S[length(sol$S)])
		EE_CinTB <- (sol$Itc[length(sol$Itc)] + 
			sol$Ict[length(sol$Ict)]) /
			(sol$It[length(sol$It)] + sol$Itc[length(sol$Itc)] + 
			sol$Ict[length(sol$Ict)])
		return(c(Ro, EE_C, EE_TB, EE_CinS, EE_CinTB, EE_N))	
		}
	)
}


get_acute_summary_stats_DD = function(params){
	with(as.list(params),{
		# # Ro
		R0_TB = (beta_t*K)/(mu+alpha_t)
		EE_TB = K/(b1*2)*( (b*b1-mu-alpha_t)/r - 
			(1+b1)/R0_TB + sqrt(((b*b1-mu-alpha_t)/r - 
			(1+b1)/R0_TB)^2 + (4*b1)/R0_TB*(1-1/R0_TB)))
		EE_suscept = K*1/R0_TB;
		R0_aTB = (beta_u*EE_suscept) / (mu+alpha_u+gamma_u) +
			(beta_tu*EE_TB)/(mu+alpha_tu+gamma_tu)
		 
		# get EE co-infection conditions
		times <- seq(0, 90000, 1)
		S0 <- EE_suscept - 1
  		It0 <- EE_TB
  		Itu0 <- 0
  		Rtu0 <- 0
  		Iu0 <- 1
  		Ru0 <- 0
  		x0 <- c(S = S0, It = It0, Itu = Itu0, Rtu = Rtu0, Iu = Iu0, Ru = Ru0)

		sol <- as.data.frame(
		    ode(x0, times, rhs_acuteDD, params, method = "ode45"))
		
		MaxA <- max(sol$Iu + sol$Itu)[1]
		TimeMax<- sol$time[which.max(sol$Iu + sol$Itu)]
		MaxR <- max(sol$Rtu + sol$Ru)[1]
		EE_FinalN <- sol$S[length(sol$S)] + sol$It[length(sol$S)] +
			sol$Itu[length(sol$S)] + sol$Rtu[length(sol$S)] + 
			sol$Iu[length(sol$S)] + sol$Ru[length(sol$S)]
		EE_R <- (sol$Rtu[length(sol$S)] + sol$Ru[length(sol$S)]) / EE_FinalN
		EE_TB<- (sol$It[length(sol$S)] + sol$Itu[length(sol$S)] +
			sol$Rtu[length(sol$S)]) / EE_FinalN
		EE_I <- (sol$Itu[length(sol$S)] + sol$Iu[length(sol$S)]) / EE_FinalN
		EE_RinnoTB <- sol$Ru[length(sol$S)] / (sol$S[length(sol$S)] + 
			sol$Iu[length(sol$S)] + sol$Ru[length(sol$S)])
		EE_RinTB <- sol$Rtu[length(sol$S)] / (sol$It[length(sol$S)] + 
			sol$Itu[length(sol$S)] + sol$Rtu[length(sol$S)])
		
		return(c(R0_aTB, MaxA, TimeMax, MaxR, EE_FinalN, EE_R,
			EE_TB, EE_I, EE_RinnoTB, EE_RinTB))
		}
	)
}



get_acute_summary_stats_FD = function(params){
	with(as.list(params),{
		R0_TB = (beta_t*K)/(mu+alpha_t)
		EE_TB = K/(b1*2)*( (b*b1-mu-alpha_t)/r - 
			(1+b1)/R0_TB + sqrt(((b*b1-mu-alpha_t)/r - 
			(1+b1)/R0_TB)^2 + (4*b1)/R0_TB*(1-1/R0_TB)))
		EE_suscept = K*1/R0_TB
		TotalPop_TB = EE_suscept + EE_TB
		
		# Ro
		R0_aTB = (beta_u*EE_suscept)/(mu+alpha_u+gamma_u)*
			(1/TotalPop_TB) + (beta_tu*EE_TB)/
			(mu+alpha_tu+gamma_tu)*(1/TotalPop_TB)
			
		# get EE co-infection conditions
		times <- seq(0, 90000, 1)
		S0 <- EE_suscept - 1
  		It0 <- EE_TB
  		Itu0 <- 0
  		Rtu0 <- 0
  		Iu0 <- 1
  		Ru0 <- 0
  		x0 <- c(S = S0, It = It0, Itu = Itu0, Rtu = Rtu0, Iu = Iu0, Ru = Ru0)

		sol <- as.data.frame(ode(x0, times, rhs_acuteFD, params, method = "ode45"))
		
		MaxA <- max(sol$Iu + sol$Itu)[1]
		TimeMax<- sol$time[which.max(sol$Iu + sol$Itu)]
		MaxR <- max(sol$Rtu + sol$Ru)[1]
		EE_FinalN <- sol$S[length(sol$S)] + sol$It[length(sol$S)] +
			sol$Itu[length(sol$S)] + sol$Rtu[length(sol$S)] + 
			sol$Iu[length(sol$S)] + sol$Ru[length(sol$S)]
		EE_R <- (sol$Rtu[length(sol$S)] + sol$Ru[length(sol$S)]) / EE_FinalN
		EE_TB<- (sol$It[length(sol$S)] + sol$Itu[length(sol$S)] +
			sol$Rtu[length(sol$S)]) / EE_FinalN
		EE_I <- (sol$Itu[length(sol$S)] + sol$Iu[length(sol$S)]) / EE_FinalN
		EE_RinnoTB <- sol$Ru[length(sol$S)] / (sol$S[length(sol$S)] + 
			sol$Iu[length(sol$S)] + sol$Ru[length(sol$S)])
		EE_RinTB <- sol$Rtu[length(sol$S)] / (sol$It[length(sol$S)] + 
			sol$Itu[length(sol$S)] + sol$Rtu[length(sol$S)])
		
		return(c(R0_aTB, MaxA, TimeMax, MaxR, EE_FinalN, EE_R,
			EE_TB, EE_I, EE_RinnoTB, EE_RinTB))	
		}
	)
}
			