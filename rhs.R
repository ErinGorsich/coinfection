rhs_acuteFD = function(times, x, params){
	with(as.list(c(x, params)), {
		S <- x[1]; It <- x[2]; Itu <- x[3];
		Rtu <- x[4]; Iu <- x[5]; Ru <- x[6]
		
		N <- S + It + Itu + Rtu + Iu + Ru
		Nb <- S + b1*It + b2*Itu + b1*Rtu + b3*Iu + Ru
  
  		#Calculate derivates
  		dSdt = b*Nb*(1-((r)/b)*(N/K))- mu*S - beta_t*(It + Itu + Rtu)*S - 
  		    beta_u*(Itu + Iu)*S/N
  		dItdt = beta_t*(It + Itu + Rtu)*S - beta_tu*(Itu + Iu)*It/N - mu*It - 
  		    alpha_t*It
  		dItudt = beta_tu*(Itu + Iu)*It/N - gamma_tu*Itu  - mu*Itu - alpha_tu*Itu
  		dRtudt = gamma_tu*Itu + beta_t*(It + Itu + Rtu)*Ru - mu*Rtu - 
  		    alpha_t*Rtu
  		dIudt = beta_u*(Itu + Iu)*S/N- mu*Iu - alpha_u*Iu - gamma_u*Iu
  		dRudt = gamma_u*Iu - mu*Ru - beta_t*(It + Itu + Rtu)*Ru
		
		out = list(c(dSdt, dItdt, dItudt, dRtudt, dIudt, dRudt))
		return(out)
		}
	)
}

rhs_chronicFD = function(times, x, params){
	with(as.list(c(x, params)), {		
		S <- x[1]; It <- x[2]; Ic <- x[3]
		Itc <- x[4]; Ict <- x[5]

		N = S + It + Ic + Itc + Ict;
    		Nb = S + b1*It + b4*Itc + b5*Ic + b4*Ict;
  
  		#Calculate derivates
  		dSdt = b*Nb*(1-((r)/b)*(N/K)) - mu*S - beta_t*(It+Itc+Ict)*S -
		    beta_c*(Ic+Itc+Ict)*S/N
   		dItdt = beta_t*(It+Itc+Ict)*S - beta_tc*(Ic+Itc+Ict)*It/N - mu*It -
   		    alpha_t*It
   		dIcdt = beta_c*(Ic+Itc+Ict)*S/N - beta_ct*(It+Itc+Ict)*Ic - mu*Ic - 
   		    alpha_c*Ic
   		dItcdt = beta_tc*(Ic+Itc+Ict)*It/N - mu*Itc - alpha_tc*Itc
   		dIctdt = beta_ct*(It+Itc+Ict)*Ic -mu*Ict - alpha_ct*Ict
      		
		out = list(c(dSdt, dItdt, dIcdt, dItcdt, dIctdt))
		return(out)
		}
	)
}


rhs_acuteDD = function(times, x, params){
	with(as.list(c(x, params)), {		
		S <- x[1]; It <- x[2]; Itu <- x[3];
		Rtu <- x[4]; Iu <- x[5]; Ru <- x[6]
		
		N <- S + It + Itu + Rtu + Iu + Ru
		Nb <- S + b1*It + b2*Itu + b1*Rtu + b3*Iu + Ru
  
  		#Calculate derivates
  		dSdt = b*Nb*(1-((r)/b)*(N/K))- mu*S - beta_t*(It + Itu + Rtu)*S -
  		    beta_u*(Itu + Iu)*S 
  		dItdt = beta_t*(It + Itu + Rtu)*S - beta_tu*(Itu + Iu)*It - mu*It -
  		    alpha_t*It
  		dItudt = beta_tu*(Itu + Iu)*It - gamma_tu*Itu  - mu*Itu - alpha_tu*Itu
  		dRtudt = gamma_tu*Itu + beta_t*(It + Itu + Rtu)*Ru - mu*Rtu - 
  		    alpha_t*Rtu
  		dIudt = beta_u*(Itu + Iu)*S - mu*Iu - alpha_u*Iu - gamma_u*Iu
  		dRudt = gamma_u*Iu - mu*Ru - beta_t*(It + Itu + Rtu)*Ru
		
		out = list(c(dSdt, dItdt, dItudt, dRtudt, dIudt, dRudt))
		return(out)
		}
	)
}

rhs_chronicDD = function(times, x, params){
	with(as.list(c(x, params)), {		
		S <- x[1]; It <- x[2]; Ic <- x[3]
		Itc <- x[4]; Ict <- x[5]

		N = S + It + Ic + Itc + Ict;
    		Nb = S + b1*It + b4*Itc + b5*Ic + b4*Ict;
  
  		#Calculate derivates
  		dSdt = b*Nb*(1-((r)/b)*(N/K)) - mu*S -beta_t*(It+Itc+Ict)*S - 
		    beta_c*(Ic+Itc+Ict)*S
   		dItdt = beta_t*(It+Itc+Ict)*S - beta_tc*(Ic+Itc+Ict)*It - mu*It - 
   		    alpha_t*It
   		dIcdt = beta_c*(Ic+Itc+Ict)*S - beta_ct*(It+Itc+Ict)*Ic - mu*Ic -
   		    alpha_c*Ic
   		dItcdt = beta_tc*(Ic+Itc+Ict)*It - mu*Itc - alpha_tc*Itc
   		dIctdt = beta_ct*(It+Itc+Ict)*Ic -mu*Ict - alpha_ct*Ict
   		
		out = list(c(dSdt, dItdt, dIcdt, dItcdt, dIctdt))
		return(out)
		}
	)
}



