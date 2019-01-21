b4 = 0.6
b4L = 0.5
b4U = 0.7 
b5 = 0.9
beta_c = 0.000083074/365.25
beta_ct = 0.0001/365.25
beta_tc = beta_c*2
beta_tcL = beta_c*1
beta_tcU = beta_c*4
alpha_c = 0.001/365.25
alpha_tc = 2*alpha_c
alpha_ct = alpha_tc
alpha_tcL = 1*alpha_c 
alpha_tcU = 4*alpha_c
alpha_ctL = alpha_tcL 
alpha_ctU = alpha_tcU

chronicmildparams = list(
	b4 = b4, b4L = b4L, b4U = b4U, b5 = b5,
	beta_c = beta_c, beat_ct = beta_ct, beta_tc = beta_tc, 
	beta_tcL = beta_tcL, beta_tcU = beta_tcU, alpha_c = alpha_c, 
	alpha_tc = alpha_tc, alpha_ct = alpha_ct, alpha_tcL = alpha_tcL,
	alpha_tcU = alpha_tcU, alpha_ctL = alpha_ctL, alpha_ctU = alpha_ctU)
	

beta_cf = 0.000083074 * K/365.25
beta_tcf = beta_cf*2
beta_tcLf = beta_cf*1
beta_tcUf = beta_cf*4
	
chronicmildparamsFD = list(
	b4 = b4, b4L = b4L, b4U = b4U, b5 = b5,
	beta_c = beta_cf, beat_ct = beta_ct, beta_tc = beta_tcf, 
	beta_tcL = beta_tcLf, beta_tcU = beta_tcUf, alpha_c = alpha_c, 
	alpha_tc = alpha_tc, alpha_ct = alpha_ct, alpha_tcL = alpha_tcL,
	alpha_tcU = alpha_tcU, alpha_ctL = alpha_ctL, alpha_ctU = alpha_ctU)