b2 = 0.6
b3 = 0.9
beta_u = 0.15/365.25
beta_tu = 2 * beta_u
beta_tuL = 1 * beta_u
beta_tuU = 4 * beta_u
alpha_u = 0.001
alpha_tu = 0.002
alpha_tuL = alpha_u
alpha_tuU = 4*alpha_u
gamma_u = 1/3
gamma_tu = 1/3
gamma_tuL = (1/4)*gamma_u
gamma_tuU = gamma_u

acutemildhightrans = list(
	b2 = b2, b3 = b3, beta_u = beta_u, beta_tu = beta_tu,
	beta_tuL = beta_tuL, beta_tuU = beta_tuU,
	alpha_u = alpha_u, alpha_tu = alpha_tu,
	alpha_tuL = alpha_tuL, alpha_tuU = alpha_tuU,
	gamma_u = gamma_u, gamma_tu = gamma_tu, 
	gamma_tuL = gamma_tuL, gamma_tuU = gamma_tuU)

# Frequency Dependent Transmisison
beta_uFD = (0.15/365.25)*1000
beta_tuFD = beta_uFD*2
beta_tuLFD = beta_uFD*1; 
beta_tuUFD = beta_uFD*4
	
acutemildhightransFD = list(
	b2 = b2, b3 = b3, beta_u = beta_uFD, beta_tu = beta_tuFD,
	beta_tuL = beta_tuLFD, beta_tuU = beta_tuUFD, 
	alpha_u = alpha_u, alpha_tu = alpha_tu, 
	alpha_tuL = alpha_tuL, alpha_tuU = alpha_tuU,
	gamma_u = gamma_u, gamma_tu = gamma_tu,
	gamma_tuL = gamma_tuL, gamma_tuU = gamma_tuU)