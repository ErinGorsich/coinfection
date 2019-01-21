beta_t = 0.0001/365.25; #baseline transmission rate for TB, 30% prevalence
    #for 5% TB let beta_t = 0.0000715/365.25
    #for 16% tB let beta_t = 0.0000815/365.25;
    #for 27% TB let beta_t = 0.000094/365.25;
    #for 38% TB let beta_t = 0.000113/365.25;
    #for 49% TB let beta_t = 0.000139/365.25;
    #for 60% TB let beta_t = 0.00018/365.25;
b = 0.37/365.25
b1 = 0.7
mu = (1/15)/365.25 
K = 1000
alpha_t = 0.001/365.25
r = b-mu

fixedparams = list(beta_t = beta_t, b = b, b1 = b1,
	mu = mu, K = K, alpha_t = alpha_t, r = r)
