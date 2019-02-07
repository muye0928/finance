## Monte Carlo valuation of European call option in Black-Scholes-Merton model#
import numpy as np
# parameter numbers
S0 =100 #initial index level
K = 105 # strike price
T = 1.0 # time to maturity
r = 0.05 # riskless short rate
sigma = 0.2 # volatility

I = 10000 # number of simulations

# Valuation Algorithm
z = np.random.standard_normal(I) # pseudorandom numbers
ST = S0 * np.exp((r - 0.5 * sigma ** 2) * T + sigma * np.sqrt(T) * z)
# index values at maturity
hT = np.maximum(ST - K, 0) # inner values at maturity
C0 = np.exp(-r * T) * np.sum(hT) / I # Monte Carlo estimator

# Result Output
print ("Value of the European Call Option %5.3f" % C0)



## Valuation of European call options in Black-Scholes-Merton model
# incl. Vega function and implied volatility estimation

# Analytical Black-Scholes-Merton (BSM) Formula
def bsm_call_value(S0,K,T,r,sigma):
    "valuation of European call options in bsm model.Analytics formula."

    from math import sqrt, log, exp
    from scipy import stats

    S0 = float(S0)
    d1 = (log(S0 / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * sqrt(T))
    d2 = (log(S0 / K) + (r - 0.5 * sigma ** 2) * T) / (sigma * sqrt(T))
    value = (S0 * stats.norm.cdf(d1, 0.0, 1.0)- K * exp(-r * T) * stats.norm.cdf(d2, 0.0, 1.0))
    return value

# Vega function
def bsm_vega(S0, K, T, r, sigma):
    "Vega of European  options in BSM model.(Vega:float partial derivative of BSM formula with respect to sigma)"
    from math import sqrt, log
    from scipy import stats
    S0 = float(S0)
    d1 = (log(S0 / K) + (r + 0.5 * sigma ** 2) * T / (sigma * sqrt(T))
    vega = S0 * stats.norm.cdf(d1,0.0,1.0) * sqrt(T)
    return vega


# Implied volatility function
def bsm_call_imp_vol(S0, K, T, r, C0, sigma_est, it=100):
    "Implied volatility of European call option in BSM model."
    for i in range(it):
        sigma_est -= ((bsm_call_value((S0,K,T,r,sigma_set) - CO))/bsm_vega(S0,K,T,r,sigma_est))
    return sigma_est







