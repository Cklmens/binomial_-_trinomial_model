import numpy as np
import math as ma
import matplotlib.pyplot as plt

from scipy.stats import norm


def black_scholes(S,K,T,r,sigma, style):
    d1 = (np.log(S/K) + (r  + sigma**2/2)*T) / sigma*np.sqrt(T)
    d2 = d1 - sigma* np.sqrt(T)
    if (style=='call'):
        return S * norm.cdf(d1)  - K * np.exp(-r*T)*norm.cdf(d2)
    if(style=='put'):
        return K * np.exp(-r*T)*norm.cdf(-d2)-S * norm.cdf(-d1)

def delta(S,K,T,r,sigma, style):
    d1=(ma.log(S/K)+(r+(sigma**2)/2)*T)/(sigma*ma.sqrt(T))
    if(style=='call'):
       return norm.cdf(d1)
    if(style=='put'):
       return norm.cdf(d1)-1

def gamma(S,K,T,r,sigma,style):
    d1=(ma.log(S/K)+(r+(sigma**2)/2)*T)/(sigma*ma.sqrt(T))
    return (norm.pdf(d1)/(S*sigma*ma.sqrt(T)))

def theta(S,K,T,r,sigma,style):
    d1=(ma.log(S/K)+(r+(sigma**2)/2)*T)/(sigma*ma.sqrt(T))
    d2=d1-sigma*ma.sqrt(T)
    if(style=='call'):
      return ((-(S*sigma*norm.pdf(d1))/(2*ma.sqrt(T)))-K*r*ma.exp(-r*T)*norm.cdf(d2))
    if(style=='put'):
      return (-(S*sigma*norm.pdf(d1))/(2*ma.sqrt(T))+K*r*ma.exp(-r*T)*norm.cdf(-d2))

def rho(S,K,T,r,sigma,style):
    d1=(ma.log(S/K)+(r+(sigma**2)/2)*T)/(sigma*ma.sqrt(T))
    d2=d1-sigma*ma.sqrt(T)
    if(style=='call'):
        return K*T*ma.exp(-r*T)*norm.cdf(d2)
    if(style=='put'):
        return -K*r*ma.exp(-r*T)*norm.cdf(-d2)

def vega(S,K,T,r,sigma, style):
    d1=(ma.log(S/K)+(r+(sigma**2)/2)*T)/(sigma*ma.sqrt(T))
    d2=d1-sigma*ma.sqrt(T)
    return S*norm.pdf(d1)*ma.sqrt(T)



      