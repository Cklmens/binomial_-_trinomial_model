import numpy as np
import math as ma


def multiplperiods(so,k,T,r,v,p):
    t=T/p
    u= ma.exp(v*ma.sqrt(t))
    d= 1/u
    period= np.ones((2*p+1))
    for i in range(2*p+1):
         period[i]=so*ma.pow(u, p)*ma.pow(d, i)
    return period

def payoff(bool, so, k):
     if (bool==True):
         if (so-k)>0:
             return so-k
         return 0
     else:
         if(k-so)>0:
             return k-so
         return 0
 
def payoffmultipl(so,k,T,r,v,p, cp):
    pay_off= np.ones((p+1,p+1))
    t=T/p
    u= ma.exp(v* ma.sqrt(t))
    d= 1/u
    po=(ma.exp(r*t)-d)/(u-d)

    period= np.ones((2*p+1))
    for i in range(2*p+1):
         period[i]=so*ma.pow(u, p)*ma.pow(d, i)
    
    for i in range(p+1):
         pay_off[0][i]=payoff(bool(cp),period[2*i],k) 
   
    for j in range (1,p+1):
        for i in range (p+1-j):
            pay_off[j][i]=(po*pay_off[j-1][i]+(1-po)*pay_off[j-1][i+1])*ma.exp(-r*t)         
    return pay_off

#print(multiplperiods(30,30,1,0.05,0.2,5))
#print(multiplperiods(30,32,0.5,0.5,0.3,3))
print(payoffmultipl(100,110,1,0.05,0.25,100,1))
#print(payoffmultipl(30,28,1,0.5,0.3,5,0))

