import numpy as np
import math as ma
import matplotlib.pyplot as plt


class OptionTrinomial():
    def __init__(self, S, K, T, r, v, style):
        if(S>0):
            self._assetprice=S
        else:
            print("Le sous-jacent ne peut pas etre négatif")
        
        if(K>0):
            self._Strike=K
        else:
             print("Le strike ne peut pas etre négatif")
        
        if(T>0):
            self._time=T
        else:
             print("Le Temps ne peut pas etre négatif")
        if(r>=0 and r<1):
            self._freerisk=r
        else:
             print("Le taux d'intêret sans risque doit etre compris entre 0 et 1")
        
        if(v<=1 and v>=0):
            self._volatility=v
        else:
             print("La volatilité ne peut pas etre négatif")
        if(style=="put"):
            self._style=False
        else:
            if(style=="call"):
               self._style=True
            else:
               print("error")

    def multiplperiods(self, p):
        t=self._time/p
        u= ma.exp(self._volatility*ma.sqrt(t))
        d= 1/u
        m=ma.sqrt(u*d)
        period= np.ones((2*p+1))
        for i in range(2*p+1):
            period[i]=self._assetprice*ma.pow(u, (2*p-i)/2)*ma.pow(d, i/2)
        return period
    
    def _payoff(self,bool, so, k):
     if (bool==True):
         if (so-k)>0:
             return so-k
         return 0
     else:
         if(k-so)>0:
             return k-so
         return 0
    
    def payoffmultipl(self,p,po):
        cp=self._style
        pay_off= np.ones((p+1,2*p+1))
        t=self._time/p
        h=(1-po)/2
        
        u= self._volatility*ma.sqrt(t/(2*h))
        d= 1/u

        po=1-2*h
        pu=(ma.exp(self._freerisk*t)-ma.exp(-u))/(ma.exp(u)-ma.exp(-u))-po*(1-ma.exp(-u))/(ma.exp(u)-ma.exp(-u))
        pd=(ma.exp(u)-ma.exp(self._freerisk*t))/(ma.exp(u)-ma.exp(-u))-po*(ma.exp(u)-1)/(ma.exp(u)-ma.exp(-u))

        period= np.zeros((2*p+1))
        for i in range(2*p+1):
            period[i]=self._assetprice*ma.exp(u*p)*ma.exp(-u*i)
        
        for i in range(2*p+1):
            pay_off[0][i]=self._payoff(bool(cp),period[i],self._Strike) 
    
        for j in range (1,p+1):
            for i in range (2*p+1-j*2):
                pay_off[j][i]=(pu*pay_off[j-1][i]+po*pay_off[j-1][i+1] +pd*pay_off[j-1][i+2])*ma.exp(-self._freerisk*t)   

        return pay_off[-1][0]
    def Delta(self,p,po):
        h=0.1
        if(self._style):
            op1=OptionTrinomial(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "call")
            op2=OptionTrinomial(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "call")
        else:
            op1=OptionTrinomial(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "put")
            op2=OptionTrinomial(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "put")
        return (op1.payoffmultipl(p,po)-op2.payoffmultipl(p,po))/(2*h)
    
    def Gamma(self,p,po):
        h=0.2
        if(self._style):
            op1=OptionTrinomial(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "call")
            op2=OptionTrinomial(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "call")
        else:
            op1=OptionTrinomial(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "put")
            op2=OptionTrinomial(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "put")
        return (op1.Delta(p,po)-op2.Delta(p,po))/(2*h)


    def Theta(self,p,po):
        h=0.1
        if(self._style):
            op1=OptionTrinomial(self._assetprice, self._Strike, self._time+h, self._freerisk,self._volatility, "call")
            op2=OptionTrinomial(self._assetprice, self._Strike, self._time-h, self._freerisk,self._volatility, "call")
        else:
            op1=OptionTrinomial(self._assetprice, self._Strike, self._time+h, self._freerisk,self._volatility, "put")
            op2=OptionTrinomial(self._assetprice, self._Strike, self._time-h, self._freerisk,self._volatility, "put")
        return (op1.payoffmultipl(p,po)-op2.payoffmultipl(p,po))/(2*h)

    def vega(self,p,po):
        h=0.001
        if(self._style):
            op1=OptionTrinomial(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility+h, "call")
            op2=OptionTrinomial(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility-h, "call")
        else:
            op1=OptionTrinomial(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility+h, "put")
            op2=OptionTrinomial(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility-h, "put")
        return (op1.payoffmultipl(p,po)-op2.payoffmultipl(p,po))/(2*h)
    
    def rho(self,p,po):
        h=0.001
        if(self._style):
            op1=OptionTrinomial(self._assetprice+h, self._Strike, self._time, self._freerisk+h,self._volatility, "call")
            op2=OptionTrinomial(self._assetprice-h, self._Strike, self._time, self._freerisk-h,self._volatility, "call")
        else:
         op1=OptionTrinomial(self._assetprice+h, self._Strike, self._time, self._freerisk+h,self._volatility, "put")
         op2=OptionTrinomial(self._assetprice-h, self._Strike, self._time, self._freerisk-h,self._volatility, "put")
        return (op1.payoffmultipl(p,po)-op2.payoffmultipl(p,po))/(2*h)



