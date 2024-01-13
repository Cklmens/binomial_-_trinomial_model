import numpy as np
import math as ma
import matplotlib.pyplot as plt
from scipy.stats import norm

class OptionAmericaine:
    def __init__(self,S,K,T,r,v,style):
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
        temp=[]
        for j in range(p):
          period= []
          for i in range(j+1):
            period.append(round(self._assetprice*ma.pow(u, j-i)*ma.pow(d, i),2))
          temp.append(period)
        return temp
    
    def multiplperiods_adapte(self, p):
        t=self._time/p
        u= ma.exp((self._freerisk-self._volatility**2/2)+self._volatility*ma.sqrt(t))
        d= ma.exp((self._freerisk-self._volatility**2/2)-self._volatility*ma.sqrt(t))
        period= np.ones((p+1))
        for i in range(p+1):
            period[i]=self._assetprice*ma.pow(u, p-i)*ma.pow(d, i)
        return period
    
    def actualisePayoff(self,p):
        tmp= self.multiplperiods(p)
        if(self._style):
          for i in range(len(tmp)):
            for j in range(len(tmp[i])):
                tmp[i][j]= round(max(tmp[i][j]-self._Strike,0),3)#*ma.exp(-self._freerisk*self._time/p)
        else:
          for i in range(len(tmp)):
            for j in range(len(tmp[i])):
                tmp[i][j]= round(max(self._Strike-tmp[i][j],0),3)#*ma.exp(-self._freerisk*self._time/p)
        return tmp
    
    def _payoff(self,bool, so, k):
     if (bool==True):
         if (so-k)>0:
             return so-k
         return 0
     else:
         if(k-so)>0:
             return k-so
         return 0

    def payoffmultipl(self,p):
        cp=self._style
        t=self._time/p
        u= ma.exp(self._volatility*ma.sqrt(t))
        d= 1/u
        po=(ma.exp(self._freerisk*t)-d)/(u-d)
        period=self.actualisePayoff(p)
        payoff=[]
        payoff.append(period[-1])
    
        for j in range (p-1):
            tmp=[]
            for i in range (len(payoff[-1])-1):
                tmp.append(max((po*payoff[-1][i]+(1-po)*payoff[-1][i+1])*ma.exp(-self._freerisk*t), period[p-j-2][i] ))#period[len(period)-j-2][i]
            payoff.append(tmp)       
        return payoff[-1][0]
    def Delta(self,p):
        h=0.1
        if(self._style):
            op1=OptionAmericaine(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "call")
            op2=OptionAmericaine(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "call")
        else:
            op1=OptionAmericaine(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "put")
            op2=OptionAmericaine(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "put")
        return (op1.payoffmultipl(p)-op2.payoffmultipl(p))/(2*h)
    
    def Gamma(self,p):
        h=0.2
        if(self._style):
            op1=OptionAmericaine(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "call")
            op2=OptionAmericaine(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "call")
        else:
            op1=OptionAmericaine(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "put")
            op2=OptionAmericaine(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "put")
        return (op1.Delta(p)-op2.Delta(p))/(2*h)


    def Theta(self,p):
        h=0.1
        if(self._style):
            op1=OptionAmericaine(self._assetprice, self._Strike, self._time+h, self._freerisk,self._volatility, "call")
            op2=OptionAmericaine(self._assetprice, self._Strike, self._time-h, self._freerisk,self._volatility, "call")
        else:
            op1=OptionAmericaine(self._assetprice, self._Strike, self._time+h, self._freerisk,self._volatility, "put")
            op2=OptionAmericaine(self._assetprice, self._Strike, self._time-h, self._freerisk,self._volatility, "put")
        return (op1.payoffmultipl(p)-op2.payoffmultipl(p))/(2*h)

    def vega(self,p):
        h=0.001
        if(self._style):
            op1=OptionAmericaine(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility+h, "call")
            op2=OptionAmericaine(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility-h, "call")
        else:
            op1=OptionAmericaine(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility+h, "put")
            op2=OptionAmericaine(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility-h, "put")
        return (op1.payoffmultipl(p)-op2.payoffmultipl(p))/(2*h)
    
    def rho(self,p):
        h=0.001
        if(self._style):
            op1=OptionAmericaine(self._assetprice+h, self._Strike, self._time, self._freerisk+h,self._volatility, "call")
            op2=OptionAmericaine(self._assetprice-h, self._Strike, self._time, self._freerisk-h,self._volatility, "call")
        else:
         op1=OptionAmericaine(self._assetprice+h, self._Strike, self._time, self._freerisk+h,self._volatility, "put")
         op2=OptionAmericaine(self._assetprice-h, self._Strike, self._time, self._freerisk-h,self._volatility, "put")
        return (op1.payoffmultipl(p)-op2.payoffmultipl(p))/(2*h)

