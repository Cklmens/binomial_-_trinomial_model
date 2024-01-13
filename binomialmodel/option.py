import numpy as np
import math as ma
import matplotlib.pyplot as plt
from scipy.stats import norm
class Option:
   
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
        period= np.ones((p+1))
        for i in range(p+1):
            period[i]=self._assetprice*ma.pow(u, p-i)*ma.pow(d, i)
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
 
    def payoffmultipl(self,p):
        cp=self._style
        pay_off= np.ones((p+1,p+1))
        t=self._time/p
        u= ma.exp(self._volatility*ma.sqrt(t))
        d= 1/u
        po=(ma.exp(self._freerisk*t)-d)/(u-d)

        period= np.ones((p+1))
        for i in range(p+1):
            period[i]=self._assetprice*ma.pow(u, p-i)*ma.pow(d, i)
        
        for i in range(p+1):
            pay_off[0][i]=self._payoff(bool(cp),period[i],self._Strike) 
    
        for j in range (1,p+1):
            for i in range (p+1-j):
                pay_off[j][i]=(po*pay_off[j-1][i]+(1-po)*pay_off[j-1][i+1])*ma.exp(-self._freerisk*t)         
        return pay_off[-1][0]
    
    def Delta(self,p):
        h=0.1
        if(self._style):
            op1=Option(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "call")
            op2=Option(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "call")
        else:
            op1=Option(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "put")
            op2=Option(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "put")
        return (op1.payoffmultipl(p)-op2.payoffmultipl(p))/(2*h)
    
    def Gamma(self,p):
        h=0.2
        if(self._style):
            op1=Option(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "call")
            op2=Option(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "call")
        else:
            op1=Option(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "put")
            op2=Option(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "put")
        return (op1.Delta(p)-op2.Delta(p))/(2*h)


    def Theta(self,p):
        h=0.1
        if(self._style):
            op1=Option(self._assetprice, self._Strike, self._time+h, self._freerisk,self._volatility, "call")
            op2=Option(self._assetprice, self._Strike, self._time-h, self._freerisk,self._volatility, "call")
        else:
            op1=Option(self._assetprice, self._Strike, self._time+h, self._freerisk,self._volatility, "put")
            op2=Option(self._assetprice, self._Strike, self._time-h, self._freerisk,self._volatility, "put")
        return (op1.payoffmultipl(p)-op2.payoffmultipl(p))/(2*h)

    def vega(self,p):
        h=0.001
        if(self._style):
            op1=Option(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility+h, "call")
            op2=Option(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility-h, "call")
        else:
            op1=Option(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility+h, "put")
            op2=Option(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility-h, "put")
        return (op1.payoffmultipl(p)-op2.payoffmultipl(p))/(2*h)
    
    def rho(self,p):
        h=0.001
        if(self._style):
            op1=Option(self._assetprice+h, self._Strike, self._time, self._freerisk+h,self._volatility, "call")
            op2=Option(self._assetprice-h, self._Strike, self._time, self._freerisk-h,self._volatility, "call")
        else:
         op1=Option(self._assetprice+h, self._Strike, self._time, self._freerisk+h,self._volatility, "put")
         op2=Option(self._assetprice-h, self._Strike, self._time, self._freerisk-h,self._volatility, "put")
        return (op1.payoffmultipl(p)-op2.payoffmultipl(p))/(2*h)


class Drawpricewalk(Option):
    
    def __init__(self,S,K,T,r,v,style):
        super().__init__(S,K,T,r,v,style)
    
    def multiplperiods(self, p):
        t=self._time/p
        u= ma.exp(self._volatility*ma.sqrt(t))
        d= 1/u
        period= [[self._assetprice]]
        for j in range(p):
            temp=[]
            for i in period[-1]: 
                if(temp.count(round(i*u,2))==0):
                  temp.append(round(i*u,2)) 
                if(temp.count(round(i*d,2))==0):  
                  temp.append(round(i*d,2))
            period.append(temp)         
        return period     

    def Convergence(self,p):
        err=[]
        bs=self.black_scholes()
        for i in range(1,p+1):
            err.append(bs-self.payoffmultipl(i))

        plt.plot(np.arange(0,p,1), err)
        plt.xlabel("Période")
        plt.ylabel("Err=C_bs-C")
        plt.title("Convergence numérique du modèle binomial")
    
    def  black_scholes(self):
        S=self._assetprice
        K=self._Strike
        T=self._time
        r=self._freerisk
        sigma=self._volatility
        style=self._style
        d1 = (np.log(S/K) + (r  + sigma**2/2)*T) / sigma*np.sqrt(T)
        d2 = d1 - sigma* np.sqrt(T)
        if (style==True):
            return S * norm.cdf(d1)  - K * np.exp(-r*T)*norm.cdf(d2)
        if(style==False):
            return K * np.exp(-r*T)*norm.cdf(-d2)-S * norm.cdf(-d1)
            
    def plottree(self,p,):
        t=self._time/p
        u= ma.exp(self._volatility*ma.sqrt(t))
        d= 1/u
        fig = plt.figure(figsize=[5, 5])
        for i in range(p):
            x = [1, 0, 1]
            for j in range(i):
                x.append(0)
                x.append(1)
            x = np.array(x) + i
            #y = self._assetprice+np.arange(-(i+1), i+2)[::-1]*u
            y=self._assetprice+np.concatenate((np.arange(0, i+2)[::-1]*u, np.arange(-(i+1), 0)[::-1]*d))
            plt.plot(x, y, 'bo-')
        plt.show()



#S=100       # So le prix du sous-jacent a t=0
#sigma= 0.25  # sigma: la volatilité constante
#K=110        # Strike
#r=0.05
#T=1

#opt1=Option(S,K,T,r,sigma,"call")
#opt1=Option(S,K,T,r,sigma,"call")
##print(opt1.payoffmultipl(10))
#print(opt1.Delta(100))
#print(opt1.Gamma(100))
#print(opt1.Theta(100))
#print(opt1.vega(100))
#print(opt1.rho(100))
#draw1=Drawpricewalk(S,K,T,r,sigma,"call")
#draw1.plottree(20)