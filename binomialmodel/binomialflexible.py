import numpy as np
import math as ma
import matplotlib.pyplot as plt
from scipy.stats import norm

class Optionflexible:
   
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

    
    def multiplperiods_adapte(self, p):
        t=self._time/p
        u= ma.exp((self._freerisk-self._volatility**2/2)+self._volatility*ma.sqrt(t))
        d= ma.exp((self._freerisk-self._volatility**2/2)-self._volatility*ma.sqrt(t))
        period= np.ones((p+1))
        for i in range(p+1):
          #  period[i]=self._assetprice*ma.pow(u, p-i)*ma.pow(d, i)
          period[i]=self._assetprice*ma.pow(u, (2*p-i)/2)*ma.pow(d, i/2)
        return period
    
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
    

    
    def _payoff(self,bool, so, k):
     if (bool==True):
         if (so-k)>0:
             return so-k
         return 0
     else:
         if(k-so)>0:
             return k-so
         return 0
 
    def payoffmultipl(self,p,l):
        cp=self._style
        pay_off= np.ones((p+1,p+1))
        t=self._time/p
        u= ma.exp(self._volatility*ma.sqrt(t)+l*(self._volatility*ma.sqrt(t))**2)
        d= ma.exp(-self._volatility*ma.sqrt(t)+l*(self._volatility*ma.sqrt(t))**2)
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

    def Delta(self,p,l):
        h=0.1
        if(self._style):
            op1=Optionflexible(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "call")
            op2=Optionflexible(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "call")
        else:
            op1=Optionflexible(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "put")
            op2=Optionflexible(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "put")
        return (op1.payoffmultipl(p,l)-op2.payoffmultipl(p,l))/(2*h)
    
    def Gamma(self,p,l):
        h=0.2
        if(self._style):
            op1=Optionflexible(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "call")
            op2=Optionflexible(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "call")
        else:
            op1=Optionflexible(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility, "put")
            op2=Optionflexible(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility, "put")
        return (op1.Delta(p,l)-op2.Delta(p,l))/(2*h)


    def Theta(self,p,l):
        h=0.1
        if(self._style):
            op1=Optionflexible(self._assetprice, self._Strike, self._time+h, self._freerisk,self._volatility, "call")
            op2=Optionflexible(self._assetprice, self._Strike, self._time-h, self._freerisk,self._volatility, "call")
        else:
            op1=Optionflexible(self._assetprice, self._Strike, self._time+h, self._freerisk,self._volatility, "put")
            op2=Optionflexible(self._assetprice, self._Strike, self._time-h, self._freerisk,self._volatility, "put")
        return (op1.payoffmultipl(p,l)-op2.payoffmultipl(p,l))/(2*h)

    def vega(self,p,l):
        h=0.001
        if(self._style):
            op1=Optionflexible(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility+h, "call")
            op2=Optionflexible(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility-h, "call")
        else:
            op1=Optionflexible(self._assetprice+h, self._Strike, self._time, self._freerisk,self._volatility+h, "put")
            op2=Optionflexible(self._assetprice-h, self._Strike, self._time, self._freerisk,self._volatility-h, "put")
        return (op1.payoffmultipl(p,l)-op2.payoffmultipl(p,l))/(2*h)
    
    def rho(self,p,l):
        h=0.001
        if(self._style):
            op1=Optionflexible(self._assetprice+h, self._Strike, self._time, self._freerisk+h,self._volatility, "call")
            op2=Optionflexible(self._assetprice-h, self._Strike, self._time, self._freerisk-h,self._volatility, "call")
        else:
         op1=Optionflexible(self._assetprice+h, self._Strike, self._time, self._freerisk+h,self._volatility, "put")
         op2=Optionflexible(self._assetprice-h, self._Strike, self._time, self._freerisk-h,self._volatility, "put")
        return (op1.payoffmultipl(p,l)-op2.payoffmultipl(p,l))/(2*h)

class Drawpricewalk(Optionflexible):
    
    def __init__(self,S,K,T,r,v,style):
        super().__init__(S,K,T,r,v,style)
    
    def multiplperiods(self, p,l):
        t=self._time/p
        u= ma.exp(self._volatility*ma.sqrt(t)+l*(self._volatility*ma.sqrt(t))**2)
        d= ma.exp(-self._volatility*ma.sqrt(t)+l*(self._volatility*ma.sqrt(t))**2)

        temp=[]
        for j in range(p):
          period= []
          for i in range(2*j+1):
            period.append(round(self._assetprice*ma.pow(u, j)*ma.pow(d, i),2))
          temp.append(period)
        return temp

    def Convergence(self,p,l):
        err=[]
        bs=self.black_scholes()
        for i in range(1,p+1):
            err.append(bs-self.payoffmultipl(i,l))

        plt.plot(np.arange(0,p,1), err)
        plt.xlabel("Période")
        plt.ylabel("Err=C_bs-C")
        plt.title("Convergence numérique du modèle binomial généralisé l=2")
    
    def prices(self, p,l):
        t=self._time/p
        u= ma.exp(self._volatility*ma.sqrt(t)+l*(self._volatility*ma.sqrt(t))**2)
        d= ma.exp(-self._volatility*ma.sqrt(t)+l*(self._volatility*ma.sqrt(t))**2)
        period= np.ones((2*p+1))
        for i in range(2*p+1):
            period[i]=self._assetprice*ma.pow(u, p)*ma.pow(d, i)
        return period
    
    def plottree(self,p,l):
        t=self._time/p
        fig = plt.figure(figsize=[5, 5])
        ordonnee=self.multiplperiods(p,l)
        for i in range(p-1):
            x = [1, 0, 1]
            for j in range(i):
                x.append(0)
                x.append(1)
            x = np.array(x) + i
            y=ordonnee[i+1]
            plt.plot(x, y, 'bo-')
        plt.xlabel("Période")
        plt.ylabel("Prix")
        plt.title("Arbre binomial ")
        plt.show()

   

