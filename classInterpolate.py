# -*- coding: utf-8 -*-

import numpy as np
from numpy.polynomial import polynomial as P
import matplotlib.pyplot as plt
        
class Interpolate:
    
    def solve(self,L,M,method):
        if(method=="newton"):
            return (self.Newton(L,M))
        else:
            return (self.Lagrange(L,M))

    def Lagrange(self,L,M):                                                
        n=len(L)                                                           
        w=(-1*L[0],1)                                                      
        for i in range(1,n):
            w=P.polymul(w,(-1*L[i],1))                                    
        result=np.array([0.0 for i in range(len(w)-1)])                    
        derivative=P.polyder(w)                                             
        for i in range(n):
            result+=(P.polydiv(w,(-1*L[i],1))[0]*M[i])/P.polyval(L[i],derivative)   
        s=''
        for i in range(len(result)-1, -1, -1):
            if result[i]>=0:
                s=s+(' +'+str(result[i])+'x^'+str(i))
            else:
                s=s+(' '+str(result[i])+'x^'+str(i))
        s=s+' = 0'
        x=np.linspace(-50, 50, 500)
        y=sum([x**i for i in range(len(result)-1, -1, -1)])
        plt.figure()
        plt.title('Lagrange Method')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.text(-25, max(y), s)
        plt.plot(x, y, 'b')
        plt.show()
        return s
                                             
    def Newton(self,L,M):
        n=len(L)                                                            
        mat=[[0.0 for i in range(n)] for j in range(n)]                    
        for i in range(n):                                                 
            mat[i][0]=M[i]
        for i in range(1,n):                                               
            for j in range(n-i):
                mat[j][i]=(mat[j+1][i-1]-mat[j][i-1])/(L[j+i]-L[j])
        result=np.array((mat[0][0],))                                          
        for i in range(1,n):
            prod=(-1*L[0],1)                                               
                                                                            
            for j in range(1,i):
                prod=P.polymul(prod,(-1*L[j],1))                              
            result=P.polyadd(result,np.array(prod)*mat[0][i])
        s=''
        for i in range(len(result)-1, -1, -1):
            if result[i]>=0:
                s=s+(' +'+str(result[i])+'x^'+str(i))
            else:
                s=s+(' '+str(result[i])+'x^'+str(i))
        s=s+' = 0'
        x=np.linspace(-50, 50, 500)
        y=sum([x**i for i in range(len(result)-1, -1, -1)])
        plt.figure()
        plt.title('Newton Method')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.text(-25, max(y), s)
        plt.plot(x, y, 'g')
        plt.show()
        return s


apx=Interpolate()                                                          
for method in ["newton","lagrange"]:
    solution=apx.solve([1,2,3],[0,-1,0],method)
    print(solution)
