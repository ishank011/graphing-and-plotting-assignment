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
        plt.figure()
        plt.title('Lagrange\'s Method')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        x=np.linspace(-20, 20, 500)                                           
        for i in range(n):
            temp=(P.polydiv(w,(-1*L[i],1))[0])/P.polyval(L[i],derivative)
            y=sum([temp[i]*(x**i) for i in range(len(result)-1, -1, -1)])
            if i==0:
                plt.plot(x, y, 'g--', label='Fundamental Polynomials')
            else:
                plt.plot(x, y, 'g--')
            temp*=M[i]
            result+=temp
        s=''
        for i in range(len(result)-1, -1, -1):
            if result[i]>=0:
                s=s+(' +'+str(result[i])+'x^'+str(i))
            else:
                s=s+(' '+str(result[i])+'x^'+str(i))
        s=s+' = 0'
        y=sum([result[i]*(x**i) for i in range(len(result)-1, -1, -1)])
        plt.plot(x, y, 'b', label=s)
        plt.plot(L[0], M[0], 'ro', label='Node Points')
        for i in range(1, len(L)):
            plt.plot(L[i], M[i], 'ro')
        plt.legend(loc='upper right', numpoints=1)
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
        x=np.linspace(-10, 10, 500)
        y=sum([result[i]*(x**i) for i in range(len(result)-1, -1, -1)])
        plt.figure()
        plt.title('Newton\'s Method')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.plot(x, y, 'c', label=s)
        plt.plot(L[0], M[0], 'ro', label='Node Points')
        for i in range(1, len(L)):
            plt.plot(L[i], M[i], 'ro')
        plt.legend(loc='upper right', numpoints=1)
        plt.show()
        return s


#apx=Interpolate()                                                          
#for method in ["newton","lagrange"]:
#    solution=apx.solve([1,2,3,7],[0,-1,0,5],method)
#    print(solution)
