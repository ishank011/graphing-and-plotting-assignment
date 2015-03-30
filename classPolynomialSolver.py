# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

def f(order, coefficients, x):
    s = sum([coefficients[i]*(x**(order-i)) for i in range(order+1)])
    return s

def diff(order, coefficients, x):
    s = sum([coefficients[i]*(order-i)*(x**(order-i-1)) for i in range(order)])
    return s

class PolynomialSolver(object):
    def __init__(self):
        self.mid = None
        self.x0 = None
        self.x1 = None
        self.x2 = None
        
    def BisectionSearch(self, order, coefficients, low, high, epsilon):
        plt.figure()
        plt.title('Bisection Search')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        x=np.linspace(low, high, 500)
        y=f(order, coefficients, x)
        plt.plot(x, y, 'darkslateblue')
        num = 0
        while (high-low)>=epsilon and num<500:
            f_low=f(order, coefficients, low)
            f_high=f(order, coefficients, high)
            plt.plot(low, f_low, 'gbcm'[num%4]+'o', ms=5)
            plt.plot([low, low], [f_low, 0], 'gbcm'[num%4]+'--')
            plt.plot(high, f_high, 'gbcm'[num%4]+'o', ms=5)
            plt.plot([high, high], [f_high, 0], 'gbcm'[num%4]+'--')
            self.mid = (low+high)/2
            if f(order, coefficients, low)*f(order, coefficients, self.mid) > 0:
                low = self.mid
            else:
                high = self.mid
            num += 1
        s='Root of the function x = '+str(self.mid)
        plt.plot(self.mid, f(order, coefficients, self.mid), 'rD', ms=8, label=s)
        plt.legend(loc='upper right', numpoints=1)
        plt.show()
        return self.mid

    def Secant(self, order, coefficients, epsilon):
        num, self.x0, self.x1, self.x2 = 0, 1, 2, 3
        while abs(f(order, coefficients, self.x1))>=epsilon and num<500:
            self.x2, self.x1, self.x0 = (self.x1 - ((self.x1-self.x0)*f(order, coefficients, self.x1))/(f(order, coefficients, self.x1)-f(order, coefficients, self.x0))), self.x2, self.x1
            num += 1
        return self.x1

    def SecantRF(self, order, coefficients, epsilon):
        num, self.x0, self.x1, self.x2 = 0, 1, 2, 3
        while abs(f(order, coefficients, self.x1))>=epsilon and num<500:
            self.x2, self.x1, self.x0 = (self.x1 - ((self.x1-self.x0)*f(order, coefficients, self.x1))/(f(order, coefficients, self.x1)-f(order, coefficients, self.x0))), self.x2, self.x1
            while f(order, coefficients, self.x2)*f(order, coefficients, self.x1) > 0:
                self.x2, self.x1, self.x0 = (self.x1 - ((self.x1-self.x0)*f(order, coefficients, self.x1))/(f(order, coefficients, self.x1)-f(order, coefficients, self.x0))), self.x2, self.x1            
            num += 1
        return self.x1
    
    def NewtonRaphson(self, order, coefficients, epsilon):
        num, self.x0, self.x1 = 0, 1, 2
        while abs(f(order, coefficients, self.x0))>=epsilon and num<500:
            self.x1, self.x0 = (self.x0 - (f(order, coefficients, self.x0)/diff(order, coefficients, self.x0))), self.x1
            num += 1
        return self.x0

    def solve(self, order, coefficients, method, low, high, epsilon = 1e-6):
        if method == 'bisection':
            return self.BisectionSearch(order, coefficients, low, high, epsilon)
        elif method == 'secant':
            return self.Secant(order, coefficients, epsilon)
        elif method == 'secantrf':
            return self.SecantRF(order, coefficients, epsilon)
        elif method == 'newtonraphson':
            return self.NewtonRaphson(order, coefficients, epsilon)

#ps = PolynomialSolver()
#for method in ['bisection', 'secant', 'secantrf', 'newtonraphson'] :
#    solution = ps.solve(3, [5, -6, 11, -8], method, -5, 7, 0.000001)
#    print (solution)
