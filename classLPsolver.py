# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

class LPsolver():

    def solve(self,a,b,c):
        if len(a)==2:
            plt.figure()
            plt.xlabel('x')
            plt.ylabel('y')
            fill_x=[0]
            fill_y=[0]
            for i in range(len(b)):
                if c[i][0]==0:
                    plt.axhline(y=b[i]/c[i][1], color = 'gbrycm'[i], label='Constraint: '+str(i+1))
                elif c[i][1]==0:
                    plt.axvline(x=b[i]/c[i][0], color = 'gbrycm'[i], label='Constraint: '+str(i+1))
                else:
                    x = np.array([b[i]/c[i][0], 0])
                    y = np.array([0, b[i]/c[i][1]])
                    plt.plot(x, y, 'gbrycm'[i]+'-', label='Constraint: '+str(i+1))
                    
        l_b=len(b)
        l_a=len(a)
        arr=np.array(c)
        arr=np.concatenate((arr,np.identity(len(b))),1)
        b=np.array(b)
        x=[0]*l_a

        arr=np.insert(arr,l_a+l_b,b,1)
        B=np.array([0]*l_b)
        C=np.array(a+[0]*l_b)
        bx=list(range(l_a,l_a+l_b))
        while(1):
            maxpos=[]
            minpos=[]
            for i in range(l_a+l_b):
                maxpos.append(C[i]-np.dot(B,arr[:,i]))

            c=maxpos.index(max(maxpos))

            if maxpos[c]<=0:
                break
            for i in range(l_b):
                if arr[i,c]!=0:
                    minpos.append(arr[i,-1]/arr[i,c])
                else:
                    minpos.append(10000000000)
            q=filter(lambda x:x>0,minpos)
            if q==[]:
                return "unbounded function"
                break
            r=minpos.index(min(q))

            B[r]=C[c]
            bx[r]=c
            arr[r,:]=arr[r,:]/arr[r,c]
            for i in range(l_b):
                if(i!=r):
                    arr[i,:]=arr[i,:]-(arr[i,c]*arr[r,:])
            if len(a)==2:
                if 0 in bx:
                    m=arr[bx.index(0)][-1]
                else:
                    m=0
                if 1 in bx:
                    n=arr[bx.index(1)][-1]
                else:
                    n=0
                z=int(a[0]*m + a[1]*n)
                plt.text(m+0.25, n+0.25, 'Z='+str(z))
                plt.plot(m, n, 'ko')
                fill_x.append(m)
                fill_y.append(n)

        for i in range(len(bx)):
            if bx[i]<l_a:
                x[bx[i]]=arr[i,-1]
        x=list(map(int,x))
        if len(a)==2:
            plt.plot(0, 0, 'ko', label='Candidate solution after each iteration')
            plt.plot(x[0], x[1], 'rD', label='Optimal Solution')
            plt.fill(fill_x, fill_y, color='aqua')
            plt.legend(loc='upper right', numpoints=1)
            plt.show()
        x.append(sum([x[i]*a[i] for i in range(l_a)]))
        x[-1]=(x[-1],)
        return x

#sol=LPsolver()
#solution=sol.solve([3,2],[18,42,24],[[2,1],[2,3],[3,1]])
#print(solution)
