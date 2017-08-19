#!/usr/bin/env python3
# ----------------------------------------- #
#
import numpy as np
from scipy.linalg import solve
from scipy.sparse import diags
from scipy.optimize import bisect

L = 16
beta = bisect(f=(lambda x: np.cosh(x)*np.cos(x)+1), a=5, b=10)
k = beta/L
print(k,k*L,np.cosh(k*L)*np.cos(k*L))

R = (np.cos(k*L)+np.cosh(k*L))/(np.sin(k*L)+np.sinh(k*L))
print(R)

def y(x):
 return 0.5*((np.cosh(k*x)-np.cos(k*x))-R*(np.sinh(k*x)-np.sin(k*x)))

print(y(L))

n = L
x = np.linspace(0.5,n-0.5,n)
A4 = np.array([np.ones(n-2),-4*np.ones(n-1),6*np.ones(n),-4*np.ones(n-1),np.ones(n-2),])
offset = [-2,-1,0,1,2]
A4 = diags(A4,offset).toarray()
A4[0,[0,1]] = [25,-50/9]
A4[1,[0,1]] = [-2,53/9]
# A4[-1,[-1,-2]] = [1,-2]
# A4[-2,[-1,-2]] = [-2,5]
A4[-1,[-1,-2,-3]] = [24/29,-2+10/29,24/29]
A4[-2,[-1,-2,-3]] = [-2+5/29,4+19/29,-4+5/29]

def percent_err(wrong,right):
    print(wrong/right-1)
def err(wrong,right):
    print(wrong-right)

print('---')
b = np.matmul(A4,y(x))
for i in range(n):
    print(i,x[i],y(x[i]),b[i]/k**4,b[i]/k**4-y(x[i]))
