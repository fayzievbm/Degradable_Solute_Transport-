# -*- coding: utf-8 -*-
"""
Created on Wed May  8 17:44:36 2024

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
import math
#import copy
import modul2
import modulamega

#to'r parametrlari
nn=51
#L=41
tm=3001#3001
tmax=range(1,tm) #vaqt
tau=3
h=0.01

#model parametrlari
tetta=0.3
v=1e-4
c00=0.01
D=0#.5e-6

ca0=4e-3
cp0=8e-3
ca1=2e-4
cp1=3e-3

betp0=3e-4
betar=1.2e-4#/1.2
betaa=1.0e-4#/0.8e-4
betad=1.2e-4#0.5e-4

lambdae=1e-5
lambdaea=2e-5
lambdaep=1.5e-5
w=0.001
#Procedure Nach_Dan;{Boshlang'ich va chegaraviy shartlar}
ca=np.zeros((nn,tm), dtype=float)
cp=np.zeros((nn,tm), dtype=float)
c=np.zeros((nn,tm), dtype=float)
alpha=np.zeros(nn, dtype=float)
bet=np.zeros(nn, dtype=float)
F=np.zeros((nn,tm), dtype=float)
c[0,:]=c00
x=np.zeros(nn, dtype=float)
for i in range(nn):
    x[i]=i*h
#for j in range(tm):
  # c[0,j]=c00*(1-math.exp(-w*j*tau))
   
def bettap(cp,cij):
   if  cp<=cp1:
     return betp0*cij-lambdaep*cp
   elif cp<cp0:
     return betp0*cp1/cp*cij-lambdaep*cp
   else: 
     return 0.0
#  /// bettaa ni fisoblash funksiyasi
def bettaa(ca,cij):
   if  ca<=ca1:
     return betar*cij-lambdaea*ca
   elif ca<ca0:
     return betaa*cij-betad*ca-lambdaea*ca
   else: 
     return 0.0

#Boshlandi
for j in tmax:
     for i in range(nn):
         ca[i,j]=ca[i,j-1]+tau*bettaa(ca[i,j-1],c[i,j-1])
         cp[i,j]=cp[i,j-1]+tau*bettap(cp[i,j-1],c[i,j-1])
         if ca[i,j]>ca0:
             ca[i, j]= ca0
         if cp[i,j]>cp0:
             cp[i, j]= cp0
     #progonka()
         A = (v+w*D*math.exp(-w*x[i]))*tau*h+tau*D*(1+math.exp(-w*x[i]))
         B = h*h+tau*h*(v+w*D*math.exp(-w*x[i]))+2*tau*D*(1+math.exp(-w*x[i]))+tau*h*h*lambdae
         E = tau*D*(1+math.exp(-w*x[i]))
     alpha[1] = 0
    
     bet[1] = c00
     for i in range(1, nn - 1):
        F[i, j] = (h * h ) * c[i, j - 1] - 1 / tetta * (ca[i, j] - ca[i, j - 1]) *h*h - 1 / tetta * (
                 cp[i, j] - cp[i, j - 1]) *h*h
        bet[i + 1] = (F[i, j] + A * bet[i]) / (B - A * alpha[i])
        alpha[i + 1] = E / (B - A * alpha[i])
        c[nn-1,:] = 0
        
     for i in reversed(range(nn - 1)):
         
         c[i, j] = alpha[i + 1] * c[i + 1, j] + bet[i + 1]
         
       
    
#grafiklar

c=c/c00
nk=51
x = np.linspace(0, h*(nk-1), num=nk, endpoint=True)
modulamega.grafikamega(x,c,1e3*ca,1e3*cp,tm,tau,nn)
