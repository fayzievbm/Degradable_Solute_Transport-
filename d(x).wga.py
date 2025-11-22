# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 20:45:48 2024

@author: User
"""
# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import matplotlib.pyplot as plt
import copy
import modul2
import math
#to'r parametrlari
nn=41
L=41
tm=3001#3001
tmax=range(1,tm) #vaqt
tau=3
h=0.01


#model parametrlari
tetta=0.3
v=1e-4
c00=0.01
D=5e-6

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
w=0

#Procedure Nach_Dan;{Boshlang'ich va chegaraviy shartlar}
ca=np.zeros((nn,tm), dtype=float)
cp=np.zeros((nn,tm), dtype=float)
c=np.zeros((nn,tm), dtype=float)
alpha=np.zeros(nn, dtype=float)
bet=np.zeros(nn, dtype=float)
F=np.zeros((nn,tm), dtype=float)
c[0,:]=c00
x=np.zeros(nn+1, dtype=float)
for i in range(nn):
    x[i]=i*h
#for j in range(tm):
  # c[0,j]=c00*(1-math.exp(-w*j*tau))
    #c[0,j]=c00#*(1+math.sin(w*j*tau))
   # if abs(c[0,j]-0.02)<0.0001:
   #  print(j)
#print(np.max(c[0,:]))
#print(np.min(c[0,:]))
#  /// bettap ni fisoblash funksiyasi
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

# def progonka():
#     A = tau * v * h + tau * D
#     B = h * h + tau * v * h + 2 * tau * D+lambdae * h * h * tau
#     E = tau * D
#     alpha[1] = 0
#     bet[1] =c00#*(1+math.sin(w*j*tau))
#     for i in range(1, nn - 1):
#         #F[i, j] = (h * h - lambdae*h * h*tau) * c[i, j - 1] - 1 / tetta * (ca[i, j] - ca[i, j - 1]) * h * h - 1 / tetta * (
#                 #cp[i, j] - cp[i, j - 1]) * h * h
#         F[i, j] = (h * h ) * c[i, j - 1] - 1 / tetta * (ca[i, j] - ca[i, j - 1]) * h * h - 1 / tetta * (
#                       cp[i, j] - cp[i, j - 1]) * h * h
#         bet[i + 1] = (F[i, j] + A * bet[i]) / (B - A * alpha[i])
#         alpha[i + 1] = E / (B - A * alpha[i])
#         c[nn-1, j] = 0
#     #c[nn - 2,:]=bet[nn-1]/(1-alpha[nn-1])
#     # (F+A*bet[nn-1])/(B-A*alpha[nn-1])
   
#     for i in reversed(range(nn - 1)):
        
#         c[i, j] = alpha[i + 1] * c[i + 1, j] + bet[i + 1]
#         #c[i,j]=bet[i+1]/(1-alpha[i-1])


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
    # bet[1]=c00#(1-math.exp(-w*j*tau))
     bet[1] = c00
     for i in range(1, nn - 1):
        F[i, j] = (h * h ) * c[i, j - 1] - 1 / tetta * (ca[i, j] - ca[i, j - 1]) *h*h - 1 / tetta * (
                 cp[i, j] - cp[i, j - 1]) *h*h
        bet[i + 1] = (F[i, j] + A * bet[i]) / (B - A * alpha[i])
        alpha[i + 1] = E / (B - A * alpha[i])
        #c[nn-1, j] = 0
        c[L-1, j] = bet[L-1]/(1-alpha[L-1])
        
     for i in reversed(range(nn - 1)):
         
         c[i, j] = alpha[i + 1] * c[i + 1, j] + bet[i + 1]
         
       
     '''
     for i in range(1,nn):
         c[i,j]=(tetta*h*c[i,j-1]+v*tau*c[i-1,j]-h*(ca[i,j]-ca[i,j-1]+cp[i,j]-cp[i,j-1]))/(tetta*h+v*tau)
     '''
#grafiklar

c=c/c00
ca11=copy.copy(ca)
cp11=copy.copy(cp)
c11=copy.copy(c)

#1-qiymatlar yakun


#2-qiymatlar uchun
ca2=np.zeros((nn,tm), dtype=float)
cp2=np.zeros((nn,tm), dtype=float)
c2=np.zeros((nn,tm), dtype=float)

#to'r parametrlari
import numpy as np
import matplotlib.pyplot as plt
#import copy
import modul2

#to'r parametrlari
#nn=61
#L=41
tm=3001#3001
tmax=range(1,tm) #vaqt
tau=3
h=0.01

#model parametrlari
tetta=0.3
v=1e-4
c00=0.01
D=5e-6

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
w=1
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
    #c[0,j]=c00#*(1+math.sin(w*j*tau))
   # if abs(c[0,j]-0.02)<0.0001:
   #  print(j)
#print(np.max(c[0,:]))
#print(np.min(c[0,:]))
#  /// bettap ni fisoblash funksiyasi
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

# def progonka():
#     A = tau * v * h + tau * D
#     B = h * h + tau * v * h + 2 * tau * D+lambdae * h * h * tau
#     E = tau * D
#     alpha[1] = 0
#     bet[1] =c00#*(1+math.sin(w*j*tau))
#     for i in range(1, nn - 1):
#         #F[i, j] = (h * h - lambdae*h * h*tau) * c[i, j - 1] - 1 / tetta * (ca[i, j] - ca[i, j - 1]) * h * h - 1 / tetta * (
#                 #cp[i, j] - cp[i, j - 1]) * h * h
#         F[i, j] = (h * h ) * c[i, j - 1] - 1 / tetta * (ca[i, j] - ca[i, j - 1]) * h * h - 1 / tetta * (
#                       cp[i, j] - cp[i, j - 1]) * h * h
#         bet[i + 1] = (F[i, j] + A * bet[i]) / (B - A * alpha[i])
#         alpha[i + 1] = E / (B - A * alpha[i])
#         c[nn-1, j] = 0
#     #c[nn - 2,:]=bet[nn-1]/(1-alpha[nn-1])
#     # (F+A*bet[nn-1])/(B-A*alpha[nn-1])
   
#     for i in reversed(range(nn - 1)):
        
#         c[i, j] = alpha[i + 1] * c[i + 1, j] + bet[i + 1]
#         #c[i,j]=bet[i+1]/(1-alpha[i-1])

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
    # bet[1]=c00#(1-math.exp(-w*j*tau))
     bet[1] = c00
     for i in range(1, nn - 1):
        F[i, j] = (h * h ) * c[i, j - 1] - 1 / tetta * (ca[i, j] - ca[i, j - 1]) *h*h - 1 / tetta * (
                 cp[i, j] - cp[i, j - 1]) *h*h
        bet[i + 1] = (F[i, j] + A * bet[i]) / (B - A * alpha[i])
        alpha[i + 1] = E / (B - A * alpha[i])
        #c[nn-1, j] =0  
        c[L-1, j] = bet[L-1]/(1-alpha[L-1])
     for i in reversed(range(nn - 1)):
         
         c[i, j] = alpha[i + 1] * c[i + 1, j] + bet[i + 1]
         
       
     '''
     for i in range(1,nn):
         c[i,j]=(tetta*h*c[i,j-1]+v*tau*c[i-1,j]-h*(ca[i,j]-ca[i,j-1]+cp[i,j]-cp[i,j-1]))/(tetta*h+v*tau)
     '''
#grafiklar

c=c/c00

ca2=copy.copy(ca)
cp2=copy.copy(cp)
c2=copy.copy(c)

#2-qiymatlar yakun


#3-qiymatlar uchun
ca3=np.zeros((nn,tm), dtype=float)
cp3=np.zeros((nn,tm), dtype=float)
c3=np.zeros((nn,tm), dtype=float)

#to'r parametrlari
import numpy as np
import matplotlib.pyplot as plt
#import copy
import modul2

#to'r parametrlari
#to'r parametrlari
#nn=51
#L=41
tm=3001#3001
tmax=range(1,tm) #vaqt
tau=3
h=0.01

#model parametrlari
tetta=0.3
v=1e-4
c00=0.01
D=5e-6

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
w=10
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
    #c[0,j]=c00#*(1+math.sin(w*j*tau))
   # if abs(c[0,j]-0.02)<0.0001:
   #  print(j)
#print(np.max(c[0,:]))
#print(np.min(c[0,:]))
#  /// bettap ni fisoblash funksiyasi
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

# def progonka():
#     A = tau * v * h + tau * D
#     B = h * h + tau * v * h + 2 * tau * D+lambdae * h * h * tau
#     E = tau * D
#     alpha[1] = 0
#     bet[1] =c00#*(1+math.sin(w*j*tau))
#     for i in range(1, nn - 1):
#         #F[i, j] = (h * h - lambdae*h * h*tau) * c[i, j - 1] - 1 / tetta * (ca[i, j] - ca[i, j - 1]) * h * h - 1 / tetta * (
#                 #cp[i, j] - cp[i, j - 1]) * h * h
#         F[i, j] = (h * h ) * c[i, j - 1] - 1 / tetta * (ca[i, j] - ca[i, j - 1]) * h * h - 1 / tetta * (
#                       cp[i, j] - cp[i, j - 1]) * h * h
#         bet[i + 1] = (F[i, j] + A * bet[i]) / (B - A * alpha[i])
#         alpha[i + 1] = E / (B - A * alpha[i])
#         c[nn-1, j] = 0
#     #c[nn - 2,:]=bet[nn-1]/(1-alpha[nn-1])
#     # (F+A*bet[nn-1])/(B-A*alpha[nn-1])
   
#     for i in reversed(range(nn - 1)):
        
#         c[i, j] = alpha[i + 1] * c[i + 1, j] + bet[i + 1]
#         #c[i,j]=bet[i+1]/(1-alpha[i-1])


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
    # bet[1]=c00#(1-math.exp(-w*j*tau))
     bet[1] = c00
     for i in range(1, nn - 1):
        F[i, j] = (h * h ) * c[i, j - 1] - 1 / tetta * (ca[i, j] - ca[i, j - 1]) *h*h - 1 / tetta * (
                 cp[i, j] - cp[i, j - 1]) *h*h
        bet[i + 1] = (F[i, j] + A * bet[i]) / (B - A * alpha[i])
        alpha[i + 1] = E / (B - A * alpha[i])
        #c[nn-1, j] =0
        c[L-1, j] = bet[L-1]/(1-alpha[L-1])
        
     for i in reversed(range(nn - 1)):
         
         c[i, j] = alpha[i + 1] * c[i + 1, j] + bet[i + 1]
         
       
     '''
     for i in range(1,nn):
         c[i,j]=(tetta*h*c[i,j-1]+v*tau*c[i-1,j]-h*(ca[i,j]-ca[i,j-1]+cp[i,j]-cp[i,j-1]))/(tetta*h+v*tau)
     '''
#grafiklar

c=c/c00

ca3=copy.copy(ca)
cp3=copy.copy(cp)
c3=copy.copy(c)

#3-qiymatlar yakun
nk=41
x = np.linspace(0, h*(nk-1), num=nk, endpoint=True)
#4-rasm
#modul2.grafik3d(x,c11,1e3*ca11,1e3*cp11,c2,1e3*ca2,1e3*cp2,c3,1e3*ca3,1e3*cp3,tm,tau,nk)
#modul2.grafikdamega(x,c11,1e3*ca11,1e3*cp11,c2,1e3*ca2,1e3*cp2,c3,1e3*ca3,1e3*cp3,tm,tau,nk)

#5-rasm
#4-rasm
#modul2.grafik3d(x,c11,1e3*ca11,1e3*cp11,c2,1e3*ca2,1e3*cp2,c3,1e3*ca3,1e3*cp3,tm,tau,nk)
#5-rasm
#modul2.grafik3ca1(x,c11,1e3*ca11,1e3*cp11,c2,1e3*ca2,1e3*cp2,c3,1e3*ca3,1e3*cp3,tm,tau,nk)
#6-rasm
#modul.grafik3ca0(x,c11,1e3*ca11,1e3*cp11,c2,1e3*ca2,1e3*cp2,c3,1e3*ca3,1e3*cp3,tm,tau,nk)
#7-rasm
#modul2.grafik3v(x,c11,1e3*ca11,1e3*cp11,c2,1e3*ca2,1e3*cp2,c3,1e3*ca3,1e3*cp3,tm,tau,nk)
#modul.grafik3cp0(x,c11,1e3*ca11,1e3*cp11,c2,1e3*ca2,1e3*cp2,c3,1e3*ca3,1e3*cp3,tm,tau,nk)
#8-rasm
#modul.grafik3lep(x,c11,1e3*ca11,1e3*cp11,c2,1e3*ca2,1e3*cp2,c3,1e3*ca3,1e3*cp3,tm,tau,nk)
#modul.grafik3bp0(x,c11,1e5*ca11,1e5*cp11,c2,1e5*ca2,1e5*cp2,c3,1e5*ca3,1e5*cp3,tm,tau,nn)
#9-rasm
tt = np.linspace(0, tau*(tm-1), num=tm, endpoint=True)
N1=1
#modul.grafik3v(x,c11,1e3*ca11,1e3*cp11,c2,1e3*ca2,1e3*cp2,c3,1e3*ca3,1e3*cp3,tm,tau,nk)
modul2.grafik3bp0t(N1,tt,c11,1e3*ca11,1e3*cp11,c2,1e3*ca2,1e3*cp2,c3,1e3*ca3,1e3*cp3,tm)
#np.savetxt('c1.csv',c11)
#np.savetxt('c2.csv',c2)
#np.savetxt('c3.csv',c3)
#np.savetxt('ca1.csv',ca11)
#np.savetxt('ca2.csv',ca2)
#np.savetxt('ca3.csv',ca3)
print('Ok')
#u=np.genfromtxt('24.11.2022 /u1='+str(k)+'.csv')
#modul.grafik3le(x,c11,1e5*ca11,1e5*cp11,c2,1e5*ca2,1e5*cp2,c3,1e5*ca3,1e5*cp3,tm,tau,nn)
