#Analytical solution of the NPKE for insertions.
#Authors: Carlos Cruz (cacl.nucl@gmail.com), Gilberto Espinosa (gpmex
from itertools import combinations
import numpy as np
import math

#***************************************************************************
# Function S_(m,n) 
def Suma(m, L):
   s = 0
   for k in list(combinations(L,m)):
       s = s+np.prod(np.array(k))
   return round(s,8)

#******************************************
# Function S_(i,m,n)
def Suma_i(i,m,L):
   L_i=L[:]                                              #Creates a copy of the list where the decay lambdas are storage.
   L_i.remove(L[i])                                      #Removes the decay lambda in the i-position
   s=0
   for k in list(combinations(L_i,m)):
       s = s+np.prod(np.array(k))
   return round(s,8)

#************************************************
# Function Polyn_coeff_P(L,C_P,C_P_d, rho,Betas,l)
def Polyn_coeff_P(L,C_P, C_P_d, rho, Betas,l):           #The function returns two vectors with the coefficients of C_p and C_p_d
   s_1, s_2, s_3, s_4 = 0, 0, 0, 0                       #Variables used for sums
   bet_tot = sum(Betas)                                  #β
   u = (rho-bet_tot)/l                                   #u = (ρ − β) /l
   C_P.append(1)                                         # Coefficient of x^(K+1)
   C_P.append(Suma(1,L)-u)                               # Coefficient of x^K
   for i in range(len(L)):
       s_1 = s_1 + L[i]*Betas[i]
   C_P.append(Suma(2,L)-u*Suma(1,L)-(1/l)*s_1)           #Coefficient of x^(K−1)
   for i in range(3,len(L)+1):
       s_5 = 0
       for j in range(len(L)):
           s_5 = s_5 + L[j]*Betas[j]*Suma_i(j,i-2,L)
       C_P.append(Suma(i,L)-u*Suma(i-1,L)-(1/l)*s_5)     #Coefficients of x^(K−2),..., x^1
   for i in range(len(L)):
       s_3 = s_3 + L[i]*Betas[i]*Suma_i(i,len(L)-1,L)
   C_P.append(-u*Suma(len(L),L)-(1/l)*s_3)               #Constant coefficient
   for k in range(len(C_P)-1):
       C_P_d.append(C_P[k]*(len(L)+1-k))                 #Coefficients of dP(x)/dx

#****************************************************
# Function Polyn_coeff_H(L, C_H, rho, Betas, l, C_init)
def Polyn_coeff_H(L, C_H, rho, Betas, l, C_init):        #Functions that returns the coefficients of the polynomial H given in Eq. (15)
     bet_tot = sum(Betas)                                #β 
     u = (rho-bet_tot)/l                                 #u = (ρ − β) /l
     s,a = 0, 1
     for i in range(len(C_init)):
         a = L[i]*C_init[i]                              #Product of lambda by the initial conditions of C_i(t)
         s = s + a
     C_H.append(s)                                       # Coefficient of x^k
     for j in range(2,len(L)+1):
         b = 0
         for i in range(len(L)):
             b = b + L[i]*C_init[i]*Suma_i(i,j-1, L)       
         C_H.append(b)

#****************************************************
#Function Polyn_coeff_Q(L, C_q,rho,Betas,l)
def Polyn_coeff_Q(L, C_q,rho,Betas,l):                  #Functions that returns the coefficients of the polynomial H given in Eq. (16)
     bet_tot = np.sum(np.array(Betas))                  #Total Beta
     u = (rho-bet_tot)/l                                #u = (ρ − β) /l
     C_q.append(1)                                      #First coefficient of Q(s) equals to 1
     for j in range(len(L)):
        C_q.append(Suma(j+1,L))

#Insertion solution
def Polynomial_evaluation(Coefficients, value):
    a =0
    for i in range(len(Coefficients)):
        a = a+Coefficients[i]*(value**(len(Coefficients)-1-i))
    return(a)


def n_solution(L,Betas,Lambda_m,n_0,C_init,rho,t):
                                                        #List for the coefficients
    Coefficients_P = [ ]
    Coefficients_dP = [ ]
    Coefficients_H = [ ]
    Coefficients_Q = [ ]

                                                        #Building the coefficients
    Polyn_coeff_P(L,Coefficients_P, Coefficients_dP, rho, Betas,Lambda_m)
    Polyn_coeff_H(L,Coefficients_H, rho, Betas, Lambda_m, C_init)
    Polyn_coeff_Q(L,Coefficients_Q,rho,Betas,Lambda_m)

                                                        #Computing the roots
    roots = list(np.roots(Coefficients_P))              #Computing the roots of P(s)

                                                        #Computing the sum
    a=0
    for i in range(len(roots)):
        derivative = Polynomial_evaluation(Coefficients_dP,roots[i])
        Q_valuated = Polynomial_evaluation(Coefficients_Q,roots[i])
        H_valuated = Polynomial_evaluation(Coefficients_H,roots[i])
        a=a+((n_0*Q_valuated+H_valuated)/derivative)*np.exp(roots[i]*t)
        
    return a



#************************Execution***************************************
#************************************************************************
#************************************************************************

#***********************Input Parameters ********************************
# L = list with lambda constants
# Betas = list with fractions of the precursors
# Lambda_m is the prompt neutron generation time
# Beta_total is computed by default
# Reactivity can be expressed in dollars or numerically

L=[0.0127, 0.0317, 0.115, 0.311,1.40, 3.87]
Betas = [0.000285,0.0015975,0.00141,0.0030525,0.00096,0.000195]
Lambda_m=0.0005
Beta_total = sum(Betas) #This line does not require modifications
rho = -Beta_total       #Reactivity expressed in dollars
time = 10

#********************* Initial conditions ******************************
# n_0= neutron density at t=0
# C_init = list with the initial conditions for the precursors
n_0=1
C_init = [ ]

#************* Initial conditions given by n(0)b_k/(Lambda_m*lambda_k)
for k in range(len(Betas)):
    C_init.append((Betas[k]/(L[k]*Lambda_m))*n_0)

#*************** Calling the analytical solution ************************
result = n_solution(L,Betas,Lambda_m,n_0, C_init,rho,time)
print("Neutron density:",result)
