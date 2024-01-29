# Analytical Solution of the Neutron Point Kinetic Equations 
The present repository contains the Python codes that were developed to solve the Neutron Point Kinetics Equations, using the Laplace transform and the Heaviside's Theorem. This solution was reported in the paper *A new simplified analytical solution to solve the neutron point kinetics equations using the Laplace transform method*, publised in the the journal of *Computer Physics Communications* in the following [link](https://www.sciencedirect.com/science/article/abs/pii/S0010465522002831?via%3Dihub). 

The programs are licensed under a Creative Commons Attribution 4.0 International License: http://creativecommons.org/licenses/by/4.0/

Authors: Carlos-Antonio Cruz-López (cacl.nucl@gmail.com), Gilberto Espinosa-Paredes (gepe@xanum.uam.mx)

Mathematical and algorithmical generalities of the codes are described in the following lines with the purpose to provide some insight of the developed work. Nevertheless, a more detailed and precise discussion is provided in the submitted article.
## Financial Support.
The authors appreciate the financial support received from the Consejo Nacional de Ciencia y Tecnología, CONAHCYT, under the program “Estancias Posdoctorales por México, 2022”, with the project entitled: “Desarrollo de modelos fenomenológicos energéticos de orden fraccional, para la optimización y simulación en reactores nucleares de potencia”, by which the present development was possible.

## 1. Mathematical description of the problem
The Neutron Point Kinetic Equations (NPKE) with $K$ groups of precursors of delayed neutrons can be written as follows:

$$\frac{dn\left(t\right)}{dt}=\frac{\rho\left(t\right)-\beta}{\Lambda}n\left(t\right)+\sum_{k=1}^{K}{\lambda_kC_k(t)}, \tag{1}$$

$$\frac{dC_k\left(t\right)}{dt}=\frac{\beta_k}{\Lambda}n\left(t\right)-\lambda_kC_k\left(t\right),\ \ 1\le\ k\le\ K, \tag{2}$$
where $n(t)$ denotes the neutron density and $\rho(t)$ is the reactivity. $C(t),\lambda_k,\beta_k$ denote the concentration, the decay constant and the fraction of the $k$-group of precursors of the delayed neutrons, respectively. $\Lambda$ represents the prompt neutron generation time and $\beta$ is given by:

$$\sum_{k=1}^{K}\beta_k=\beta. \tag{3}$$
It is possible to develop an analytical solution for the case of a constant reactivity, i.e. $\rho(t)=\rho_0$ (Lewins, 1978, p. 60), and particularly it is possible to use the Laplace transform for such task. The present codes were developed using this last technique, as well as a particular and non-trivial simplification of the solution. 

## 2. Laplace transform of the system.
After appliying the Laplace transform on both sides of Eq. (1) and Eq. (2), and manipuling the resultant expression, the following equations is obtained for the neutron density:

$$\widetilde{n}\left(s\right)=\frac{n\left(0\right)+\Omega(s)}{s-\frac{\rho-\beta}{\Lambda}-\frac{1}{\Lambda}\Psi(s)} \tag{4}$$

where:
$$\Omega(s)=\sum_{k=1}^{K}\frac{\lambda_kC_k(0)}{s+\lambda_k} \tag{5}$$
and:
$$\Psi(s)=\sum_{k=1}^{K}\frac{\lambda_k\beta_k}{s+\lambda_k} \tag{6}.$$
And for the concentration of the precursors of the delayed neutrons:
$${\widetilde{C}}_k\left(s\right)=\frac{\beta}{\Lambda}\frac{\widetilde{n}\left(s\right)}{s+\lambda_k}+\frac{C\left(0\right)}{s+\lambda_k} \tag{7}$$

One of the main novelties of the analytical solution developed in the mentioned paper consists of writting the last expressions as:
$$\widetilde{n}\left(s\right)=\left(\frac{n\left(0\right)Q\left(s\right)}{P\left(s\right)}+\frac{H\left(s\right)}{P\left(s\right)}\right) \tag{8}$$
where the polynomials $Q(s), H(s)$ and $P(s)$ are given by:
$$Q\left(s\right)=\prod_{k=1}^{K}{(s+\lambda_k)},\ H\left(s\right)=\sum_{k=1}^{K}{C_k(0)}\prod_{j=1,j\neq k}^{K}{(s+\lambda_k}) \tag{9}$$
and:
$$P\left(s\right)=s\prod_{k=1}^{K}{(s+\lambda_k)}-\frac{\rho-\beta}{\Lambda}\prod_{k=1}^{K}\left(s+\lambda_k\right)-\frac{1}{\Lambda}\sum_{k=1}^{K}{\lambda_k\beta_k}\prod_{j=1,j\neq k}^{K}{(s+\lambda_j}). \tag{10}$$

## 3. Analytical solutions.
The Laplace transform system given in Eq. (7) and (8) can be solved using the Heaviside's Theorem Expansion (Arfken et al., 2013, p. 1045), which states that:

>**Theorem 1**:
For two polynomials $L(s)$ and $M(s)$, where the degree of $L(s)$ is less than the degree of $M(s)$, it follows that:
$$\frac{L\left(s\right)}{M(s)}=\sum_{i}\frac{L\left(m_i\right)}{M\prime(m_i)}\frac{1}{s-m_i} \tag{11}$$
where $m_i$ denotes the roots of $M(s)$ and $M'(s)= dM(s)/ds$.
 
Using the last theorem, the analytical solutions of the neutron density is given by:
$$n\left(t\right)=\sum_{k=1}^{K+1}{\frac{n\left(0\right)Q\left(p_k\right)+H(p_k)}{P\prime(p_k)}\exp(p_kt)} \tag{12}$$
where $p_k, 1\leq k \leq K+1$ are the roots of the polynomial $P(s)$ given in Eq. (10), and $P'(s)=dP(s)/ds$. The analytical solution for the precursors of the delayed neutrons is given by:
$$C_k\left(t\right)=\frac{\beta_k}{\Lambda}\sum_{j=1}^{K+1}{\frac{n\left(0\right)Q\left(p_j\right)+H(p_j)}{P\prime(p_j)}\frac{\exp{\left(p_jt\right)}-\exp(-\lambda_kt)}{p_j+\lambda_k}}+C_k\left(0\right)\exp(-\lambda_kt) \tag{13}$$

## 4. Simplification of the polynomials.
One of the most important contribution of our work (Cruz-López et al., 2023) consists of simplifying the Polynomials given in Eq. (9) and Eq. (10), which can be conveniently written as:
$$P\left(s\right)=s^{K+1}+\left(S_{1,K}-u\right)s^K+\left(S_{2,K}-uS_{1,K}-\frac{1}{\Lambda}\sum_{i=1}^{K}{\lambda_i\beta_i}\right)s^{K-1}$$
$$+\sum_{i=3}^{K}{\left(S_{i,K}-uS_{i-1,K}-\frac{1}{\Lambda}\sum_{j=1}^{K}{\lambda_j\beta_jS_{i-2,K-1}^j}\right)s^{K+1-i}}-uS_{K,K}-\frac{1}{\Lambda}\sum_{k=1}^{K}{\lambda_k\beta_kS_{K-1,K-1}^k}; \tag{14}$$
and:
$$H\left(s\right)=\left(\sum_{k=1}^{K}{\lambda_kC_k\left(0\right)}\right)s^{K-1}+\sum_{j=2}^{K}\left(\sum_{i=1}^{K}{\lambda_iC_i\left(0\right)}S_{j-1,K-1}^i\right)s^{K-j} ;\tag{15}$$
and:
$$Q\left(s\right)=s^K+\sum_{j=1}^{K}{S_{j,K}s^{K-j}}, \tag{16}$$
where $u=(\rho-\beta)/\Lambda$ and:
$$S_{m,n}=\sum_{k_1=1}^{n-m+1}\sum_{k_2=k_1+1}^{n-m+2}\cdots\sum_{k_m=k_{m-1}+1}^{n}{\lambda_{k_1}\lambda_{k_2}\cdots\lambda_{k_m}}, \tag{17}$$
$$S_{m,n}^i=\sum_{k_1=1,\ k_1\neq i}^{n-m+1}{\ \sum_{k_2=k_1+1,k_2\neq i}^{n-m+2}\cdots}\sum_{k_m=k_{m-1}+1,\ k_m\neq i}^{n}{\lambda_{k_1}\lambda_{k_2}\cdots\lambda_{k_m}}. \tag{18}$$

## 5. Algorithmical implementation.
## 5.1 Sums
The first step in the algorithmical implementation consists of building a code for the sum $S_{m,n}$ given in Eq. (18). Such expression can be interpreted in combinatorial terms, which allows programming it in a straighforward way. For example for the case of $m=2$, it follows that:
$$S_{2,n}=\sum_{k_1=1}^{n-2+1}\sum_{k_2=k_1+1}^{n-2+2}{\lambda_{k_1}\lambda_{k_2}} \tag{19}$$
which can be expanded as:
$$S_{2,n}=\sum_{k_1=1}^{n-1}\sum_{k_2=k_1+1}^{n}{\lambda_{k_1}\lambda_{k_2}}=\sum_{k_2=2}^{n}{\lambda_1\lambda_{k_2}}+\sum_{k_2=3}^{n}{\lambda_2\lambda_{k_2}}+\ldots+\sum_{k_2=n}^{n}{\lambda_{n-1}\lambda_{k_2}} \tag{20}$$
which can be understood as the sum of all combinations of two elements of the set $\Omega=\set{\lambda_1,\lambda_2,\ldots,\lambda_n}$. Such observation can be extended to the general case, which can be proved using the Vieta's formulas (Vinberg, 2003, p. 89). One of the main computational contributions of our work was provided a simple algorithm for the Eq. (19), which is given in the following code:

### Code 1 ###
```Python
# Function S_(m,n)

from itertools import combinations
import numpy as np
def Suma (m, L):
   s = 0
   for k in list(combinations(L,m)):
       s = s+np.prod(np.array(k))
   return round(s,8)
```
## 5.2 Shifted Sums.
As it can be observed, the sums given in Eq. (18) have an important restriction in each sum, which implies a disadvantages in terms of the execution's time. We developed an elementary way to improve this step, defining a new set of elements where the one that is given in the $i$-position is removed. For example, the element in the third position of the following set will be ommited:
$$\Omega=\set{\lambda_1,\lambda_2,\lambda_3,\lambda_4,\ldots,\lambda_{n-1},\lambda_n} \tag{21}$$
$$\downarrow$$
$$\Omega_3=\set{\lambda_1,\lambda_2,\lambda_4,\ldots,\lambda_{n-1},\lambda_n}. \tag{22}$$
Once that such element was removed, it is possible to use the Eq. (17) instead of the Eq. (18), using a new set of index as follows:
$$\Omega_3=\set{\lambda_1,\lambda_2,\lambda_4,\ldots,\lambda_{n-1},\lambda_n}\rightarrow \Psi=\set{\psi_{1,3},\psi_{2,3},\psi_{3,3},\ldots,\psi_{n-1,3},\psi_{n,3}}. \tag{23}$$
Finally, instead of using the Eq. (18), it is possible to use the Eq. (17) to the set $\Psi$, avoinding in such way the disadvantage restriction that as mentioned before. In the following code such procedure is implemented:
### Code 2 ###
```Python
# Function S_(i,m,n)

from itertools import combinations
import numpy as np
def Suma_(i,m, L):
   L_i=L[:]     #Creates a copy of the list where the decay lambdas are storage.
   L_i.remove(L[i])    #Removes the decay lambda in the i-position
   s=0
   for k in list(combinations(L_i,m)):
       s = s+np.prod(np.array(k))
   return round(s,8)
```
## 5.3 Polynomials 
The implementation of the Polynomials is given in the functions **Polyn_coeff_P**, **Polyn_Coeff_H**, and **Polyn_coeff_Q**, which require the following arguments:
>Polyn_coeff_P $\leftarrow \set{\Omega, \mathcal{P}, \mathcal{P}^\prime,\rho,\mathcal{B},\Lambda}$
>
>Polyn_Coeff_H $\leftarrow\set{\Omega,\mathcal{H},\rho,\mathcal{B},\Lambda,\left[C_1(0),C_2(0),\ldots,C_n(0)\right]}$
>
>Polyn_Coeff_Q $\leftarrow \set{\Omega,\mathcal{Q},\rho,\mathcal{B}}$

where:
> [!IMPORTANT]
> - $\Omega=\set{\lambda_{1}\lambda_{2},\ldots,\lambda{n}},$ is the set that contains all the decay constants of the precursors of the delayed neutrons, which is denoted by **L** in the Pytho's code.
> - $\mathcal{P}$ is a list (or vector) whose elements are the coefficients of the $P(s)$ polynomial given in Eq. (14), which is denoted by the **C_P** variable in the codes.
> - $\mathcal{P}^\prime$ is a list (or vector) that contains the coefficients of the derivative of $P(s)$, given in Eq. (14), that is represented by **C_P_d**.
> - $\rho$ is the reactivity, which is considered as constant and it is represented by **rho**.
> - $\mathcal{B}$ is a list (or vector) whose elements are the fractions $\beta_{k},1\leq k \leq K$, of the precursors of the delayed neutrons. This is denoted by the variable **Betas** in the code.
> - $\Lambda$ has the meaning described in Section 2, and it is given by the **l** variable.
> - $\mathcal{H}$ is a list (or vector) whose elements are the coefficients of the $H(s)$ polynomial given in Eq. (15), and it is denoted by **C_H**.
> - $[C_1(0),C_2(0),\ldots,C_n(0)]$ is a vector whose entries are the initial conditions of the groups of precursors of delayed neutrons. Such vector is denoted in the code as **C_init**.
> - $Q(s)$ is a list (or vector) whose elements are the coefficients of the polynomial given in Eq.(9), which is given by **C_q** in the code. 

### 5.3.1 Codes that compute the polynomials coefficients.
The following three codes computes the coefficients of the polynomials that are used in the analytical solution.
### Code 3: Polyn_coeff_P
As it was mentioned in **Section 5.3**, this code admits as inputs the decay constants and the $\beta$ values (as lists), as well as the reactivity $\rho$ and $\Lambda$. Two lists are generated by the code which contains the coefficients of $P(s)$ and its derivative. This code requieres the definition of **Code 1** and **Code 2**, which appears as **Suma** and **Suma_i**, respectively.
```Python
# Function Polyn_coeff_P(L,C_P,C_P_d, rho,Betas,l)

import numpy as np
def Polyn_coeff_P(L,C_P, C_P_d, rho, Betas,l):           #The function returns two vectors with the coefficients of C_p and C_p_d
   s_1, s_2, s_3, s_4 = 0, 0, 0, 0                       #Variables used for sums
   bet_tot = np.sum(np.array(Betas))                     #β
   u = (rho-bet_tot)/l                                   #u = (ρ − β) /l
   C_P.append(1)                                         # Coefficient of x^(K+1)
   C_P.append(Suma(1,L)-u)                               # Coefficient of x^K
   for i in range(len(L)):
       s_1 = s_1 + L[i] ∗ Betas[i]
   C_P.append(Suma(2,L)-u*Suma(1,L)-(1/l)*s_1)           #Coefficient of x^(K−1)
   for i in range(3,len(L)+1):
       s_5 = 0
       for j in range(len(L)):
           s_5 = s_5 + L[j]*Betas[j]*Suma_i(j,i − 2, L)
       C_P.append(Suma(i,L)-u*Suma(i-1,L)-(1/l)*s_5)     #Coefficients of x^(K−2),..., x^1
   for i in range(len(L)):
       s_3 = s_3 + L[i]*Betas[i]*Suma_i(i,len(L) − 1, L)
   C_P.append(-u*Suma(len(L),L)-(1/l)*s_3)               #Constant coefficient
   for k in range(len(C_P)-1):
       C_P_d.append(C_P[k]*(len(L)+1-k))                 #Coefficients of dP(x)/dx
```
### Code 4: Polyn_coeff_H
This code generates the coefficients of the $H(s)$ polynomial that is given by Eq.(15). It requires the same input of the **Polyn_coeff_P**, but in addition it needs the initial conditions of the precursors of the delayed neutrons given as a vector $[C_1(0),C_2(0),...,C_n(0)]$. 
> [!IMPORTANT]
> -The **Polyn_coeff_P** is the only code related to the polynomials that requieres the use of initial conditions. It will be very important for cases with non-constant reactivities, as it will be discussed in Section 6.2.

```Python

def Polyn_coeff_H(L, C_H, rho, Betas, l, C_init):              #Functions that returns the coefficients of the polynomial H given in Eq. (15)
     bet_tot = np.sum(np.array(Betas))                         #β = n
     u = (rho-bet_tot)/l                                       #u = (ρ − β) /l
     s,a = 0, 1
     for i in range(len(C_init)):
         a = L[i] ∗ C_init[i]                                  #Product of lambda by the initial conditions of C_i(t)
         s = s + a
     C_H.append(s)                                             # Coefficient of x^k
     for j in range(2,len(L)+1):
         b = 0
         for i in range(len(L)):
             b = b + L[i] ∗ C_init[i] ∗ Suma_i(i,j − 1, L)       
     C_H.append(b)
```
### Code 5: Polyn_coeff_Q
This code provides the coefficients of the $Q(s)$ polynomial given in Eq. (16). It requires the same input nuclear parameters as the **Polyn_coeff_P**, as well as the **Code 1**. 
```Python

def Polyn_coeff_Q(L, C_q,rho,Betas,l):              #Functions that returns the coefficients of the polynomial H given in Eq. (16)
     bet_tot = np.sum(np.array(Betas))              #Total Beta
     u = (rho-bet_tot)/l                            #u = (ρ − β) /l
     C_q.append(1)                                  #First coefficient of Q(s) equals to 1
     for j in range(len(L)):
        C_q.append(Suma(j+1,L))
```

## 5.4 Flow diagram of the Python code's
A flow diagram of the dependence of the codes is provided in the following image. It is worth mentioning that this scheme does not show the discretization method described in Section **6.2**
<details><summary>CLICK HERE to expand the diagram.</summary>
<p>

![image](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE/blob/main/Imagen%201.jpg)


</p>
</details>

## 6. AnalyticNPKE-Insertion.py
The code AnalyticNPKE-Ramp solves the system given in Eq. (13) and Eq. (14). It is provided in the text files of the present repository.
### 6.1 An application. 
The AnalyticNPKE-Insertion code will be used to reproduce the data reported by Nahla (2010, p. 1626). For such scenario the input is the following:
> [!IMPORTANT]
>  Nuclear parameter | Value     | Nuclear parameter | Value           |
| ------------- | ------------- | -------------     | --------------  |
| $\lambda_1$   | 0.000125591s  | lambda            | 0.0769478s^(-1) |
| PNL           | 0.975         | beta              |0.00645          |
| rho           |0.002          |Lambda             |0.00005s         |









