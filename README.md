# Analytical Solution of the Neutron Point Kinetic Equations 
The present repository contains the Python codes that were developed to solve the Neutron Point Kinetics Equations (NPKE), using the Laplace transform and the Heaviside's Theorem. This solution was reported in the paper *A New Simplified Analytical Solution to Solve the Neutron Point Kinetics Equations Using the Laplace Transform Method*, publised in the journal *Computer Physics Communications*. That work can be found in the following [link](https://www.sciencedirect.com/science/article/abs/pii/S0010465522002831?via%3Dihub). 

The programs are licensed under a Creative Commons Attribution 4.0 International License: http://creativecommons.org/licenses/by/4.0/

Authors: Carlos-Antonio Cruz-López (cacl.nucl@gmail.com), Gilberto Espinosa-Paredes (gepe@xanum.uam.mx)

Mathematical and algorithmical generalities of the codes are described in the following lines with the purpose to provide some insight of the developed work. Nevertheless, a more detailed discussion is provided in the submitted article.
## Software specifications.
The AnalyticNPKE codes were written in the Python programming language in its version 3. The codes require the following libraries:

- [x] itertools
- [x] numpy
- [x] math

The itertools and the math libraries are usually included in the installation of the Python 3. Nevertheless, the numpy library needs to be installed because it is not a native library. We suggest tu use the Anaconda distribuitor in order to carry out the installation in a straighforward way. Such distribuitor can be consulted [here](https://www.anaconda.com/download).

## Financial Support.
The authors appreciate the financial support received from the Consejo Nacional de Humanidades, Ciencia y Tecnología (CONAHCYT), under the program *Estancias Posdoctorales por México, 2022*, with the project entitled: *Desarrollo de modelos fenomenológicos energéticos de orden fraccional, para la optimización y simulación en reactores nucleares de potencia*, by which the present development was possible.
## Index of the Repository
1. [Mathematical description of the problem](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#1-mathematical-description-of-the-problem).
1. [Laplace transform of the system](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#2-laplace-transform-of-the-system).
1. [Analytical Solutions](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE/blob/main/README.md#3-analytical-solutions).
1. [Simplification of the Polynomials](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE?tab=readme-ov-file#4-simplification-of-the-polynomials).
2. [Algorithmical Implementation](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE?tab=readme-ov-file#5-algorithmical-implementation).
   - [5.1 Sums](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE?tab=readme-ov-file#51-sums)
   - [5.2 Shifted sums](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE/blob/main/README.md#52-shifted-sums)
   - [5.3 Polynomials](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE/blob/main/README.md#53-polynomials)
     - [5.3.1 Codes that compute the polynomials coefficients.](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE/blob/main/README.md#531-codes-that-compute-the-polynomials-coefficients)
   - [5.4 Evaluation of Polynomials](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE/blob/main/README.md#54-evaluation-of-the-polynomials)
   - [5.5 Flow diagram of the codes](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#55-flow-diagram-of-the-python-codes)
1. [AnalyticNPKE-Insertion.py](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#6-analyticnpke-insertionpy)
   - [6.1 Example of an application of the AnalayticNPKE-Insertion.py](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#61-application)
1. [AnalyticNPKE-Ramp.py](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#7-analyticnpke-ramppy)
   - [7.1 Example of an application of the AnalyticNPKE-Ramp.py](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#71-application)
1. [AnalyticNPKE-Feedback.py](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#8-analyticnpke-feedbackpy)
   - [8.1 Example of an Application of the AnalyticNPKE-Feedback.py](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#81-example-of-an-application-of-the-analyticnpke-feedbackpy)
1. [References](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE?tab=readme-ov-file#9-references)

   
## 1. Mathematical description of the problem
The Neutron Point Kinetic Equations with $K$ groups of precursors of delayed neutrons can be written as follows (Duderstadt and Hamilton, 1976, p. 238):

$$\frac{dn\left(t\right)}{dt}=\frac{\rho\left(t\right)-\beta}{\Lambda}n\left(t\right)+\sum_{k=1}^{K}{\lambda_kC_k(t)}, \tag{1}$$

$$\frac{dC_k\left(t\right)}{dt}=\frac{\beta_k}{\Lambda}n\left(t\right)-\lambda_kC_k\left(t\right),\ \ 1\le\ k\le\ K, \tag{2}$$
where $n(t)$ denotes the neutron density and $\rho(t)$  the reactivity. $C(t),\lambda_k,\beta_k$ are the concentration, the decay constant and the fraction of the $k$-group of precursors of the delayed neutrons, respectively. $\Lambda$ represents the prompt neutron generation time and $\beta$ is given by:

$$\sum_{k=1}^{K}\beta_k=\beta. \tag{3}$$
It is possible to develop an analytical solution for the case of a constant reactivity, i.e. $\rho(t)=\rho_0$ (Lewins, 1978, p. 60), and particularly it is possible to use the Laplace transform for such task. The present codes were developed using this last technique, as well as a particular and non-trivial simplifications of the involved Polynomials. 

## 2. Laplace transform of the system.
After appliying the Laplace transform on both sides of Eq. (1) and Eq. (2), and manipuling the resultant expression, the following equations is obtained for the neutron density:

$$\widetilde{n}\left(s\right)=\frac{n\left(0\right)+\Omega(s)}{s-\frac{\rho-\beta}{\Lambda}-\frac{1}{\Lambda}\Psi(s)} \tag{4}$$

where:
$$\Omega(s)=\sum_{k=1}^{K}\frac{\lambda_kC_k(0)}{s+\lambda_k} \tag{5}$$
and:
$$\Psi(s)=\sum_{k=1}^{K}\frac{\lambda_k\beta_k}{s+\lambda_k} \tag{6}.$$
And for the concentration of the precursors of the delayed neutrons:
$${\widetilde{C}}_k\left(s\right)=\frac{\beta}{\Lambda}\frac{\widetilde{n}\left(s\right)}{s+\lambda_k}+\frac{C\left(0\right)}{s+\lambda_k}. \tag{7}$$

One of the main novelties of the analytical solution developed in the mentioned paper consists of writting the last expressions as:
$$\widetilde{n}\left(s\right)=\left(\frac{n\left(0\right)Q\left(s\right)}{P\left(s\right)}+\frac{H\left(s\right)}{P\left(s\right)}\right) \tag{8}$$
where the polynomials $Q(s), H(s)$ and $P(s)$ are given by:
$$Q\left(s\right)=\prod_{k=1}^{K}{(s+\lambda_k)},\ H\left(s\right)=\sum_{k=1}^{K}{C_k(0)}\prod_{j=1,j\neq k}^{K}{(s+\lambda_k}) \tag{9}$$
and:
$$P\left(s\right)=s\prod_{k=1}^{K}{(s+\lambda_k)}-\frac{\rho-\beta}{\Lambda}\prod_{k=1}^{K}\left(s+\lambda_k\right)-\frac{1}{\Lambda}\sum_{k=1}^{K}{\lambda_k\beta_k}\prod_{j=1,j\neq k}^{K}{(s+\lambda_j}). \tag{10}$$

## 3. Analytical solutions.
The Laplace transform system given in Eq.(8) can be solved using the Heaviside's Theorem Expansion (Arfken et al., 2013, p. 1045), which states that:

>**Theorem 1**:
For two polynomials $L(s)$ and $M(s)$, where the degree of $L(s)$ is less than the degree of $M(s)$, it follows that:
$$\frac{L\left(s\right)}{M(s)}=\sum_{i}\frac{L\left(m_i\right)}{M\prime(m_i)}\frac{1}{s-m_i} \tag{11}$$
where $m_i$ denotes the roots of $M(s)$ and $M'(s)= dM(s)/ds$.
 
Using the last theorem, the analytical solutions of the neutron density is given by:
$$n\left(t\right)=\sum_{k=1}^{K+1}{\frac{n\left(0\right)Q\left(p_k\right)+H(p_k)}{P\prime(p_k)}\exp(p_kt)} \tag{12}$$
where $p_k, 1\leq k \leq K+1$ are the roots of the polynomial $P(s)$ given in Eq. (10), and $P'(s)=dP(s)/ds$. The analytical solution for the precursors of the delayed neutrons can be found, on the other hand, using the Convolution's theorem, after which the following expression is obtained:
$$C_k\left(t\right)=\frac{\beta_k}{\Lambda}\sum_{j=1}^{K+1}{\frac{n\left(0\right)Q\left(p_j\right)+H(p_j)}{P\prime(p_j)}\frac{\exp{\left(p_jt\right)}-\exp(-\lambda_kt)}{p_j+\lambda_k}}+C_k\left(0\right)\exp(-\lambda_kt). \tag{13}$$

As it can be observed, the solutions are expressed as a linear combination of the exponential functions.
## 4. Simplification of the polynomials.
One of the most important contribution of our work (Cruz-López et al., 2023) consists of simplifying the Polynomials given in Eq. (9) and Eq. (10), using Theory of Equations instead of a Matrix Approach, as other authors did. After a new procedure using the Viéte's formula (Vinberg, 2003), it is possible to rewrite the polynomials as:

$$P\left(s\right)=s^{K+1}+\left(S_{1,K}-u\right)s^K+\left(S_{2,K}-uS_{1,K}-\frac{1}{\Lambda}\sum_{i=1}^{K}{\lambda_i\beta_i}\right)s^{K-1}$$

$$+\sum_{i=3}^{K}{\left(S_{i,K}-uS_{i-1,K}-\frac{1}{\Lambda}\sum_{j=1}^{K}{\lambda_j\beta_jS_{i-2,K-1}^j}\right)s^{K+1-i}}-uS_{K,K}-\frac{1}{\Lambda}\sum_{k=1}^{K}{\lambda_k\beta_kS_{K-1,K-1}^k}; \tag{14}$$

and:

$$H\left(s\right)=\left(\sum_{k=1}^{K}{\lambda_kC_k\left(0\right)}\right)s^{K-1}+\sum_{j=2}^{K}\left(\sum_{i=1}^{K}{\lambda_iC_i\left(0\right)}S_{j-1,K-1}^i\right)s^{K-j} ;\tag{15}$$

as well as:

$$Q\left(s\right)=s^K+\sum_{j=1}^{K}{S_{j,K}s^{K-j}}, \tag{16}$$

where $u=(\rho-\beta)/\Lambda$, and:

$$ S_{m,n}=\sum_{k_1=1}^{n-m+1}\sum_{k_2=k_1+1}^{n-m+2} \cdots \sum_{k_m=k_{m-1}+1}^{n} {\lambda_{k_1}\lambda_{k_2}\cdots\lambda_{k_m}}, \tag{17}$$

$$S_{m,n}^i=\sum_{k_1=1,\ k_1\neq i}^{n-m+1}{\ \sum_{k_2=k_1+1,k_2\neq i}^{n-m+2}\cdots}\sum_{k_m=k_{m-1}+1,\ k_m\neq i}^{n}{\lambda_{k_1}\lambda_{k_2}\cdots\lambda_{k_m}}. \tag{18}$$

## 5. Algorithmical implementation.
### 5.1 Sums
The first step in the algorithmical implementation consists of building a code for the sum $S_{m,n}$ given in Eq. (18). Such expression can be interpreted in combinatorial terms, which allows programming it in a straighforward way. For example, for the case of $m=2$, it follows that:
$$S_{2,n}=\sum_{k_1=1}^{n-2+1}\sum_{k_2=k_1+1}^{n-2+2}{\lambda_{k_1}\lambda_{k_2}} \tag{19}$$
which can be expanded as:
$$S_{2,n}=\sum_{k_1=1}^{n-1}\sum_{k_2=k_1+1}^{n}{\lambda_{k_1}\lambda_{k_2}}=\sum_{k_2=2}^{n}{\lambda_1\lambda_{k_2}}+\sum_{k_2=3}^{n}{\lambda_2\lambda_{k_2}}+\ldots+\sum_{k_2=n}^{n}{\lambda_{n-1}\lambda_{k_2}} $$

$$=\lambda_1\lambda_2+\lambda_1\lambda_3+\ldots+\lambda_1\lambda_2+\lambda_2\lambda_3+\lambda_2\lambda_4+\ldots+\lambda_2\lambda_n +\ldots+\lambda_{n-1}\lambda_n, \tag{20}$$

which can be understood as the sum of all combinations of two elements of the set $\Omega=\set{\lambda_1,\lambda_2,\ldots,\lambda_n}$. Such important observation, that was part of the most relevants findings in our work, can be extended to the general case, which can be proved using the Vieta's formulas (Vinberg, 2003, p. 89). The corresponding algorithm for the Eq. (19) is given in the following code:

### Code 1 ###
```Python
# Function S_(m,n)

from itertools import combinations
import numpy as np
def Suma (m, L):
   s = 0
   for k in list(combinations(L,m)):             #Builds the combinations of the set of the Lambda constants
       s = s+np.prod(np.array(k))                #Carries out the product of the combinations. 
   return round(s,8)
```
### 5.2 Shifted Sums.
As it can be observed, the sums given in Eq. (18) have an important restriction in each sum, which implies a disadvantages in terms of the execution's time. We developed an elementary way to improve this step, defining a new set of elements where the one that is given in the $i$-position is removed. For example, the element in the third position of the following set will be ommited:
$$\Omega=\set{\lambda_1,\lambda_2,\lambda_3,\lambda_4,\ldots,\lambda_{n-1},\lambda_n} \tag{21}$$
$$\downarrow$$
$$\Omega_3=\set{\lambda_1,\lambda_2,\lambda_4,\ldots,\lambda_{n-1},\lambda_n}. \tag{22}$$

Once that such element was removed, it is possible to define a new set of index as follows:

$$\Omega_3=\set{\lambda_1,\lambda_2,\lambda_4,\ldots,\lambda_{n-1},\lambda_n}\rightarrow \Psi=\set{\psi_{1,3},\psi_{2,3},\psi_{3,3},\ldots,\psi_{n-2,3},\psi_{n-1,3}}. \tag{23}$$

Finally, instead of using the Eq. (18), it is possible to apply the Eq. (17) to the set $\Psi$, because the element has been removed and there is not need to check, in each sum, that the index is different from it. We avoid, in this way, the disadvantage restriction that as mentioned before. This procedure is implemented in the following code:
### Code 2 ###
```Python
# Function S_(i,m,n)

from itertools import combinations
import numpy as np
def Suma_(i,m, L):
   L_i=L[:]                                 #Creates a copy of the list where the decay lambdas are storage.
   L_i.remove(L[i])                         #Removes the decay lambda in the i-position
   s=0
   for k in list(combinations(L_i,m)):      #Builds the combinations of the set of the Lambda constants
       s = s+np.prod(np.array(k))           #Carries out the product of the combinations. 
   return round(s,8)
```
### 5.3 Polynomials 
The implementation of the Polynomials is given in the functions **Polyn_coeff_P**, **Polyn_Coeff_H**, and **Polyn_coeff_Q**, which require the following arguments:
>Polyn_coeff_P $\leftarrow \set{\Omega, \mathcal{P}, \mathcal{P}^\prime,\rho,\mathcal{B},\Lambda}$
>
>Polyn_Coeff_H $\leftarrow\set{\Omega,\mathcal{H},\rho,\mathcal{B},\Lambda,\left[C_1(0),C_2(0),\ldots,C_n(0)\right]}$
>
>Polyn_Coeff_Q $\leftarrow \set{\Omega,\mathcal{Q},\rho,\mathcal{B}}$

where:
> [!IMPORTANT]
> - $\Omega=\set{\lambda_{1}\lambda_{2},\ldots,\lambda{n}},$ is the set that contains all the decay constants of the precursors of the delayed neutrons, which is denoted by **L** in the Python's code.
> - $\mathcal{P}$ is a list (or vector) whose elements are the coefficients of the $P(s)$ polynomial given in Eq. (14), which is denoted by the **C_P** variable in the codes.
> - $\mathcal{P}^\prime$ is a list (or vector) that contains the coefficients of the derivative of $P(s)$, given in Eq. (14), that is represented by **C_P_d**.
> - $\rho$ is the reactivity, which is considered as constant and it is represented by **rho**.
> - $\mathcal{B}$ is a list (or vector) whose elements are the fractions $\beta_{k},1\leq k \leq K$, of the precursors of the delayed neutrons. This is denoted by the variable **Betas** in the code.
> - $\Lambda$ has the meaning described in Section 2, and it is given by the **l** variable.
> - $\mathcal{H}$ is a list (or vector) whose elements are the coefficients of the $H(s)$ polynomial given in Eq. (15), and it is denoted by **C_H**.
> - $[C_1(0),C_2(0),\ldots,C_n(0)]$ is a vector whose entries are the initial conditions of the groups of precursors of delayed neutrons. Such vector is denoted in the code as **C_init**.
> - $Q(s)$ is a list (or vector) whose elements are the coefficients of the polynomial given in Eq. (9), which is given by **C_q** in the code. 

### 5.3.1 Codes that compute the polynomials coefficients.
The following three codes computes the coefficients of the polynomials that are used in the analytical solution.
### Code 3: Polyn_coeff_P
As it was mentioned in **Section 5.3**, this code admits as inputs the decay constants and the $\beta$ values (as lists), as well as the reactivity $\rho$ and the parameter $\Lambda$. Two lists are generated by the code which contains the coefficients of $P(s)$ and its derivative. This code requieres the definition of **Code 1** and **Code 2**, which appears as **Suma** and **Suma_i**, respectively.
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
This code generates the coefficients of the $H(s)$ polynomial that is given by Eq. (15). It requires the same input of the **Polyn_coeff_P**, but in addition it needs the initial conditions of the precursors of the delayed neutrons given as a vector $[C_1(0),C_2(0),...,C_n(0)]$. 
> [!IMPORTANT]
> -The **Polyn_coeff_H** is the only code related to the polynomials that requieres the use of initial conditions. It will be very important for cases with non-constant reactivities, as it will be discussed in Section 6.2.

```Python

def Polyn_coeff_H(L, C_H, rho, Betas, l, C_init):              #Functions that returns the coefficients of the polynomial H given in Eq. (15)
     bet_tot = np.sum(np.array(Betas))                         #β 
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
     bet_tot = np.sum(np.array(Betas))              #β
     u = (rho-bet_tot)/l                            #u = (ρ − β) /l
     C_q.append(1)                                  #First coefficient of Q(s) equals to 1
     for j in range(len(L)):
        C_q.append(Suma(j+1,L))
```
### 5.4 Evaluation of the Polynomials.
Since the analytical solution requires the evaluation of the polynomials, it is necessary to introduce a brief routine that carries-out this procedure. This routine requires the coefficients of the polynomial as well as the real number where they will be evaluated. The following code contains this function:
### Code 6. Polynomial_evaluation.
```Python

def Polynomial_evaluation(Coefficients, value):   #Function that evaluates the polynomials.
    a =0
    for i in range(len(Coefficients)):
        a = a+Coefficients[i]*(value**(len(Coefficients)-1-i))
    return(a)
```

### 5.5 Flow diagram of the Python code's
A flow diagram of the dependence of the codes is provided in the following image. It is worth mentioning that this scheme does not show the discretization method described in Section **6.2**
<details><summary>CLICK HERE to expand the diagram.</summary>
<p>

![image](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE/blob/main/Imagen%201.jpg)


</p>
</details>

## 6. AnalyticNPKE-Insertion.py
The code **AnalyticNPKE-Insertion** solves the system given in Eq. (13) and Eq. (14) including all the past Codes, as well as a routine that integrates them. It is provided in the text files of the present repository.

>[!WARNING]
> The AnalyticNPKE-Insertion.py code only can be used for cases with constant reactivities. For linear-time reactivities see the AnalyticNPKE-Ramp.py code.

### 6.1 Example of an application of the AnalayticNPKE-Insertion.py
In the present example, the **AnalyticNPKE-Insertion** code will be used to reproduce part of the data reported by Nahla (2010, p. 1626). For such scenario the input parameters are the following: 
|Nuclear parameter | Value  ($\mathrm{s^{-1}}$)| Nuclear parameter | Value           |
| ------------- | ------------- | -------------     | --------------  |
| $\lambda_1$   |0.0127         | $\beta_1$         |0.000285         |
| $\lambda_2$   |0.0317         | $\beta_2$         |0.0015975        |
| $\lambda_3$   |0.115          | $\beta_3$         |0.00141          |
| $\lambda_4$   |0.311          | $\beta_4$         |0.0030525        |
| $\lambda_5$   |1.40           | $\beta_5$         |0.00096          |
| $\lambda_6$   |3.87           | $\beta_6$         |0.000195         |

and $\beta=0.0075$ and $\Lambda=0.0005 \mathrm{s}$. A negative reactivity given by $\rho=-1$ dollar will be used as well as a time of $t=10$ seconds. The following part of the code is considered as the "Input" section. The last nuclear parameters are introduced in such part, as well as the initial conditions given by:

$$n(0)=1, \ \ C_k(0)=\frac{\beta_k n(0)}{\Lambda \cdot \lambda_k}, 1\leq k \leq K.$$

### Input:

<details><summary>CLICK HERE to expand the input of the application of the AnalyticNPKE-Insertion.py</summary>
<p>


```Python

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
```
</p>
</details>

### Output
```Python
Neutron density: 0.23611064482555608
```
which coincides in the first seven digits with the data reported by Nahla (2010, p. 8). 

> [!IMPORTANT]
> The initial conditions for the concentration of the precursors can be modified by hand, instead of using a loop. It is only necessary to modify the following lines in the code:

```Python
C_init = [ ]   #In this part it is necessary to introduce the particular initial conditions 
```
as well as ignore the following lines:

```Python
#************* Initial conditions given by n(0)b_k/(Lambda_m*lambda_k)
for k in range(len(Betas)):
    C_init.append((Betas[k]/(L[k]*Lambda_m))*n_0)
```

## 7. AnalyticNPKE-Ramp.py

Even when the analytical solution was developed assuming a constant reactivity, it also can be used when the reactivity is a linear function of time, i.e., $\rho(t)=\gamma t$, where $\gamma$ is a constant. Such scenario is known as the "ramp reactivity case". To solve it, it is necessary to discretize the time in small intervals, assuming a constant value of the reactivity in each of them, given by:

$$\bar{\rho}=\frac{\ \rho\left(t_n\right)+\rho(t_{n-1})}{2} \tag{24}$$

where the lower and upper times are defined as:

$$t_{n}=\Delta t \cdot n = h \cdot n \tag{25}$$

$$t_{n-1}=\Delta t \cdot (n-1) = h \cdot (n-1). \tag{26}$$
Aditionally, it is necessary to update the initial conditions at the end of each step, considering the following relationships:
$$n_{t_{n-1}}(t_n)=n_{t_n}(0),\ \ \ \ C_{k_{t_{n-1}}}(t_n)=C_{k_{t_n}}(0), \tag{27}$$

where $h$ is the time step.  In other words, the value of the variables at the end of an interval are the initial conditions for the next one. 
The code **AnalyticNPKE-Ramp.py**, that is provided in the repository, includes the last methodology. This code is similar to **AnalyticNPKE-Insertion.py**, considering the following modifications:

> + Since the polynomials depends on the reactivity, they must to be updated in each step.
> + It is necessary to solve the equations related to the precursors because the vector $C=[C_1(0),C_2(0),...,C_m(0)]$ is required in each time step.
> + The solution is provided as a vector who contains, not only the Target time, but also the solution for each time step.
> + 
### 7.1 Example of an application of the AnalayticNPKE-Ramp.py 
The **AnalyticNPKE-Ramp.py** will be used to reproduce data reported by Nahla (2010, p. 9), For such scenario the input parameters are the following: 
|Nuclear parameter | Value  ($\mathrm{s^{-1}}$)| Nuclear parameter | Value           |
| ------------- | ------------- | -------------     | --------------  |
| $\lambda_1$   |0.0127         | $\beta_1$         |0.000266         |
| $\lambda_2$   |0.0317         | $\beta_2$         |0.001491         |
| $\lambda_3$   |0.115          | $\beta_3$         |0.001316         |
| $\lambda_4$   |0.311          | $\beta_4$         |0.002849         |
| $\lambda_5$   |1.40           | $\beta_5$         |0.000896         |
| $\lambda_6$   |3.87           | $\beta_6$         |0.000182         |

and $\beta=0.007$ and $\Lambda= 0.00002 \mathrm{s}$. A negative reactivity given by $\rho$=0.1 dollar will be used as well as a time of $t=2$ seconds.
>[!WARNING]
> We concluded that an adequate time step is given by $h=0.001 \mathrm{s}$, because with this value several reference values were reproduced. Even more,  the execution's time is acceptable, carrying out 2000 iterations in nearly 5 seconds on a 3.80 GHz Desktop Computer, under a Windows 11 operative system. Greater time steps can lead to errors.

The corresponding input and outputs are given in the following sections:

### Input
<details><summary>CLICK HERE to expand the input of the application of AnalyticNPKE-Ramp.py</summary>
<p>


```Python
#***********************Input Parameters ********************************
#************************************************************************
# L = list with lambda constants
# Betas = list with fractions of the precursors
# Lambda_m is the prompt neutron generation time
# Beta_total is computed by default
# Reactivity can be expressed in dollars or numerically
#Gamma slope of the linear reactivity expressed in dollars

L=[0.0127, 0.0317, 0.115, 0.311,1.40, 3.87]
Betas = [0.000266,0.001491,0.001316,0.002849,0.000896,0.000182]
Lambda_m=0.00002
gamma = 0.1                                            

#********************* Initial conditions *******************************
# n_0= neutron density at t=0
# C_init = list with the initial conditions for the precursors
n_0=1

#************* Initial conditions given by n(0)b_k/(Lambda_m*lambda_k)
for k in range(len(Betas)):
    C_init.append((Betas[k]/(L[k]*Lambda_m))*n_0)

#*************************************************************************
#********************* Discretization Input ******************************
#*************************************************************************
    
Target = 2                                            # Target is the desired time
step = 0.001                                          # Recommended value

#************************************************************************
```
</p>
</details>

### Outputs
Two outputs are generated. The first one is provided by the Python's shell, where the time and the neutron density are printed. The second output is given as a txt file whose name is "solutions_final", and which is generated in the same folder where the code is saved. 

#### Output 1. Python's shell.
<details><summary>CLICK HERE to expand the Output of the application of AnalyticNPKE-Ramp.py</summary>
<p>
 
```Python
Only the last 30 values are showed by brevity. But the Python's shell prints all the results.
...
1.97 1.330516908276276
1.971 1.3307711779807945
1.972 1.3310255646844884
1.973 1.331280068451131
1.974 1.3315346893445346
1.975 1.3317894274285604
1.976 1.3320442827671142
1.977 1.3322992554241664
1.978 1.3325543454637243
1.979 1.3328095529498534
1.98 1.3330648779466536
1.981 1.3333203205182909
1.982 1.3335758807289764
1.983 1.3338315586429637
1.984 1.3340873543245553
1.985 1.334343267838108
1.986 1.3345992992480376
1.987 1.334855448618793
1.988 1.3351117160148882
1.989 1.335368101500862
1.99 1.3356246051413316
1.991 1.3358812270009428
1.992 1.336137967144411
1.993 1.336394825636477
1.994 1.3366518025419405
1.995 1.3369088979256796
1.996 1.3371661118525853
1.997 1.3374234443876065
1.998 1.3376808955957533
1.999 1.3379384655420707
2.0 1.3381961542916707
```
</p>
</details>

#### Output 2. A ".txt" file.

Unlike the Python's shell that only contains the results for the neutron density, the "Solutions_final.txt" contains the solution of the precursors of the delayed neutrons. An example of such file can be reviewed in the following link: [Example of an output file](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE/blob/main/Example_of_output_ramp.txt)
## 8. AnalyticNPKE-Feedback.py
The analytic solution can also be applied to non-linear cases where the reactivity depends on the neutro nensity $n(t)$. Particularly in the cases where such relationship is given by:
$$\frac{d\rho\left(t\right)}{dt}=a-bn(t), \tag{27}$$
where $at$ and $b$ represents the impressed reactivity and the shutdown coefficient, respectively (Nahla, 2010, p. 1626). The procedure is simialr to the one described in Section 7, i.e., it is necessary to discretize the time in small intervals, assuming that reactivity is constant in each of them, which is computed approximating Eq. (27) as follows:
$$\frac{d\rho\left(t\right)}{dt}=\approx\frac{\rho\left(t+h\right)-\rho\left(t\right)}{h}\approx\ a-bn(t) \tag{28}$$
where $h$ denotes the time step using for the approximation. Considering the case of $t=0$, Eq. (28) is reduced to:
$$\rho(h)=\rho(0)+(a-bn(0)). \tag{29}$$
**AnalyticNPKE-Feedback.py**, that is provided in this repository, includes this last methodology. Essentially it has the same structure as the **AnalyticNPKE-Ramp.py**, but it updates the reactivity function in a different way. 
### 8.1 Example of an Application of the AnalyticNPKE-Feedback.py
As an example of application, the code will be used to reproduce data that is reported by Nahla (2010, p. 1626). For such task the following data will be used:
|Nuclear parameter | Value  ($\mathrm{s^{-1}}$)| Nuclear parameter | Value           |
| ------------- | ------------- | -------------     | --------------  |
| $\lambda_1$   |0.0124         | $\beta_1$         |0.00021          |
| $\lambda_2$   |0.0305         | $\beta_2$         |0.00141          |
| $\lambda_3$   |0.111          | $\beta_3$         |0.00127          |
| $\lambda_4$   |0.301          | $\beta_4$         |0.00255          |
| $\lambda_5$   |1.13           | $\beta_5$         |0.00074          |
| $\lambda_6$   |3.0            | $\beta_6$         |0.00027          |

with $\beta=0.00645$, $\Lambda=5 \times 10^{-5} \mathrm{s}$, and $a=0.01$ and $b=10^{-11}$. In this scenario is assummed that $\rho(0)=0$. The corresponding input and outputs are given in the following sections:

### Input

<details><summary>CLICK HERE to expand the input of the application of AnalyticNPKE-Feedback.py</summary>
<p>

```Python
#***********************Input Parameters ********************************
#************************************************************************
# L = list with lambda constants
# Betas = list with fractions of the precursors
# Lambda_m is the prompt neutron generation time
# Beta_total is computed by default
# a, b, parameters related to the Feedback reactivity

L=[0.0124, 0.0305, 0.111, 0.301, 1.13, 3]
Betas = [0.00021, 0.00141, 0.00127, 0.00255, 0.00074, 0.00027]
Lambda_m=0.00005
a = 0.01
b = 10**(-11)
#********************* Initial conditions *******************************
# n_0= neutron density at t=0
# C_init = list with the initial conditions for the precursors
n_0=1

#************* Initial conditions given by n(0)b_k/(Lambda_m*lambda_k)
for k in range(len(Betas)):
    C_init.append((Betas[k]/(L[k]*Lambda_m))*n_0)

#*************************************************************************
#********************* Discretization Input ******************************
#*************************************************************************
    
Target = 2                                            # Target is the desired time
step = 0.001                                          # Recommended value
#************************************************************************
```
</p>
</details>

### Outputs
As in the case of **AnalyticNPKE-Insertion.py**, there are two outputs. The one provided by the Python's shell where is given the maximum of the neutron density (also called "peak of the neutron density"), as well as the time where it is reached and the corresponding concentration of the precursors of the delayed neutrons. The second one consists of a ".txt" file where all the computed points are contained. 
### Output 1

```Python

Neutron peak: 20119823293.181976
Time of the peak: 1.106
Concentration of the precursors related to the neutron density peak
[1996512972.7807598, 13401179130.68308, 12054651574.611992, 24129095447.662556, 6908405306.398891, 2446315046.9759207]
```
which coincides with the one reported by Nahla (2010, p. 1626). 
### Output 2.
The "Solutions_feedback_final.txt" contains the solution of the precursors of the delayed neutrons for each step of calculation. An example of such file can be reviewed in the following link: [Example of an output file](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE/blob/main/Example_of_output_ramp.txt)

## 9. References.

1. Arfken, G. B., Weber, H. J., Harris, F. E. 2013. Mathematical Methods for Physicists, Seventh Edition, p. 1045.
2. Cruz-López, C.-A., Espinosa-Paredes, G., François, J.-L. 2023. A New Simplified Analytical Solution to Solve the Neutron Point Kinetics Equations Using the Laplace Transform Method. Computer Physics Communications, Vol. 238, 108564 https://doi.org/10.1016/j.cpc.2022.108564
3. Duderstadt, J. J., Hamilton, L. J. 1976. Nuclear Reactor Analysis. John Wiley & Sons, Inc.
4. Lewins, J. 1978. Nuclear Reactor Kinetics and Control. Pergamon Press.
5. Nahla, A. A. 2010. Analytical Solution to Solve the Point Rector Kinetic Equations. Nuclear Engineering and Design, 240, 1622-1629. https://doi.org/10.1016/j.nucengdes.2010.03.003 
6. Vinberg, E. B. 2003. A Course in Algebra. Graduated Studies in Mathematics, Vol. 56. American Mathematical Society. Providence, Rhode Island.  

