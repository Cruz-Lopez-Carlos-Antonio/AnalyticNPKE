# Analytical Solution of the Neutron Point Kinetic Equations 
The present repository contains the Python codes that were developed to solve the Neutron Point Kinetics Equations, using the Laplace transform and the Heaviside's Theorem. Such solution was reported in the paper *A new simplified analytical solution to solve the neutron point kinetics equations using the Laplace transform method*, publised in the the journal of *Computer Physics Communications* in the following [link](https://www.sciencedirect.com/science/article/abs/pii/S0010465522002831?via%3Dihub). 

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
$$\widetilde{n}\left(s\right)=\left(\frac{n\left(0\right)Q\left(s\right)}{P\left(s\right)}+\frac{H\left(s\right)}{P\left(s\right)}\right) \tag{7}$$

