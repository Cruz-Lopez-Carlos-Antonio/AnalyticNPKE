# Analytical Solution of the Neutron Point Kinetic Equations 
The present repository contains the Python codes that were developed to solve the Neutron Point Kinetics Equations, using the Laplace transform and the Heaviside's Theorem. Such solution was reported in the paper *A new simplified analytical solution to solve the neutron point kinetics equations using the Laplace transform method*, publised in the the journal of *Computer Physics Communications* in the following [link](https://www.sciencedirect.com/science/article/abs/pii/S0010465522002831?via%3Dihub). 

The programs are licensed under a Creative Commons Attribution 4.0 International License: http://creativecommons.org/licenses/by/4.0/

Authors: Carlos-Antonio Cruz-López (cacl.nucl@gmail.com), Gilberto Espinosa-Paredes (gepe@xanum.uam.mx)

Mathematical and algorithmical generalities of the codes are described in the following lines with the purpose to provide some insight of the developed work. Nevertheless, a more detailed and precise discussion is provided in the submitted article.
## Financial Support.
The authors appreciate the financial support received from the Consejo Nacional de Ciencia y Tecnología, CONAHCYT, under the program “Estancias Posdoctorales por México, 2022”, with the project entitled: “Desarrollo de modelos fenomenológicos energéticos de orden fraccional, para la optimización y simulación en reactores nucleares de potencia”, by which the present development was possible.

## 1. Mathematical description of the problem
The Neutron Point Kinetic Equations (NPKE) with $K$ groups of precursors of delayed neutrons can be written as follows:
$$\frac{dC_k\left(t\right)}{dt}=\frac{\beta_k}{l}n\left(t\right)-\lambda_kC_k\left(t\right),\ \ 1\le\ k\le\ K&&
