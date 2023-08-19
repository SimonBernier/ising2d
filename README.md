# ising2d
Itensor analysis of two-dimension transverse field Ising (2D-TFI) model

criticality calculates correlations at the critical point of the 2D-TFI model

fidelity calculates the quantum fidelity of 2D-TFI, aimed at finding the critical point and critical exponent via scaling

gap calculates the energy gap of 2D-TFI, aimed at finding the critical point and critical exponents via scaling

ising-tebd is a test to implement 2D TEBD with smart ordering of swap gates

movingFront is aimed at calculating the energy, spin correlations and von Neumann entanglement entropy during an inhomogeneous quench in the 2D-TFI with a quench front moving at constant velocity. In this version, the program keeps track of the instantaneous ground state. The time evolution is carried out using a 4th order time dependent variational principle.

superluminal is similar to movingFront except the program does not keep track of the instantaneous ground state. It is in principle a faster but less detailed implementation of moving quench fronts in the 2D-TFI model

uniform calculates the energy, spin correlations and von Neumann entanglement entropy during a homogeneous quench in the 2D-TFI. It is aimed at comparing the effieciency of our inhomogenous technique to uniform quenches. The time evolution is carried out using a 4th order time dependent variational principle.

vCrit calculates the von Neumann entanglement entropy after a local perturbation in the ground state of the 2D-TFI model. It is aimed at calculating the speed of light of the model. The time evolution is carried out using a 4th order time dependent variational principle.
