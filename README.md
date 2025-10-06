# 2D Transverse Field Ising Model Analysis

Tensor network simulations of the two-dimensional transverse field Ising (2D-TFI) model using ITensor. This project explores equilibrium properties, quantum quenches, and entanglement dynamics in strongly correlated quantum systems.

## 🔬 Overview

This repository contains numerical implementations for studying quantum phase transitions and out-of-equilibrium dynamics in the 2D transverse field Ising model. The simulations use tensor network methods, specifically the Time-Dependent Variational Principle (TDVP), to efficiently simulate quantum many-body systems.

### Physical Model

The 2D transverse field Ising Hamiltonian:

```
H = -J Σ σˣᵢ σˣⱼ - Σ hᵢ σᶻᵢ
```

where J is the coupling strength, hᵢ is the local transverse field, and the sums run over nearest neighbors and all sites respectively.

---

## 📂 Project Structure

### Core Modules

#### 1. **gap** - Equilibrium Properties
Calculates equilibrium properties of the 2D-TFI model, including:
- Energy gap between ground state and first excited state (critical point identification via gap scaling)
- Ground state energy
- Phase transition characterization

**Key Features:**
- Finite-size scaling analysis
- Critical exponent extraction
- Quantum phase transition mapping

---

#### 2. **vCrit** - Speed of Excitations
Determines the effective speed of light in the quantum system using entanglement dynamics.

**Method:**
- Introduces local perturbation to ground state wavefunction
- Evolves disturbed state using 4th order TDVP
- Tracks von Neumann entropy growth
- Extracts propagation velocity from entanglement spreading

**Physical Insight:** The speed at which entanglement propagates reveals fundamental information about the model's causal structure and emergent light cone.

---

#### 3. **uniform** - Homogeneous Quantum Quenches
Simulates uniform (homogeneous) quenches in the 2D-TFI model.

**Protocol:**
- Initial state: Ground state in gapped phase
- Quench: Instantaneous change to critical point
- Evolution: 4th order TDVP time evolution

**Observables:**
- Energy density evolution
- Spin correlations
- Von Neumann entanglement entropy dynamics

**Purpose:** Provides baseline for comparing with inhomogeneous quench dynamics.

---

#### 4. **movingFront** - Inhomogeneous Quenches ⭐
The main focus of this project - simulates spatially inhomogeneous quenches with moving quench fronts.

**Key Features:**
- Quench front moves at constant velocity through the system
- Tracks instantaneous ground state during evolution
- 4th order TDVP for accurate time evolution
- Detailed tracking of non-equilibrium dynamics

**Observables:**
- Energy density distribution
- Spin correlation functions
- Von Neumann entanglement entropy
- Spatial structure of excitations

**Physical Motivation:** Studies how quantum information and excitations propagate when parameters change locally and spread through the system - relevant for understanding causality and out-of-equilibirum dynamics in quantum systems.

---

### Additional Modules

- **criticality**: Correlation functions at the critical point
- **fidelity**: Quantum fidelity calculations for critical point detection
- **ising-tebd**: Test implementation of 2D Time-Evolving Block Decimation (TEBD) with optimized swap gate ordering
- **superluminal**: Streamlined version of movingFront without instantaneous ground state tracking (faster but less detailed)

---

## 🛠️ Technical Implementation

**Framework:** [ITensor](https://itensor.org/) - C++ library for tensor network calculations

**Key Algorithms:**
- 4th order Time-Dependent Variational Principle (TDVP)
- Matrix Product State (MPS) representation
- Finite-size scaling analysis
- Entanglement entropy calculations

**Language:** C++

---

## 📊 Applications

This work contributes to understanding:
- **Quantum phase transitions** in 2D systems
- **Non-equilibrium dynamics**
- **Entanglement propagation** and quantum information spreading
- **Kibble-Zurek mechanism** in quantum quenches
- **Emergent causality** in quantum many-body systems

---

## 🎓 Academic Context

This project was developed as part of graduate research in computational quantum many-body physics, focusing on:
- Tensor network methods for 2D systems
- Out-of-equilibrium quantum dynamics
- Quantum criticality and universality
- Numerical methods for strongly correlated systems

---

## 📚 Related Physics Concepts

- **Quantum Phase Transitions:** The 2D-TFI model exhibits a quantum phase transition between ordered and disordered phases
- **Entanglement Entropy:** Quantifies quantum correlations and serves as a diagnostic for phase transitions
- **TDVP:** Variational method for time evolution that preserves the tensor network structure
- **Kibble-Zurek Dynamics:** Universal scaling behavior when crossing phase transitions

---

## 🔗 Dependencies

- ITensor library
- C++ compiler with C++11 support
- LAPACK/BLAS libraries (typically handled by ITensor)

---

## 📖 Usage

Each module can be compiled independently. Typical workflow:

1. Set physical parameters (system size, coupling strengths, quench protocol)
2. Initialize ground state or desired initial state
3. Run simulation (equilibrium calculation or time evolution)
4. Output observables for analysis

---

## 📈 Future Directions

Potential extensions:
- Higher-dimensional generalizations
- Different lattice geometries
- Additional observables (mutual information, correlation length)
- Floquet engineering with periodic driving
- Connections to machine learning and tensor network quantum states

---

## 📧 Contact

**Simon Bernier**
- Email: simon.bernier@mail.mcgill.ca
- LinkedIn: [simon-bernier-6701a9285](https://www.linkedin.com/in/simon-bernier-6701a9285)

---

## 📝 Citation

If you use this code in your research, please cite appropriately and feel free to reach out with questions or collaboration ideas.

---

*This project demonstrates expertise in: computational physics, tensor networks, C++ programming, quantum many-body systems, and numerical methods for strongly correlated systems.*
