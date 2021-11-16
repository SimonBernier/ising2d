Roadmap to building the code.

Code works for time independent Hamiltonian, which means the ground state just accumulates a phase and stays in ground state. Previous work shows that the gap in energy closes near 2.0 (1D) and 2.3(2D with Ny=10).

Add time dependent Hamiltonian which goes from h = 4 to h = 2 thereby closing the gap. The time step will likely have to be small, but first, compare the DensityMatrix and Fit methods.

Important to get this right because we can then apply everything we develop here for a simple 2D system for the more complicated Heisenberg Hamiltonian.