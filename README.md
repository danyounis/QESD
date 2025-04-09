# QESD
Quantum Entanglement Sudden Death

Source code for our publication:
S. Xie, D. Younis, and J. H. Eberly, "Evidence for unexpected robustness of multipartite entanglement against sudden death from spontaneous emission", Phys. Rev. Research **5**, L032015 (2023).

Primarily optimization routines and callers for the many objective functions therein.

A description of the taskx.f08 Fortran files is provided. The purpose of each
program is to compute the entanglement measure of a given mixed state, and they
are to be cross-referenced with Songbo's notes located in the previous directory.

Each program expects an input deck to be located in the current working directory (CWD)
with a fixed name corresponding to the task. Input decks are generally saved along with
their respective output datasets.

----------

task2: We have an analytic formula for the entanglement vs. time of the |GHZ> and |W> state mixture,
and, as a test of the method, we seek to reproduce it by numerical optimization.

task3: We have a |GHZ> state and know *a priori* that sudden death occurs at t~0.64. The task3.f08
program is actually a generalization of Songbo's notes; the state mixture there is obtained by
setting "Theta: 0.25 pi" in the input deck. Note the occurence of entanglement sudden death (ESD)
depends on the tuning of this state-mixture capital-Theta parameter, see the datasets:
task3/r3-x/Theta0.5rad.h5 and /Theta1.4rad.h5.

task4: The generalized |W> state, i.e. including global mixture parameter Theta.
ESD has, at the time of writing, not been observed, regardless of Theta's value.

task5m: The |W-bar> state, which is the spin-flipped |W> state. No Theta variable present.

task5v(1/2): The generalized |W-bar> state.
The |W-bar> state of task5m ("m" for "minus") is recovered by taking cos(Theta)=1/sqrt(3).
v1: The v-vectors of the Nelder-Mead simplex vertices are all initialized randomly.
v2: One v-vector is initialized randomly, and the rest are calculated by taking a fixed step in
each coordinate direction, vi = v0 + h*ei where h: step-size and ei: i-th unit vector.

task6: The |eee> + |W> state. (Sigma-1)

task7: The |eee> + |W-bar> state. (Sigma-2)

----------

task5-pso.tar.gz archive

These codes apply the particle swarm optimization algorithm to attempt to solve task5.
v1: Nelder-Mead exterior optimization / Particle Swarm interior optimization
v2: Particle Swarm exterior optimization / Nelder-Mead interior optimization

----------

Calculation codes ending in "a2" refer to a pure minimization procedure outlined
in "notes/Task 5 - 65 parameters.pdf". Note, the exact dimensionality of the
parameter-space varies from state to state, and is given by the formula:
    no. parameters = r×(2×k-r+1)
where r = density matrix rank and k = order + rank (where order >= 3).
Evidently, the greater the density matrix rank, the more challenging the minimization task.

We discovered that the "a2" algorithm is far superior to the nested max-min version.
All one must do is hard-code the eigenvalue/vector expressions for a particular state, and tune the
k parameter (k=rank+3 is the recommended minimum, which yields good results, but one can go higher).

task5a2: The generalized |W-bar> state, minimization over 65 parameters.
See "notes/Task 5 - 65 parameters.pdf".

task6a2: Pure minimization algorithm for the |Sigma-1> (|eee> + |W>) state.
task4a2: Pure minimization algorithm for the generalized |W> state.
