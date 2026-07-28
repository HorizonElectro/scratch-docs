# scratch

```
Alex Laut
Feb 2026
```

## Introduction

`scratch` solves the electric field integral equation (EFIE) for scattering from arbitrary
perfectly-electric-conducting (PEC) surfaces, by the method of moments (MoM) with
Rao–Wilton–Glisson (RWG) basis functions [1].

## The physical problem

An incident field $\mathbf{E}^i$ illuminates a PEC surface $S$. The boundary condition — the
tangential electric field vanishes on a perfect conductor — relates the (unknown) induced surface
current $\mathbf{J}$ on $S$ to $\mathbf{E}^i$ through the retarded mixed-potential integral
equation:

$$
\hat{\mathbf n}\times\mathbf E^i(\mathbf r) = \hat{\mathbf n}\times\Big[j\omega\mathbf A(\mathbf r) + \nabla\Phi(\mathbf r)\Big], \qquad \mathbf r \in S
$$

with the vector and scalar potentials

$$
\mathbf A(\mathbf r) = \mu \int_S \mathbf J(\mathbf r')\, G(\mathbf r, \mathbf r')\, dS'
$$

$$
\Phi(\mathbf r) = \frac{j}{\omega\varepsilon} \int_S \big[\nabla_s'\cdot \mathbf J(\mathbf r')\big]\, G(\mathbf r, \mathbf r')\, dS'
$$

(the scalar potential follows from the surface continuity equation $\nabla_s\cdot\mathbf J = -j\omega\sigma$, i.e. $\sigma = \frac{j}{\omega}\nabla_s\cdot\mathbf J$), and the free-space Green's function

$$
G(\mathbf r, \mathbf r') = \frac{e^{-jk|\mathbf r - \mathbf r'|}}{4\pi|\mathbf r - \mathbf r'|}, \qquad k = \frac{2\pi}{\lambda} = \frac{\omega}{c}.
$$

## Discretization: RWG basis functions

$S$ is triangulated, and the current is expanded over the mesh's interior edges,

$$
\mathbf J(\mathbf r) = \sum_{n=1}^{N} I_n\, \mathbf f_n(\mathbf r),
$$

using the Rao–Wilton–Glisson basis function of edge $n$ [1]. Each interior edge is shared by two
triangles $T_n^+, T_n^-$ with respective free (non-shared) vertices $\mathbf r_n^+, \mathbf r_n^-$
and areas $A_n^+, A_n^-$; with $l_n$ the shared edge length,

$$
\mathbf f_n(\mathbf r) =
\begin{cases}
\dfrac{l_n}{2A_n^+}\big(\mathbf r - \mathbf r_n^+\big), & \mathbf r \in T_n^+ \\[8pt]
\dfrac{l_n}{2A_n^-}\big(\mathbf r_n^- - \mathbf r\big), & \mathbf r \in T_n^- \\[6pt]
0, & \text{otherwise},
\end{cases}
\qquad
\nabla_s\cdot \mathbf f_n(\mathbf r) =
\begin{cases}
\dfrac{l_n}{A_n^+}, & \mathbf r \in T_n^+ \\[6pt]
-\dfrac{l_n}{A_n^-}, & \mathbf r \in T_n^-.
\end{cases}
$$

$\mathbf f_n$ is continuous across the shared edge, divergence-free everywhere else, and represents
a current of uniform total flow $l_n$ across that edge — the RWG basis exists precisely so surface
current can be expanded edge-by-edge without introducing spurious line charges at triangle
boundaries.

## The linear system

Testing the EFIE with the same basis, $\langle \mathbf f_m, \cdot\rangle$ (evaluated here at each
testing function's own $T_m^+/T_m^-$ centroids $\mathbf r_{c,m}^{\pm}$), turns the integral equation
into a dense linear system

$$
\mathsf Z \mathbf I = \mathbf V,
$$

$$
Z_{mn} = l_m\left[ j\omega\left(\frac{\mathbf A_{mn}(\mathbf r_{c,m}^+)\cdot \boldsymbol\rho_{c,m}^+}{2} + \frac{\mathbf A_{mn}(\mathbf r_{c,m}^-)\cdot \boldsymbol\rho_{c,m}^-}{2}\right) + \Phi_{mn}(\mathbf r_{c,m}^-) - \Phi_{mn}(\mathbf r_{c,m}^+) \right]
$$

$$
V_m = l_m\left(\frac{\mathbf E^i(\mathbf r_{c,m}^+)\cdot \boldsymbol\rho_{c,m}^+}{2} + \frac{\mathbf E^i(\mathbf r_{c,m}^-)\cdot \boldsymbol\rho_{c,m}^-}{2}\right)
$$

where $\mathbf A_{mn}, \Phi_{mn}$ denote $\mathbf A, \Phi$ due to source basis $\mathbf f_n$ alone,
and $\boldsymbol\rho_{c,m}^+ = \mathbf r_{c,m}^+ - \mathbf r_m^+$,
$\boldsymbol\rho_{c,m}^- = \mathbf r_m^- - \mathbf r_{c,m}^-$ are each testing triangle's
centroid-to-free-vertex vectors [1, eq. 17–20]. $N$ (the number of interior mesh edges) unknown
current coefficients $I_n$ are solved for via GMRES.

## Radiated field reconstruction

Given the solved $I_n$ (and hence $\mathbf J$ and, via continuity, $\sigma$), the scattered field at
any point in free space — on the mesh itself or at a remote probe — follows from the same
mixed-potential form, now evaluated off the boundary (so no tangential projection is needed):

$$
\mathbf E(\mathbf r) = -j\omega\mathbf A(\mathbf r) - \nabla\Phi(\mathbf r).
$$

## References

1. S. M. Rao, D. R. Wilton, and A. W. Glisson, "Electromagnetic scattering by surfaces of arbitrary
   shape," *IEEE Transactions on Antennas and Propagation*, vol. AP-30, no. 3, pp. 409–418, May 1982.
2. R. F. Harrington, *Field Computation by Moment Methods*. Macmillan, 1968 (reissued, IEEE Press,
   1993) — the general method-of-moments framework RWG basis functions are built on.
3. W. C. Gibson, *The Method of Moments in Electromagnetics*, 2nd ed. CRC Press, 2014 — a modern
   treatment of RWG-based MoM, including implementation-level detail (matrix assembly, radiation
   integrals, fast solvers).
