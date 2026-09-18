# Literature review: three-dimensional finite-element TDGL

## Scope and reading method

This review supports a production-quality MATLAB solver for a superconducting body coupled self-consistently to its electromagnetic surroundings. The local library contains legally accessible author manuscripts, arXiv preprints, and institutional-repository copies. The acquisition manifest is in [`papers/README.md`](../papers/README.md), and machine-readable citations are in [`papers/references.bib`](../papers/references.bib). Sources without a legal open full text are recorded by DOI and abstract/catalog page rather than copied.

The literature separates into four strands that are not interchangeable:

1. TDGL analysis and discretization on the superconducting domain;
2. three-dimensional computational TDGL, mostly structured finite differences;
3. compatible finite elements for Maxwell fields and nonsmooth topology; and
4. exterior-domain methods for magnetoquasistatic or wave Maxwell problems.

No single paper supplies the complete method required here. The recommended solver therefore combines complex nodal elements for the order parameter, an exact-sequence edge-element Maxwell discretization, a mixed gauge constraint, and an explicitly modelled vacuum exterior.

## Method comparison at a glance

“GL” below means the gauge-covariant order-parameter equation coupled to the stated magnetic/current equation. A dash means the paper does not specify that component because it studies a different subproblem.

| Source | Equations and gauge | Spatial discretization | Boundary/exterior treatment | Time integration |
|---|---|---|---|---|
| Hong et al. 2023 | TDGL, temporal gauge | complex Lagrange $H^1$ for $\psi$; second-family Nédélec $H(\mathrm{curl})$ for $\mathbf A$ | bounded sample; natural superconducting boundary data; no all-space exterior | nonlinear backward step; Newton with block preconditioner |
| Gao–Sun 2018 | TDGL, Lorenz gauge; $\boldsymbol\sigma=\nabla\times\mathbf A$ mixed form | Lagrange $H^1$, first-family Nédélec $H(\mathrm{curl})$, Raviart–Thomas $H(\mathrm{div})$ | simply connected bounded Lipschitz/polyhedral domain | linearized backward Euler |
| Li–Zhang 2015/2017 | TDGL, Lorenz gauge recast by 2-D Hodge potentials | nodal scalar FEM after Hodge decomposition | curved/nonconvex polygons; only simply connected equivalence | linearized implicit step |
| Du–Gunzburger–Peterson 1992 | stationary GL, gauge-equivalent formulations | continuous quadratic nodal FEM | 2-D periodic cell | equilibrium iteration/minimization |
| Du 1994/1998 | TDGL; gauge-covariant/discrete gauge-invariant variants | Galerkin FEM or link/difference formulation | bounded model domains | semidiscrete and fully discrete schemes |
| Gunter–Kaper–Leaf 2002 | TDGL with link-variable gauge consistency | structured differences | superconducting region with thin insulating setting | four methods from explicit to fully nonlinear implicit |
| Winiecki–Adams 2002 | 3-D TDGL + MQS normal/self current | Cartesian link-variable differences | finite wire; self-field iterated by Biot–Savart | semi-implicit Crank–Nicolson and fractional steps |
| Sadovskyy et al. 2015 | 3-D TDGL; $\nabla\cdot\mathbf A=0$-type convention; scalar-potential current constraint | structured GPU finite differences/link variables | open covariant Neumann or quasiperiodic; usually prescribed applied $\mathbf A$ | implicit Crank–Nicolson, linearized cubic, Jacobi iteration |
| Oripov–Anlage 2020 | 3-D TDGL in S + Maxwell in V; scalar potential set by chosen formulation | COMSOL general nodal PDE on free tetrahedra | no-normal-current S–V interface; applied field on remote vacuum boundary | COMSOL BDF, direct MUMPS, maximum nondimensional step 1 |
| Lara et al. 2020 | 3-D TDGL, temporal gauge | Cartesian link variables | field imposed on sample boundary; no demagnetizing exterior | explicit/iterative evolution as reported by implementation |
| Doria et al. 2006 | stationary 3-D GL | structured variational discretization | de Gennes condition on insulating spheres | energy minimization, not TDGL time stepping |
| Alstrøm et al. 2011 | 2-D TDGL, temporal gauge | COMSOL nodal FEM | de Gennes/no-current condition; applied field on sample boundary | implicit COMSOL integration |
| Chen–Dai 2001 | 2-D TDGL with adaptive error control | continuous P1 Galerkin FEM | bounded domain | semi-implicit evolution with space/time estimators |
| Gao–Li–Sun 2014 | TDGL Galerkin formulation | conforming nodal FEM under paper assumptions | bounded domain | linearized Crank–Nicolson, second order |
| Feischl–Tran 2016/2017 | MQS eddy current + LLG, not TDGL | interior 3-D FEM/Nédélec plus BEM traces | exact unbounded electromagnetic exterior | linear systems per implicit time step; convergence/error analysis |
| Bao et al. 2023 | time-harmonic full-wave Maxwell, no TDGL | 3-D edge FEM | spherical DtN radiation boundary with truncation estimator | frequency-domain, no transient step |
| Chen–Guo–Zou 2022 | Maxwell interface problem, no TDGL | immersed Nédélec-type $H(\mathrm{curl})$ and de Rham spaces | unfitted material interface | static/frequency-domain linear solve |
| Xue et al. 2024 | review of several TDGL variants/applications | surveys finite-difference/FEM practices | surveys pinning, SRF, and material configurations | review, not a new integrator |

This table is deliberately explicit about omissions: a paper that has excellent time integration but prescribes $\mathbf A$ does not validate an exterior self-field solve, and a rigorous exterior Maxwell paper does not validate superconducting dynamics.

## Governing-model sources

### Schmid (1966) and Gor'kov–Eliashberg (1968)

**Citations.** A. Schmid, “A time dependent Ginzburg–Landau equation and its application to the problem of resistivity in the mixed state,” *Physik der kondensierten Materie* 5, 302–317 (1966), [DOI: 10.1007/BF02422669](https://doi.org/10.1007/BF02422669). L. P. Gor'kov and G. M. Eliashberg, “Generalization of the Ginzburg–Landau equations for non-stationary problems in the case of alloys with paramagnetic impurities,” *Soviet Physics JETP* 27, 328–334 (1968), [catalog record](https://www.jetp.ras.ru/cgi-bin/e/index/e/27/2/p328?a=list).

These are the historical foundations for dissipative, gauge-covariant TDGL near the critical temperature. The order parameter evolves by a covariant time derivative and couples to electric and magnetic potentials; the electromagnetic equation contains superconducting and normal current. They justify the phenomenological model but do not prescribe a modern 3-D finite-element space, exterior truncation, or robust nonlinear solver. Their microscopic validity is restricted: conventional gapless/dirty superconductors near $T_c$ are the natural regime, and quantitative low-temperature dynamics or strong nonequilibrium physics needs a more microscopic theory.

### Xue et al. (2024): applications review

**Citation.** C. Xue, Q.-Y. Wang, H.-X. Ren, A. He, and A. V. Silhanek, “Case studies on time-dependent Ginzburg–Landau simulations for superconducting applications,” [arXiv:2403.03729](https://arxiv.org/abs/2403.03729) ([local PDF](../papers/xue-et-al-2024-tdgl-review.pdf)).

This modern review organizes applications around vortex ratchets/diodes, pinning and critical current, and superconducting RF structures. It is useful for locating contemporary validation targets and for seeing how material inhomogeneity is represented in large-scale TDGL. It is not a new edge-element FEM formulation, and its broad case-study scope does not resolve the exterior-field or topology questions by itself.

### Full electromagnetic-wave TDGL

**Citations.** C. Fan and M. Ozawa, “The time-dependent Ginzburg–Landau–Maxwell equations,” Kyoto University preprint/analysis page (2017), [repository search](https://repository.kulib.kyoto-u.ac.jp/). Y. Tsutsumi and H. Kasai, “The time-dependent Ginzburg–Landau–Maxwell equations,” *Nonlinear Analysis* 37, 187–216 (1999), [DOI: 10.1016/S0362-546X(98)00043-1](https://doi.org/10.1016/S0362-546X(98)00043-1).

These works retain the displacement-current term, schematically

$$
\epsilon(\partial_{tt}\mathbf A+\nabla\partial_t\phi)
+\sigma(\partial_t\mathbf A+\nabla\phi)
+\nabla\times\mu^{-1}\nabla\times\mathbf A+\mathbf J_s=\mathbf J_{\rm ext}.
$$

They establish that “TDGL–Maxwell” can mean a hyperbolic electromagnetic-wave system, not merely eddy-current coupling. For vortex motion, DC transport, and dimensions much smaller than the electromagnetic wavelength, the magnetoquasistatic (MQS) omission of the first term is normally the appropriate baseline. RF radiation, resonant cavities, or wavelength-scale devices require the full-wave branch, its radiation boundary treatment, and a different time-step restriction. The two models will be kept distinct in this project.

## TDGL finite-element formulations

### Hong, Ma, Xu, and Chen (2023): direct edge-element TDGL

**Citation.** Q. Hong, L. Ma, J. Xu, and L.-Q. Chen, “An efficient iterative method for dynamical Ginzburg–Landau equations,” *Journal of Computational Physics* 474, 111794 (2023), [DOI: 10.1016/j.jcp.2022.111794](https://doi.org/10.1016/j.jcp.2022.111794), [arXiv](https://arxiv.org/abs/2207.01425), [local PDF](../papers/hong-ma-xu-2022-nedelec-tdgl.pdf).

The paper works in temporal gauge and discretizes the complex order parameter with $H^1$-conforming Lagrange elements and the magnetic vector potential with lowest-order, second-family Nédélec edge elements in $H(\mathrm{curl})$. Backward time stepping produces a nonlinear system solved by Newton iteration; a block preconditioner is developed. Numerical tests include nonsmooth/reentrant geometries and show the practical advantage of a curl-conforming vector potential.

**Relevance.** This is the closest direct foundation for the superconducting part of the proposed 3-D method: it preserves the natural Maxwell regularity and avoids falsely assuming $\mathbf A\in H^1$ on nonconvex domains.

**Limitations.** The electromagnetic problem is posed on the computational/sample domain with prescribed boundary data, not on an unbounded vacuum exterior. Temporal gauge is convenient for magnetic-only tests but does not by itself provide a physical scalar potential for voltage-driven transport. Multiply connected harmonic fields and S–N proximity interfaces are not the main subject.

### Gao and Sun (2018): mixed Lorenz-gauge formulation

**Citation.** H. Gao and W. Sun, “Analysis of linearized Galerkin-mixed FEMs for the time-dependent Ginzburg–Landau equations of superconductivity,” *Advances in Computational Mathematics* 44, 923–949 (2018), [DOI: 10.1007/s10444-017-9568-2](https://doi.org/10.1007/s10444-017-9568-2), [arXiv](https://arxiv.org/abs/1508.05601), [local PDF](../papers/gao-sun-2016-galerkin-mixed-tdgl.pdf).

The authors introduce $\boldsymbol\sigma=\nabla\times\mathbf A$ and use a de Rham-compatible triple: Lagrange $H^1$ elements for $\psi$, first-family Nédélec $H(\mathrm{curl})$ elements for $\boldsymbol\sigma$, and Raviart–Thomas $H(\mathrm{div})$ elements for $\mathbf A$. A linearized backward-Euler scheme is proved unconditionally stable with optimal error estimates in two and three dimensions.

**Relevance.** It supplies a rigorous mixed alternative, exposes magnetic flux as a primary variable, and is attractive for local current/flux conservation.

**Limitations.** The analysis assumes a simply connected bounded Lipschitz/polyhedral domain and a particular Lorenz-gauge formulation. It does not solve the all-space magnetic field, add transport terminals, or resolve cohomology modes of multiply connected bodies.

### Li and Zhang (2015, 2017): nonsmooth-domain warning

**Citations.** B. Li and Z. Zhang, “A new approach for numerical simulation of the time-dependent Ginzburg–Landau equations,” *Journal of Computational Physics* 303, 238–250 (2015), [DOI: 10.1016/j.jcp.2015.09.049](https://doi.org/10.1016/j.jcp.2015.09.049), [arXiv](https://arxiv.org/abs/1410.3746), [local PDF](../papers/li-zhang-2014-hodge-tdgl.pdf). B. Li and Z. Zhang, “Mathematical and numerical analysis of the time-dependent Ginzburg–Landau equations in nonconvex polygons,” *Mathematics of Computation* 86, 1579–1608 (2017), [DOI: 10.1090/mcom/3177](https://doi.org/10.1090/mcom/3177), [arXiv](https://arxiv.org/abs/1410.3547), [local PDF](../papers/li-zhang-2017-nonconvex-tdgl.pdf).

These papers demonstrate that the magnetic potential in a nonconvex domain need not have $H^1$ regularity. Ordinary nodal-vector FEM can converge to a spurious solution. In two dimensions they use Hodge decomposition to replace the vector potential by scalar potentials and employ stable nodal elements and linearized time stepping.

**Relevance.** They are decisive evidence against assembling three independent nodal components of $\mathbf A$ for arbitrary 3-D geometries.

**Limitations.** The scalar Hodge reduction is inherently two-dimensional. In three dimensions the Hodge component is vector-valued; in multiply connected domains nontrivial harmonic fields appear and the simplified system is no longer equivalent unless cycle constraints are supplied. The papers motivate Nédélec/mixed elements rather than provide the final 3-D algorithm.

### Du, Gunzburger, and Peterson (1992); Du (1994, 1998)

**Citations.** Q. Du, M. D. Gunzburger, and J. S. Peterson, “Solving the Ginzburg–Landau equations by finite-element methods,” *Physical Review B* 46, 9027–9034 (1992), [DOI: 10.1103/PhysRevB.46.9027](https://doi.org/10.1103/PhysRevB.46.9027), [local PDF](../papers/du-gunzburger-peterson-1992-gl-fem.pdf). Q. Du, “Finite element methods for the time-dependent Ginzburg–Landau model of superconductivity,” *Computers & Mathematics with Applications* 27(12), 119–133 (1994), [DOI: 10.1016/0898-1221(94)90091-4](https://doi.org/10.1016/0898-1221(94)90091-4). Q. Du, “Discrete gauge invariant approximations of a time dependent Ginzburg–Landau model of superconductivity,” *Mathematics of Computation* 67, 965–986 (1998), [DOI: 10.1090/S0025-5718-98-00954-5](https://doi.org/10.1090/S0025-5718-98-00954-5).

The 1992 paper analyzes stationary GL minimization and finite-element approximation, including gauge-equivalent formulations and periodic-cell magnetization calculations. The later work develops time-dependent Galerkin and discrete gauge-invariant approximations. These sources establish the variational/free-energy viewpoint and gauge invariance as a numerical property to test.

**Limitations.** The benchmark computations are chiefly two-dimensional, periodic, or bounded-domain. Continuous nodal vector potentials from early formulations are not robust enough for general 3-D nonconvex Maxwell fields. They also predate modern block preconditioning and explicit exterior coupling.

### Chen and Dai (2001); Gao, Li, and Sun (2014)

**Citations.** Z. Chen and S. Dai, “Adaptive Galerkin methods with error control for a dynamical Ginzburg–Landau model in superconductivity,” *SIAM Journal on Numerical Analysis* 38, 1961–1985 (2001), [DOI: 10.1137/S0036142998349102](https://doi.org/10.1137/S0036142998349102). H. Gao, B. Li, and W. Sun, “Optimal error estimates of linearized Crank–Nicolson Galerkin FEMs for the time-dependent Ginzburg–Landau equations,” *SIAM Journal on Numerical Analysis* 52, 1183–1202 (2014), [DOI: 10.1137/130918678](https://doi.org/10.1137/130918678).

Chen–Dai combine continuous finite elements, semi-implicit evolution, and dual-weighted/a posteriori error control in two dimensions. Gao–Li–Sun establish optimal error estimates for a linearized Crank–Nicolson method under a mild time-step condition rather than a mesh-dependent CFL restriction.

**Relevance.** They support residual adaptivity and second-order semi-implicit time integration after a conservative backward-Euler baseline is validated.

**Limitations.** Neither delivers a complete 3-D $H(\mathrm{curl})$ exterior-Maxwell implementation, and the 2-D estimators cannot simply be transplanted to edge elements.

## Three-dimensional computational TDGL

### Gunter, Kaper, and Leaf (2002)

**Citation.** D. O. Gunter, H. G. Kaper, and G. K. Leaf, “Implicit integration of the time-dependent Ginzburg–Landau equations of superconductivity,” *SIAM Journal on Scientific Computing* 23, 1943–1958 (2002), [DOI: 10.1137/S1064827500375473](https://doi.org/10.1137/S1064827500375473), [arXiv](https://arxiv.org/abs/math/9906176), [local PDF](../papers/gunter-kaper-leaf-2002-implicit-tdgl.pdf).

Four schemes ranging from explicit to nonlinearly implicit are compared for vortex equilibration, with a superconductor embedded in an insulating setting. The paper demonstrates the stability benefit and nonlinear-solve cost of implicit treatment.

**Relevance.** It motivates a robust fully implicit reference integrator against which cheaper schemes can be checked.

**Limitations.** The spatial method and test configuration do not solve the general tetrahedral, multiply connected, all-space Maxwell problem.

### Winiecki and Adams (2002)

**Citation.** T. Winiecki and C. S. Adams, “A fast semi-implicit finite-difference method for the TDGL equations,” *Journal of Computational Physics* 179, 127–139 (2002), [DOI: 10.1006/jcph.2002.7047](https://doi.org/10.1006/jcph.2002.7047), [arXiv](https://arxiv.org/abs/cond-mat/0106466), [local PDF](../papers/winiecki-adams-2002-semi-implicit-tdgl.pdf).

This is a practical three-dimensional, gauge-consistent link-variable finite-difference method. A semi-implicit Crank–Nicolson/fractional-step treatment permits much larger time steps than explicit updates. It includes a current-carrying wire and iterated Biot–Savart self-field calculation while neglecting displacement current.

**Relevance.** It supplies useful transport and self-field benchmarks and reinforces the MQS model choice.

**Limitations.** Cartesian differences do not represent curved/nonconvex geometries efficiently. The global Biot–Savart iteration is costly and is not a substitute for a scalable exterior FEM/BEM solve.

### Sadovskyy et al. (2015)

**Citation.** I. A. Sadovskyy et al., “Stable large-scale solver for Ginzburg–Landau equations for superconductors,” *Journal of Computational Physics* 294, 639–654 (2015), [DOI: 10.1016/j.jcp.2015.04.002](https://doi.org/10.1016/j.jcp.2015.04.002), [arXiv](https://arxiv.org/abs/1409.8340), [local PDF](../papers/sadovskyy-et-al-2015-large-scale-tdgl.pdf).

The paper develops a GPU-oriented 3-D structured-grid solver using link variables, Landau/Coulomb-like gauge conventions, a scalar-potential Poisson equation from current continuity, and implicit Crank–Nicolson with a linearized cubic term and Jacobi iteration. Material inhomogeneity enters through a local critical-temperature coefficient; insulating inclusions use no-current internal boundaries. Quasiperiodic conditions and thermal noise are supported.

**Relevance.** It is a strong performance and bulk-pinning reference, and its treatment of current conservation is useful for transport.

**Limitations.** The magnetic potential is often prescribed rather than solved with the finite specimen and vacuum, so demagnetizing/self-screening effects are outside its usual operating mode. A Cartesian grid also makes thin curved layers expensive.

### Oripov and Anlage (2020): explicit superconductor–vacuum coupling

**Citation.** B. Oripov and S. M. Anlage, “Time-dependent Ginzburg–Landau treatment of rf magnetic vortices in superconductors: Vortex semiloops in a spatially nonuniform magnetic field,” *Physical Review E* 101, 033306 (2020), [DOI: 10.1103/PhysRevE.101.033306](https://doi.org/10.1103/PhysRevE.101.033306), [legal APS accepted manuscript](https://link.aps.org/accepted/10.1103/PhysRevE.101.033306), [local PDF](../papers/oripov-anlage-2020-rf-vortices.pdf).

This COMSOL study solves coupled TDGL in a three-dimensional superconducting region and Maxwell equations in an adjoining vacuum region. The applied field is placed on the remote vacuum boundary. At the S–vacuum interface it uses no normal order-parameter/current flow and electromagnetic matching. It uses free tetrahedra, a direct MUMPS solve, and COMSOL's backward differentiation formula time integrator. Validation includes a superconducting sphere in a uniform field and screening beneath a point magnetic dipole; it then computes RF-driven vortex semiloops.

**Relevance.** This paper gives unusually direct evidence for the project's central physical constraint. Its single-domain comparison fixes the applied field on the superconducting surface and misses equatorial field enhancement, whereas the two-domain vacuum model agrees with the analytical sphere solution.

**Limitations.** The formulation is implemented through general COMSOL nodal PDE variables rather than a documented Nédélec exact-sequence discretization. It uses a finite vacuum truncation and does not address multiply connected harmonic modes or FEM–BEM. Its chosen boundary representation for $\mathbf A$ must therefore be rederived in the proposed mixed $H(\mathrm{curl})$ weak form rather than copied componentwise.

### Lara et al. (2020)

**Citation.** A. Lara et al., “Three-dimensional time-dependent Ginzburg–Landau simulations of superconducting samples,” *Low Temperature Physics* 46, 316–324 (2020), [DOI: 10.1063/10.0000861](https://doi.org/10.1063/10.0000861), [arXiv](https://arxiv.org/abs/2001.07971), [local PDF](../papers/lara-et-al-2020-3d-tdgl.pdf).

The authors use 3-D link-variable finite differences in temporal gauge and demonstrate vortex structures in finite samples. The implementation applies magnetic data at the sample boundary and discusses topology limitations.

**Relevance.** The examples are valuable for qualitative 3-D vortex visualization.

**Limitations.** Applying the ambient field on the superconductor boundary suppresses the distinction between applied and self-consistent screening fields and omits finite-sample demagnetization. Holes/multiply connected domains are not fully supported. This is exactly the shortcut the present solver must avoid.

### Doria et al. (2006)

**Citation.** M. M. Doria, A. R. de C. Romaguera, and W. A. M. Morgado, “Three-dimensional Ginzburg–Landau simulation of a vortex line displaced by a zigzag of pinning spheres,” *Pramana* 66, 119–127 (2006), [DOI: 10.1007/BF02704958](https://doi.org/10.1007/BF02704958), [arXiv](https://arxiv.org/abs/cond-mat/0503691), [local PDF](../papers/doria-et-al-2006-3d-pinning.pdf).

This stationary 3-D GL study models insulating spherical cavities with de Gennes/no-normal-supercurrent conditions and analyzes vortex detachment and free-energy changes.

**Relevance.** It supplies a geometry/material benchmark for insulating inclusions and vortex-line extraction.

**Limitations.** It is an equilibrium GL calculation rather than dissipative TDGL with a dynamic exterior field.

### Alstrøm et al. (2011)

**Citation.** T. S. Alstrøm, M. P. Sørensen, N. F. Pedersen, and S. Madsen, “Magnetic flux lines in complex geometry type-II superconductors studied by the time dependent Ginzburg–Landau equation,” *Acta Applicandae Mathematicae* 115, 63–74 (2011), [DOI: 10.1007/s10440-010-9580-8](https://doi.org/10.1007/s10440-010-9580-8), [repository copy](https://www2.imm.dtu.dk/pubdb/edoc/imm5949.pdf), [local PDF](../papers/kaper-et-al-2010-complex-geometry-tdgl.pdf).

This work uses COMSOL/nodal finite elements in two-dimensional complex sample geometries, temporal gauge, and implicit time integration to compare vortex configurations with experiment.

**Relevance.** It shows why unstructured meshes are attractive for experimentally fabricated shapes.

**Limitations.** The formulation imposes applied magnetic data and a gauge condition on the sample boundary rather than solving a surrounding vacuum. It is 2-D and does not address $H(\mathrm{curl})$ regularity or topological harmonic modes.

## Maxwell-compatible finite elements and exterior domains

### Nédélec (1980) and Monk (2003)

**Citations.** J.-C. Nédélec, “Mixed finite elements in $\mathbb R^3$,” *Numerische Mathematik* 35, 315–341 (1980), [DOI: 10.1007/BF01396415](https://doi.org/10.1007/BF01396415). P. Monk, *Finite Element Methods for Maxwell's Equations*, Oxford University Press (2003), [publisher page](https://global.oup.com/academic/product/finite-element-methods-for-maxwells-equations-9780198508885).

Nédélec elements place tangential degrees of freedom on edges and form the compatible $H^1\to H(\mathrm{curl})\to H(\mathrm{div})\to L^2$ sequence with nodal, face, and cell spaces. This structure represents curl fields, avoids spurious Maxwell modes, and provides the right continuity across material interfaces. Monk gives the comprehensive variational and approximation theory.

**Relevance.** These are the mathematical basis for custom tetrahedral edge-element assembly.

**Limitations.** They are Maxwell references, not coupled TDGL solvers; nonlinear complex order-parameter coupling, gauge fixing, and time integration must be added.

### Feischl and Tran (2016): FEM–BEM exterior coupling

**Citation.** M. Feischl and T. Tran, “The Eddy Current–LLG Equations—Part I: FEM–BEM Coupling,” [arXiv:1602.00744](https://arxiv.org/abs/1602.00744), [local PDF](../papers/feischl-tran-2016-fem-bem-eddy-current.pdf). The combined journal treatment is “The Eddy Current–LLG Equations: FEM–BEM Coupling and A Priori Error Estimates,” *SIAM Journal on Numerical Analysis* 55, 1786–1819 (2017), [DOI: 10.1137/16M1065161](https://doi.org/10.1137/16M1065161).

The paper couples interior three-dimensional Nédélec finite elements to boundary elements for an unbounded eddy-current exterior and proves convergence under realistic interface regularity.

**Relevance.** It is the closest transferable method for replacing a large artificial vacuum box in the MQS production solver. It also clarifies the correct electromagnetic trace variables.

**Limitations.** It treats eddy currents, not the nonlinear TDGL order parameter. Implementing singular boundary integrals and scalable dense/hierarchical algebra in MATLAB is a later phase.

### Du and Wu (1999): artificial boundaries for superconductivity

**Citation.** Q. Du and X. Wu, “A numerical method for the time-dependent Ginzburg–Landau equations in three dimensions,” *SIAM Journal on Numerical Analysis* 36, 1573–1597 (1999), [DOI: 10.1137/S0036142997330317](https://doi.org/10.1137/S0036142997330317).

This work emphasizes that, in three dimensions, physical electromagnetic variables for a finite sample belong to the whole space. It derives exact and approximate artificial boundary treatments.

**Relevance.** It directly supports the project constraint: an applied field imposed at the superconducting surface is not a faithful finite-sample self-field model. A surrounding exterior and a controlled far-field closure are required.

**Limitations.** The implementation context predates current edge-element and preconditioning practice, so the boundary machinery must be recast in a modern compatible discretization.

### Bao et al. (2022): adaptive full-wave Maxwell with DtN

**Citation.** G. Bao, M. Zhang, X. Jiang, P. Li, and X. Yuan, “An adaptive finite element DtN method for Maxwell's equations,” *East Asian Journal on Applied Mathematics* 13, 610–645 (2023), [DOI: 10.4208/eajam.2022-289](https://doi.org/10.4208/eajam.2022-289), [arXiv:2202.09203](https://arxiv.org/abs/2202.09203), [local PDF](../papers/bao-et-al-2022-adaptive-fem-dtn-maxwell.pdf).

The method uses curl-conforming elements, a spherical Dirichlet-to-Neumann map, and an a posteriori estimator that separates FE and DtN truncation errors for an unbounded radiation problem.

**Relevance.** It is a good template if a future full electromagnetic-wave branch needs a controlled nonreflecting boundary.

**Limitation/incompatibility.** A time-harmonic radiation DtN map is not the correct closure for the baseline MQS magnetostatic/eddy-current exterior. It must not be imported without changing the physical model.

### Chen, Guo, and Zou (2022): unfitted material interfaces

**Citation.** L. Chen, R. Guo, and J. Zou, “A family of immersed finite element spaces and applications to three-dimensional H(curl) interface problems,” [arXiv:2205.14127](https://arxiv.org/abs/2205.14127), [local PDF](../papers/chen-guo-zou-2022-immersed-hcurl.pdf).

The authors construct immersed $H(\mathrm{curl})$ spaces on unfitted tetrahedral meshes and prove optimal approximation for Maxwell interface problems.

**Relevance.** This is a possible later route for moving inclusions or geometry sweeps without remeshing.

**Limitations.** It is substantially more complex than fitted elements and does not solve TDGL. A boundary-fitted first implementation is lower risk.

## Evidence-based synthesis

The evidence leads to the following design decisions.

- **Physics baseline:** use magnetoquasistatic TDGL–Maxwell with normal conductivity and current continuity. Keep the displacement-current/full-wave Maxwell model as a documented extension, not an ambiguous optional term.
- **Exterior field:** solve the vector potential in a surrounding vacuum region. Put applied-field data on a remote outer boundary or represent source coils/currents. Demonstrate convergence as the vacuum padding grows. Add FEM–BEM later for an exact unbounded MQS exterior.
- **Spaces:** use complex nodal $H^1$ elements for $\psi$, Nédélec $H(\mathrm{curl})$ elements for $\mathbf A$, and nodal gauge/scalar-potential unknowns. Consider Raviart–Thomas/BDM magnetic flux or current variables when local conservation warrants a larger mixed system.
- **Gauge/topology:** enforce Coulomb gauge weakly with a multiplier and explicit nullspace treatment. Compute a cohomology/cycle basis in multiply connected domains; gauge fixing alone does not remove harmonic fields. Retain a temporal-gauge branch only as a validation path.
- **Time:** start with fully implicit backward Euler and Newton/Picard damping. Add BDF2 or linearized Crank–Nicolson after energy, gauge, and convergence tests pass.
- **Geometry/adaptivity:** use fitted tetrahedral meshes first, with edge orientation and exact-sequence incidence matrices treated as core infrastructure. Refine around vortex cores and penetration layers; unfitted interface elements are a later option.
- **Proximity:** a GL-coefficient continuation into a normal metal is a near-$T_c$, phenomenological proximity model. Quantitative low-temperature S–N proximity requires Usadel/BdG coupling and is outside the initial solver's claims.

These choices resolve the main incompatibility in the literature: efficient bulk TDGL solvers often prescribe $\mathbf A$ or $\mathbf B$, while rigorous 3-D Maxwell FEM papers omit superconducting dynamics. The proposed method couples the strengths of both without treating a sample-boundary applied field as a self-consistent exterior solution.
