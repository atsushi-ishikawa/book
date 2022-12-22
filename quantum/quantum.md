# Basics in Theoretical and Computational Methods
* The fundamental players in the catalytic science is atoms and molecules, as the catalysis is the object that is promoting chemical reactions.
* Chemical reactions involve the formation and the breaking of the chemical bond, and because of this electrons consisting the chemical bond is also a matter of subject.
* Since the behavior of microparticles such as electrons or nucleus is described the quantum physics, we should start from such area. In other words, catalytic science (and also surface science) is a research field at the border between chemistry and physics. Most of the theoretical methods have been derived either from quantum physics or from condensed matter physics.

## 1. Electronic Structure Theory
* The behavior of the particles is governed by the equation of motion, and its classical mechanical version is known as the Newton's law.
* Since we are interested with microscopic particles like electrons or atoms, the proper description should be given by the laws of quantum mechanics.
* For this reason, we need to consider the **Schrödinger equation**, which is a quantum-mechanical version of the equation of motion.
* If the solutions of the Schrodinger equations are generated without reference to experimental data, the methods are usually called "ab initio" (latin: "from the beginning") or "first principle".

* The relativistic version of the Schrödinger equation, known as **Dirac equation**, provides more accurate description of microparticles. However, for most of chemistry problems, non-relativistic Schrödinger equation under the Born-Oppenheimer approximation is sufficient. This is mainly because the relativistic effect is large for core electrons but its effect to valence electrons are smaller. The problems for core electrons and valence electrons can be separated to each other by introducing pseudopotentials. Thus researchers often concentrate on the Schrodinger equation for valence electrons. The details of the pseudopotentials will be discussed in the later chapters.

* From above reasons, we consider the Schrödinger equation under following approximations throughout the book;
	1. Time-independent.
  2. When only valence electrons are considered, the relativistic effect is small.
  3. electron and nuclear motions are decoupled.

* It is convenient to use the Dirac's "bra-ket notation" for wave functions and multi-dimensional integrals in electronic structure theory in order to simplify the notation. The equivalences are defined as
$$
\begin{align*}
\ket{\Psi} \equiv \Psi;&\ \bra{\Psi} \equiv \Psi^{*} \\
\int{\Psi^{*}\Psi d{\bf r}} &= \braket{\Psi|\Psi} \\
\int{\Psi^{*}\hat{H}\Psi d{\bf r}} &= \braket{\Psi|\hat{H}|\Psi}
\end{align*}
$$
The ket $\ket{m}$ denotes a wave function with quantum number $m$ standing to the right of the operator, while the bra $\bra{n}$ denotes a complex conjugate wave function with quantum number $n$ standing to the left of the operator. The combined braket denotes that the whole expression should be integrated over all coordinates.

#### Born-Oppenheimer Approximation
* All physical and chemical properties of any system can be derived from its Hamiltonian.
* In this section, we derive the quantum-mechanical Hamiltonian under the Born-Oppenheimer approximation as this is the basics for the most of the chemistry problems.
* Except for hydrogen and helium, atoms have a mass that is $10^4$ to $10^5$ times heavier than the mass of an electron. Consequently, at the same kinetic energy, electrons moves $10^2$ to $10^3$ times faster than the nuclei. Hence one supposes that the electrons follow the motion of the nuclei almost instantaneously.
* Since the motion of nucleus is much less than that of electrons, we can consider the nucleus to be fixed while electrons are moving. Because of this, we can separate the electrons' degree of freedom by the nucleus' degree of freedom.
* The Schrödinger equation for the electrons for a fixed nuclear coordinate is then
$$
H_{el}(\left\{ {\bf R} \right\}) \Psi(r, \left\{ {\bf R} \right\}) = E_{el}(\left\{ {\bf R} \right\}) \Psi(r, \left\{ {\bf R} \right\})
$$
* Again, the nuclear coordinates $\left\{ {\bf R} \right\}$ are not meant to be variables but parameters i.e. the R-values does not change when solving the above equation.
* From the nucleus side, the electrons' motion is separated out and the motion of electron is not explicitly coupled to the nucleus motion. In other words, for nucleus, the electrons act as the "field" and not the freely-moving particles.
* By using the Born-Oppenheimer approximation, one can make the direct functional relationship (or "mapping") between the nuclear coordinates $\{{\bf R}\}$ and the electronic energy $E_{el}$. It is quite informative to visualize $E_{el}(\left\{ {\bf R} \right\})$ as the function of {R}, which is often called **potential energy surface** (Fig.X).
* Since the geometry or structures of molecules, bulk materials, surfaces are all described by their nuclear coordinates, this mapping provides the relationship between the molecular geometry and the energy. The energy is related with the stability of the molecules via chemical bond formation/breaking, thus it is a fundamental property for discussing the reactivity of the molecules. For this reason, the property of catalysts can be finally reduced to the atomic or molecular quantum physics within the Born-Oppenheimer approximation.
* For this reason, all the discussions in the following book are based on the Born-Oppenheimer approximation.

### 1-1. Wavefunction based theory
* Quantum chemical methods usually describe the electrons by a localized basis set (or basis function) derived from atomic orbitals.
* Historically, two localized basis set is commonly used; the Slater type functions $\exp(-\alpha r_{iA})$ or the Gaussian type functions $\exp(-\alpha r_{iA}^2)$. Here, $\alpha$ is the exponent of the function and $r_{iA}$ is the distance between the electron $i$ and the nucleus $A$. The Gaussian type function is more often used because they allow the analytic evaluation of the matrix elements necessary to perform an electronic structure calculation. The most popular electronic structure code in this line is the GAUSSIAN program. However, some code (such as Amsterdam Density Functional package) uses the Slater type functions in their calculation.

* When the N-body wave function is expressed as the form of the product of the one-electron functions, it becomes like
$$
Hartree-product
$$
This is called the *Hartree-product*.

#### Hartree-Fock method
* The Hartree product obeys the Pauli principle only to some extent, by populating each electronic state once. However, it does not take into account the anti-symmetry of the wave function.
* The Pauli principle requires that the sign of $\Psi$ changes when two electrons are exchanged.
* The simplest ansatz obeying the anti-symmetry requirement is to replace the product wave function by a single **Slater determinant**.
* The Slater determinant is constructed from the single-particle wave functions by
$$
\Psi({\bf r_1, r_2}\cdots {\bf r_N}) 
= \frac{1}{\sqrt{N!}}
\begin{vmatrix}
\psi_1({\bf r_1}) & \cdots & \psi_1({\bf r_N}) \\
\vdots            & \ddots & \vdots \\
\psi_N({\bf r_1}) & \cdots & \psi_N({\bf r_N}) \\
\end{vmatrix}
$$
* Then the expectation value of the total energy becomes
$$
\begin{split}
\braket{\Psi|H|\Psi} 
&= \sum_{i=1}^N\int d^3r \psi_i^*({\bf r}) \left(-\frac{\hbar^2}{2m}\nabla^2 + v_{ext}({\bf r}) \right) \psi_i({\bf r}) \\
+& \frac{1}{2}\sum_{i,j=1}^N\int d^3r d^3r' \frac{e^2}{|r-r'|}
\psi_i^*({\bf r})\psi_i({\bf r})\psi_j^*({\bf r'})\psi_j({\bf r'}) \\
-& \frac{1}{2}\sum_{i,j=1}^N\int d^3r d^3r' \frac{e^2}{|r-r'|}
\psi_i^*({\bf r})\psi_i({\bf r'})\psi_j^*({\bf r'})\psi_j({\bf r}) \\
+& V_{nuc-nuc}
\end{split}
$$
* The last term is called **exchange term**.
* If we minimize the above expression with respect to the $\psi_i^*$ under the constraint of normalization, we have **Hartree-Fock equation**.
$$
\begin{split}
&\left\{ -\frac{1}{2}\nabla^2 + v_{ext}({\bf r}) + v_H({\bf r}) \right\} \psi_i({\bf r}) \\
& -\sum\int\frac{e^2}{|r-r'|}\psi_j^*({\bf r'})\psi_i({\bf r'})\psi_j({\bf r}) = \epsilon_i \psi_i({\bf r})
\end{split}
$$

#### Electron correlation
* The Hartree-Fock method generates solutions to the Schrodinger equation where the real electron-electron interaction is replaced by an average interaction. In a sufficiently large basis set, the HF wave function accounts for ~99% of the total energy, but the remaining 1% is often very important for describing chemical phenomena such as chemical bond formation or breaking.
* The difference in energy between the HF and the lowest possible energy in the given basis set is called the *electron correlation energy*.
* Physically, it corresponds to the motion of the electrons being correlated, i.e. on the average they are further apart than described by the HF wave function.
* The HF method determines the energetically best one-determinant trial wave function (within the given basis set). It is therefore clear that, in order to improve on the HF results, the starting point must be a trial wave function that contains more than one Slater determinant.
* As the HF solution usually gives ~99% of the correct answer, electron correlation methods normally use the HF wave function as a starting point for improvements.

##### Excited determinants
* The starting point is usually an restricted Hartree-Fock (RHF) calculation where a solution of the Roothaan-Hall equations for a system with $N_{elec}$ electrons and $M_{basis}$ basis functions will yield $1/2 N_{elec}$ occupied MOs and $M_{basis} - 1/2 N_{elec}$ unoccupied (virtual) MOs. Except for a minimal basis, there will always be more virtual than occupied MOs.
* A whole series of determinants may be generated by replacing MOs that are occupied in the HF determinant by MOs that are unoccupied. There can be denoted according to how many occupied HF MOs have been replaced by unoccupied MOs; Slater determinants that are singly, doubly, triply, quadruply etc., excited relative to the HF determinant, up to a maximum of $N_{elec}$ excited electrons. These determinants are often referred as singles (S), doubles (D), triples (T), quadruples (Q), etc.

##### Configuration interaction
* The trial wave function is written as a linear combination of determinants with the expansion coefficient determined by requiring that the energy should be minimum. This procedure is known as **configuration interaction (CI)**. The MOs used for building the excited Slater determinants taken from a Hartree-Fock calculation, and its MO coefficients are held fixed. Subscripts S, D, etc. indicate determinants that are singly, double, etc. excited relative to the Hartree-Fock configuration.
$$
\Psi_{CI} = a_0\Psi_{HF} + \sum_{S}a_S\Psi_S + \sum_{D}a_D\Psi_D + \cdots = \sum_{i=0}a_i\Psi_i
$$
* The coefficients $a_i$ is determined by solving the CI secular equations. Solving the secular equation is equivalent to diagonalization.

##### MCSCF
* The **multi-configurational self-consistent field (MCSCF)** method can be considered as a CI where not only are the coefficients in front of the determinants but the MOs used for constructing the determinants are also optimized.
* The MCSCF optimization is iterative like the SCF procedure.
* The major problem with the MCSCF method is selecting which electron configurations are necessary for the property of interest. One of the most popular approaches is the *complete active space self-consistent field (CASSCF)* method.
* Here the selection of configurations is done by partitioning the MOs into active and inactive spaces (Fig.X). The active MOs will typically be some of the highest occupied and some of the lowest unoccupied MOs from a RHF calculation.
* The inactive MOs have either 2 or 0 electrons, i.e. always either doubly occupied or empty.
* Within the active MOs a full CI is performed and all the proper symmetry-adapted configurations are included in the MCSCF optimization.
* Which MOs to include in the active space must be decided manually, by considering the problem at hand and the computational expense.

##### Moller-Plesset perturbation theory
* The idea in the perturbation methods is that the problem at hand only differs slightly from a problem that has already been solved (exactly or approximately). The solution to the given problem should therefore in some sense be close to the solution to the already known system. This is described mathematically by defining a Hamiltonian operator that consists of two parts, a reference and a perturbation.
* Since the solutions to the unperturbed Schrödinger equation generate a complete set of functions, the unknown first-order correction to the wave function can be expanded in these functions. This is known as **Rayleigh-Schrodinger perturbation theory**.
* In order to apply the perturbation theory to the calculation of correlation energy, the unperturbed Hamiltonian must be selected. The most common choice is to take this as sum over Fock operators; this method is called **Moller-Plesset perturbation theory**. In this method, the zeroth-order wave function is the Hartree-Fock determinant, and the zeroth-order energy is just a sum of MO energies.
* By using the sum of the Fock operator as the zeroth-order Hamiltonian, the first-order energy correction becomes the (half of) electron-electron repulsion energy, and the summation of the zeroth- and first-order energies gives the Hartree-Fock energy. Therefore the second-order correction to the energy is the first contribution to the electron correlation energy. This involves a sum over doubly excited determinants, which is generated by promoting two electrons from occupied orbital $i$ and $j$ to virtual orbitals $a$ and $b$.
$$
W_2 = \sum_{i<j}^{occ}\sum_{a<b}^{vir}
			\frac{ \braket{\Psi_0|H'|\Psi_{ij}^{ab}}
		           \braket{\Psi_{ij}^{ab}|H'|\Psi_0} }
				   {E_0-E_{ij}^{ab} }
$$
* The matrix elements between the HF and a double excited state are given by two-electron integrals over MOs. The difference in total energy between two Slater determinants becomes a difference in MO energies. Thus the explicit formula for the second-order Moller-Plesset correction is
$$
E(MP2) = \sum_{i<j}^{occ}\sum_{a<b}^{vir}
				 \frac{ \braket{\phi_i\phi_j|\phi_a\phi_b}
					  - \braket{\phi_i\phi_j|\phi_b\phi_a} }
{\varepsilon_i + \varepsilon_j - \varepsilon_a - \varepsilon_b}
$$
* Once the two-electron integrals over MOs are available, the second-order energy correction can be calculated as a sum over such integrals.
* MP2 typically accounts for 80-90% of the correlation energy, and it is the most economical method for including electron correlation. In practical calculations, the MP2 energy for systems with a few hundred basis functions can be calculated at a cost similar to or less than what is required for calculating the HF energy. It is also possible to include the higher order perturbation correction; the n-th order Moller-Plesset perturbation theory is often denoted as *MPn* method, and MP2 to MP4 is used in practical calculations.
* One should note that when the reference wave function contains substantial multi-reference character (e.g. at the bond-dissociation limit), a perturbation expansion based on a single determination is no longer valid.
* There are now several different implementation for the multi-reference many-body perturbation theory; the most popular one is to use the CASSCF as a reference, denoted as the *complete active-space second-order perturbation theory (CASPT2)*.

##### Coupled cluster method
* Perturbation methods add all types of corrections (S, D, T, etc.) to the reference wave function to a given order (2, 3, etc.). The idea in **coupled cluster** methods is to include all excited determinants of a given type.
* Define an excitation operator as
$$
T=T_1 + T_2 + T_3 + \cdots + T_{N_{elec}}
$$
* The $T_i$ operator acting on the HF reference wave function $\Psi_o$ generates all i-th excited Slater determinants.
$$
T_1\Psi_0 = \sum_i^{occ}\sum_j^{vir}t_i^a\Psi_i^a \\
T_2\Psi_0 = \sum_{i<j}^{occ}\sum_{a<b}^{vir}t_{ij}^{ab}\Psi_{ij}^{ab}
$$
* In coupled cluster theory, it is customary to use the terms *amplitudes* for the expansion coefficient $t$.
* The coupled cluster wave function is defined as
$$
\Psi_{CC} = \exp(T)\Psi_0 \\
\exp(T) = 1 + T + \frac{1}{2}T^2 + \frac{1}{6}T^3 + \cdots = \sum_{k=0}^{\infty}\frac{1}{k!}T^k
$$
* Usually, the $T$ operator is truncated some of the terms. For example, using $T=T_1+T_2$ gives the coupled cluster singles and doubles (CCSD) method.
* Unfortunately, a straightforward application of the variational principle to the coupled cluster wave function i.e. variational coupled cluster approach leads a large number of excitation operators, thus it can be solved only for small systems. Instead, projecting the Schrodinger equation for the coupled cluster wave function onto the reference (such as Hartree-Fock) wave function $\Psi_0$ is often done (non-variational coupled cluster).
$$
\braket{\Psi_0|H\exp(T)|\Psi_0} = E_{CC}\braket{\Psi_0|\exp(T)\Psi_0}
$$
* Expanding out the exponential and using the fact that the Hamiltonian operator contains only one- and two-electron operators, we get
$$
E_{CC} = \braket{\Psi_0|H(1 + T_1 + T_2 + 1/2 T_1^2)|\Psi_0} \\
E_{CC} = \braket{\Psi_0|H|\Psi_0} + \braket{\Psi_0|H|T_1\Psi_0} + \braket{\Psi_0|H|T_2\Psi_0} + \frac{1}{2}\braket{\Psi|H|T_1^2\Psi_0} \\
E_{CC} = E_0 + \sum_i^{occ}\sum_a^{vir}t_i^a\braket{\Psi_0|H|\Psi_i^a} + \sum_i^{occ}\sum_a^{vir}(t_{ij}^{ab} + t_i^a t_j^b - t_i^b t_j^a)\braket{\Psi_0|H|\Psi_{ij}^{ab}}
$$
* Then the coupled cluster energy is given by
$$
E_{CC} = E_0 + \sum_{i<j}^{occ}\sum_{a<b}^{vir}(t_{ij}^{ab}+t_i^a t_j^b - t_i^b t_j^a)(\braket{\phi_i\phi_j|\phi_a\phi_b}-\braket{\phi_i\phi_j|\phi_b\phi_a})
$$
where the equations for amplitude should be solved in an iterative manner since the parameters enter in an non-linear fashion.
* If all cluster operators up to $T_N$ are included in $T$, all possible excited determinants are generated and the coupled cluster wave function is equivalent to full CI. This is, as already stated, impossible for all but the smallest systems.
* The cluster operator must therefore be truncated at some excitation level. Using $T = T_1 + T_2$ gives the *coupled cluster with singles and doubles (CCSD)* model, and $T = T_1 + T_2 + T_3$ gives the *coupled cluster with singles, doubles, and triples (CCSDT)*. However, CCSDT is only applicable to small systems.
* Instead of including $T_3$ in the cluster operator $T$, the triples contribution may be evaluated by perturbation theory and added to the CCSD results. The most commonly such approach is the *CCSD(T)* method. In this case, the triples contribution is calculated from the formula given by MP4, but using the CCSD amplitudes instead of the perturbation coefficients for the wave function correction.
* The CCSD(T) method with reasonable size of the basis set (at least better than double-zeta plus polarization function) is called *golden standard* for quantum chemistry calculations, because usually one can expect the *chemical accuracy* (at kcal/mol level accuracy) with this level.

### 1-2. Density Functional Theory
<!-- Jensen start -->
* A *function* is a prescription for producing a number from a set of *variables*. A *functional* is a prescription for producing a number from a *function*, which in turn depends on variables. A wave function and the electron density are thus *functions*, while the energy depending on a wave function or an electron density is a *functional*. In this book, we will denote a function depending on a set of variables ${\bf x}$ with $f({\bf x})$, while a functional depending on a function $f$ is denoted as $F[f]$.

#### Orbital-free DFT
* The electronic energy functional of electron density $\rho$ can be divided into three parts; kinetic energy functional $T[\rho]$, the nuclei-electron attraction energy functional $E_{ne}[\rho]$, and the electron-electron repulsion energy functional $E_{ee}[\rho]$.
* With reference to the Hartree-Fock theory, $E_{ee}[\rho]$ may be divided into Coulomb and exchange parts, $J[\rho]$ and $K[\rho]$.
* Among these energy functionals, $E_{ne}[\rho]$ and $J[\rho]$ can be interpreted by the classical electrodynamics, 
$$
E_{ne}[\rho] = -\sum_i^{N_{nuc}}\int\frac{Z_a(R_a)\rho({\bf r})}{|{\bf R}_a - {\bf r}|}d{\bf r} \\
J[\rho] = \frac{1}{2}\int\int\frac{\rho({\bf r})\rho({\bf r'})}{|{\bf r} - {\bf r'}|}d{\bf r}d{\bf r'}
$$
* Here, the factor of 1/2 in $J[\rho]$ allows the integration to be done over all space fo both ${\bf r}$ and ${\bf r'}$ variables. Unlike these energy components, exchange part $K[\rho]$ can only be interpreted by the quantum mechanics.
* Early attempts of deducing functionals for the kinetic and exchange energies considered a uniform electrons gas where it may be shown that $T[\rho]$ and $K[\rho]$ are given by
$$
T_{TF}[\rho] = C_F \int \rho^{5/3}({\bf r})d{\bf r} \\
K_D[\rho] = -C_X \int \rho^{4/3}({\bf r})d{\bf r} \\
C_F = \frac{3}{10}(3\pi^2)^{2/3}, \ C_X = \frac{3}{4}\left( \frac{3}{\pi} \right)^{1/3}
$$
* The energy functional $E_{TF}[\rho] = T_{TF}[\rho] + E_{ne}[\rho] + J[\rho]$ is known as *Thomas-Fermi theory*, while including $K_D[\rho]$ exchange part constitutes the *Thomas-Fermi-Dirac model*.
* Since $T[\rho]$ and $K[\rho]$ functionals are depending directly on the electron density, these methods are called *orbital-free DFT*, as opposed to the Kohn-Sham theory discussed in the next section.
* Unfortunately, the accuracy of the orbital-free DFT is too low to be of general use.

#### Kohn-Sham theory
* The success of modern DFT method is based on the suggestion by Kohn and Sham in 1965 that the electron kinetic energy should be calculated from an auxiliary set of orbitals used for representing the electron density.
* The main drawback in the orbital-free DFT is the poor representation of the kinetic energy, and the idea in the Kohn-Sham formalism is to split the kinetic energy functional into two parts, one which can be calculated exactly, and a small correction term.
* Assume for the moment a Hamiltonian operator of the form with $0 \le \lambda \le 1$.
$$
H_{\lambda} = {\bf T} + {\bf V}_{ext}(\lambda) + \lambda {\bf V}_{ee}
$$
* The external potential operator ${\bf V}_{ext}$ is equal to ${\bf V}_{ee}$ for $\lambda = 1$, but for intermediate $\lambda$ value it is assumed that ${\bf V}_{ext}(\lambda)$ is adjusted such that the same density is obtained for $\lambda = 1$ (the real system), for $\lambda = 0$ (a hypothetical system with non-interaction electrons) and for all intermediate $\lambda$ values.
* For the $\lambda = 0$ case, the electrons are non-interacting, and the exact solution to the Schrodinger equation is given as a Slater determinant composed of molecular orbitals $\phi_i$. Then the exact kinetic energy functional is
$$
T_{KS} = \sum_i^{N_{el}}\braket{\phi_i|-\frac{1}{2}\nabla^2|\phi_i}
$$
* The $\lambda = 1$ corresponds to interacting electrons, and Eq.X is therefore only an approximation to the real kinetic energy.
* The key to Kohn-Sham theory is to calculated the kinetic energy under the assumption of non-interacting electrons (in the sense that the HF orbitals in the wave function theory).
* In reality, the electrons are interacting, and Eq.X does not provide the total kinetic energy. However, just as HF theory provides ~99% of the correct answer, the difference between the exact kinetic energy and that calculated by assuming non-interacting orbitals is small.
* The remaining kinetic energy is absorbed into an exchange-correlation term, and a general DFT energy expression can be
$$
E_{DFT}[\rho] = T_{KS}[\rho] + E_{ne}[\rho] + J[\rho] + E_{XC}[\rho]
$$
* By equating $E_{DFT}$ to the exact energy, the above expression actually defines $E_{XC}$ i.e. it is the part that remains after subtraction of the non-interacting kinetic energy, and $E_{ee}$ and $J$ potential energy terms;
$$
E_{XC}[\rho] = (T[\rho] - T_{KS}[\rho]) + (E_{ee}[\rho] - J[\rho])
$$
* The task in developing orbital-free DFT is to derive approximations to the kinetic, exchange, and correlation energy functionals, while the corresponding task in Kohn-Sham theory is to derive approximations to the exchange-correlation energy functional only.
* Since the exchange-correlation energy is roughly a factor of 10 smaller than the kinetic energy, the Kohn-Sham theory is much less sensitive to inaccuracies in the functionals than the orbital-free DFT.
* The division of the electron kinetic energy into two parts, with the major contribution being equivalent to the Hartree-Fock energy, can be justified as follows.

#### Exchange-correlation hole
* Electrons avoid each other owing to their electric charges, and the energy associated with this repulsion is given classically by the Coulomb's law.
* Quantum mechanically, however, this repulsion must be modified to take into account that electrons have spins of 1/2. The Pauli principle states that two Fermions (particles with half-integer spin) cannot occupy the same spatial position, or equivalently, that the total wave function must be antisymmetric upon interchange of any two particles.
* These quantitative considerations can be put into quantitative terms by probability holes, or specifically exchange and correlation holes.
* To this aim, here we define the *one-electron density* $\rho_1$ and the *electron-pair density* $\rho_2$. The former is
$$
\rho_1({\bf r}_1) = N_{el}\int
\Psi^{*}({\bf r}_1; {\bf r}_2, \cdots, {\bf r}_{N_{el}})
\Psi({\bf r}_1; {\bf r}_2, \cdots, {\bf r}_{N_{el}})d{\bf r}_2 \cdots d{\bf r}_{N_{el}}
$$
* Note that the integration is done for ${\bf r}_2$ to ${\bf r}_{N_{el}}$, thus $\rho_1$ is the function of only ${\bf r}_1$. The integral is the probability of finding an electron a position ${\bf r}_1$, and the $N_{el}$ prefactor ensures that the integral value equals the number of electrons.
* The electron-pair density is
$$
\rho_2({\bf r}_1, {\bf r}_2) = N_{el}(N_{el}-1)\int
\Psi^{*}({\bf r}_1, {\bf r}_2; {\bf r}_3, \cdots, {\bf r}_{N_{el}})
\Psi({\bf r}_1, {\bf r}_2; {\bf r}_3, \cdots, {\bf r}_{N_{el}})d{\bf r}_3 \cdots d{\bf r}_{N_{el}}
$$
* Here, the integration is done for ${\bf r}_3$ to ${\bf r}_{N_{el}}$ thus $\rho_2$ has two set of variables, i.e. ${\bf r}_1$ and ${\bf r}_2$. $\rho_2$ means the probability of finding an electron at position ${\bf r}_1$ and another electron at position ${\bf r}_2$, and the $N_{el}(N_{el}-1)$ prefactor ensures that $\rho_2$ integrates to the number of electron pairs.
* Let us go back to the exchange-correlation hole issue. If electrons did not have charge or spin, the probability of finding an electron at a given position would be independent of the position of a second electron, and $\rho_2$ would be given as a simple product of two $\rho_1$, with a proper normalization factor.
$$
\rho_2^{indep}({\bf r}_1, {\bf r}_2) = \frac{N_{el}-1}{N_{el}}\rho_1({\bf r}_2)\rho_1({\bf r}_2)
= \left( 1-\frac{1}{N_{el}}\right)\rho_1({\bf r}_2)\rho_1({\bf r}_2)
$$
* Since electrons have both charge and spin, however, there is a reduced probability of finding an electron near another electron. We can write this formally in terms of a conditional probability factor $h_{xc}({\bf r}_1, {\bf r}_2)$ as
$$
\rho_2({\bf r}_1, {\bf r}_2) = \rho_1({\bf r}_1)\rho_1({\bf r}_2) + \rho_1({\bf r}_1)h_{xc}({\bf r}_1, {\bf r}_2)
$$
* The reduced probability is called the *exchange-correlation hole*, and can be written in terms of $\rho_1$ and $\rho_2$ as
$$
h_{xc}({\bf r}_1, {\bf r}_2) = \frac{\rho_2({\bf r}_1, {\bf r}_2)}{\rho_1({\bf r}_1)} - \rho_1({\bf r}_2)
$$
* The exchange-correlation hole represents the reduced probability of finding electron 2 at position ${\bf r}_2$ given that electron 1 is located at ${\bf r}_1$ (see Fig 6.1 of Jensen).
* The exchange part of $h_{xc}$ is called the Fermi hole, while the correlation part is the Coulomb hole.
* Since exchange only occurs between electrons of the same spin, the total hole can also be written in terms of individual spin contributions.
$$
h_{xc} = h_x + h_c \\
h_x = h_x^{\alpha\alpha} + h_x^{\beta\beta} \\
h_c = h_c^{\alpha\alpha} + h_c^{\beta\beta} + h_c^{\alpha\beta}
$$
* From the definitions of $\rho_1$ and $\rho_2$, it follows that the integral of $h_{xc}$ over ${\bf r}_2$ gives -1;
$$
\int h_{xc}({\bf r}_1, {\bf r}_2)d{\bf r}_2 = \int\frac{\rho_2({\bf r}_1, {\bf r}_2)}{\rho_1({\bf r}_1)}d{\bf r}_1 d{\bf r}_2 - \int \rho_1({\bf r}_2)d{\bf r}_2
= \frac{N_{el}(N_{el}-1)}{N_{el}} - N_{el} = -1
$$
* This means that the Fermi hole integrates to -1 (and always $h_x$ takes negative value). The Fermi hole describes a static reduction in the probability function corresponding one-electron. On the other hand, the Coulomb hole takes both positive and negative values, and it integrates to 0. This means that the Coulomb hole reduced the probability of finding an electron near the reference electron but increases the probability of finding it far from the reference electron.
* The exchange-correlation hle is necessary condition that the exact exchange-correlation functional should satisfy. Therefore this property is an important guideline in deriving an approximate exchange-correlation functional, as we will see them in the next sections.

#### Exchange-correlation functionals
* The difference between various DFT methods is the choice of functional form in the exchange-correlation (XC) energy.
* It can be proven that the XC potential is a unique functional, valid for all systems.
* However, an explicit functional form of this potential has been elusive, except for special cases such as an uniform electron gas.
* XC functionals have a mathematical form containing parameters. There are two main philosophies for assigning values to these parameters, either by requiring the functional to fulfill the necessary condition that the exact functional has to satisfy, or by fitting the parameters to experimental data. Often, a combination of these two approaches is needed in practice.
* It is possible to separate the XC energy $E_{xc}$ into two parts; a pure exchange energy $E_x$ and a correlation energy $E_c$, but the current trend is to construct the two parts in a combined fashion.

#### LDA
* In the *Local Density Approximation (LDA)*, it is assumed that the electron density can be treated locally as an uniform electron gas. The exchange energy $E_x$ for a uniform electron gas is given by the Dirac formula,
$$
E_x^{LDA}[\rho] = -C_x\int\rho^{4/3}({\bf r})d{\bf r} \\
\epsilon_x^{LDA} = -C_x \rho^{1/3}({\bf r})
$$
* In the more general case, where $\alpha$ and $\beta$ densities are not equal, the *Local Spin Density Approximation (LSDA)* is used;
$$
E_x^{LSDA}[\rho] = -2^{1/3}C_x\int(\rho_{\alpha}^{4/3}({\bf r}) + \rho_{\beta}^{4/3}({\bf r}))d{\bf r}
$$
* The analytical form for the correlation energy of an uniform electron gas, which is purely dynamical correlation, has been derived in the high and low density limits. For intermediate densities, the correlation energy has been determined to a high precision by quantum Monte Carlo calculation; the resulting functional is known as VWN functional, named after Vosko, Wilk, and Nasair.

#### GGA
* In the *Generalized Gradient Approximation (GGA)* methods, the first derivative of the density is included as a variable, and in addition it is required that the Fermi and Coulomb holes integrate to the required values of -1 and 0, respectively.
* One of the earliest GGA functional was proposed by A. D. Becke as a correction to the LSDA exchange energy. This is called B88 (Becke, 1988) exchange functional.
* A popular correlation functional at GGA level is LYP functional, proposed by Lee, Yang, and Parr.
* J. P. Perdew and co-workers have proposed several XC functionals based on removing spurious oscillations in the Taylor-like expansion to the first-order, and also ensuring that the exchange and correlation holes integrate to -1 and 0. The resulting functional are PW86 (Perdew-Wang 1986), PW91 (Perdew-Wang 1991), and PBE (Perdew-Burke-Ernzerhof).
* The parameters in the PBE functional are non-empirical, that is they are  not obtained by fitting to experimental data, but derived from the basic requirements for the XC functional.
* The PBE functional has been slightly modified to improve the performance for periodic systems (RPBE), while this modification actually destroys the hole condition for the exchange energy.

#### meta-GGA
* The logical extension of GGA method is to allow the EX functional to include the variables of higher order derivatives of the electron density, for example the Laplacian ($\nabla^2 \rho$) as the second-order derivative term.
* The same information is actually gained by introducing the orbital kinetic energy density $\tau$, which is
$$
\tau({\bf r}) = \frac{1}{2}\sum_i^{occ}|\nabla \phi_i({\bf r})|^2
$$
* and this approach is more often used.
* Including of either the Laplacian or the orbital kinetic energy density $\tau$ as a variable leads to so-called *meta-GGA functionals*.

#### Dispersion-corrected methods
* One of the serious shortcomings of standard DFT methods is the inability to describe the dispersion forces (part of the van der Waals-type interactions).
* Many functionals provide a purely repulsive interaction between rare gas atoms, while others describe a weak stabilization interaction, but fail to have the correct $R^{-6}$ long-range distance behavior.
* Although the dispersion is a short-ranged weak interaction, it is cumulative, and therefore becomes increasingly important as the system gets larger.
* S. Grimme has proposed to include dispersion by additive empirical terms.
* The parametrizatino in the earliest models simply included an $R^{-6}$ energy term for each atom pair, with an atom-dependent $C_6$ parameter. The method has been refined by including higher-order terms ($R^{-8}, $R^{-10}$) to better describe the medium-range dispersino and making the parameters depend on the atomic environment in terms of the number of directly bonded atoms.
* Adding such attractive energy terms has the potential problem of divergence for short interatomic distances and is conseqently often used in connection with a dampling function:
$$
\Delta E_{\rm disp} = -\sum_{n=6(8,10)}s_n\sum_{AB}^{N_{atoms}}\frac{C_n^{AB}}{R_{AB}}f_{\rm damp}(R_{AB})
$$
* A complication is that the parametrization of the dispersion correctino depends on the underlying XC functional, as different functinoals via their parametrization may include some of the short-range interaction, and this can be taken into account by a functional-dependent scaling factor $s_n$ in the above equation.
* Such dispersion corrected methods are denoted with a D/D2/D3 after the DFT acronym, for example PBE-D3.

* There is an another approach to the dispersion correction, that is making it to be directly dependent on the actual electron density. This is done by writing it as a six-dimensional integral over electron densities with an appropriate *dispersion kernel* $\Phi$ (the factor 1/2 corrects for double counting):
$$
\Delta E_{\rm disp} = \frac{1}{2}\int\rho({\bf r})\Phi({\bf r}, {\bf r}')\rho({\bf r}')d{\bf r}d{\bf r}'
$$
* Such dispersion methods are often denoted as *non-local*, as they depend on electron densities that can be far apart. The *van Vohris-Vydorv (VV10) kernel* has the form shown in
$$
\Phi({\bf r},{\bf r}') = -\frac{3}{2}\left[g({\bf r})g({\bf r}')\left\{g({\bf r})+g({\bf r}')\right\}\right]^{-1}
$$
* Here the $g({\bf r})$ or $g({\bf r}')$ functions are defined below with the $C$ and $b$ parameters chosen to provide the correct asymptotic $C_6$ coefficient and controlling the short-range dampling, respectively:
$$
g({\bf r}) = |{\bf r}-{\bf r}'|^2\sqrt{C\left(\frac{\nabla\rho({\bf r})}{\rho({\bf r})}\right) + \frac{4\pi}{3}\rho({\bf r})} + b\frac{3\pi}{2}\left(\frac{\rho({\bf r})}{9\pi}\right)^{1/6}
$$

#### Computational issue
* Once an XC functional has been selected, the computational problem of the Kohn-Sham DFT is very similar to the Hartree-Fock theory; determine a set of orthogonal orbitals that minimizes the energy.
* The orbital orthogonality constraint may be enforced by the Lagrange method, where Lagrangian is defined as
$$
L[\rho] = E_{DFT}[\rho] - \sum_{ij}^{N_{orb}}\lambda_{ij}\left(\braket{\phi_i|\phi_j}-\delta_{ij}\right)
$$
* Requiring the variation of $L$ to vanish provides a set of equations involving an effective one-electron operator $h_{KS}$, as
$$
h_{KS}\phi_i = \sum_j^{N_{orb}}\lambda_{ij}\phi_j \\
h_{KS} = -\frac{1}{2}\nabla^2 + v_{eff} \\
v_{eff} = v_{ne}({\bf r}) + \int\frac{\rho({\bf r}')}{|{\bf r} - {\bf r}'|}d{\bf r}' + v_{xc}({\bf r})
$$
* The effective potential $v_{eff}$ contain the nuclear-electron contribution $v_{ne}$, the Coulomb repulsion, and the XC potential $v_{xc}$. This term is given as the functional derivative of the XC energy with respect to the density, as
$$
\begin{align*}
v_{xc}({\bf r}) &= \frac{\delta E_{xc}[\rho]}{\delta \rho({\bf r})} \\
&= \epsilon_{xc}[\rho] + \int \rho({\bf r}')\frac{\delta \epsilon_{xc}({\bf r}')}{\delta \rho({\bf r}')}d{\bf r}'
\end{align*}
$$
* By making the matrix of the Lagrange multiplier to be a diagonal matrix, one obtains the canonical Kohn-Sham equation
$$
h_{KS}\phi_i = \epsilon_i \phi_i
$$
* Employing an expansion of the Kohn-Sham orbitals $\phi_i$ by atomic basis set gives
$$
\phi_i = \sum_{\alpha}^{M_{basis}}C_{\alpha i}\chi_{\alpha}
$$
* The variational procedure again leads to a matrix equation in the atomic orbital basis
$$
{\bf h}_{KS}{\bf C} = {\bf S}{\bf C}{\bf \epsilon} \\
h_{\alpha\beta} = \braket{\chi_\alpha | h_{KS} | \chi_\beta} \\
S_{\alpha\beta} = \braket{\chi_\alpha | \chi_\beta}
$$
* Since $V_{xc}$ functional depends on the integration variable implicitly via $\rho$, these integrals cannot be evaluated analytically. For this reason, a numerical integration is often performed thus we need to set up some grid points.

<!-- Jensen end --> 

* As stated, the density functional theory (DFT) is the one approach to solve the Schrodinger equation, which is clearly different from the wave-function-based methods discussed in the previous chapter. The fundamental object in the DFT is the electron density, and not the wave function.
* DFT originated from the work by Llewellyn H. Thomas and Enrico Fermi in 1927. In this work, using the uniform electron gas as the model, they proposed the kinetic energy as the functional of electron density.[Math Proc Camb., 1926, 21, 542–6; ZS f Phys, 1928, 48, 73–9] This is known as the Thomas-Fermi model.

* The electron density is rather expressed as a sum over single-particle states
$$
n({\bf r}) = \sum_{i=1}^N |\psi_i({\bf r})|^2
$$
* **The first Hohenberg-Kohn theorem** states that the ground-state energy from the Schrodinger equation is a unique functional of the electron density. This theorem states that there exists a one-to-one mapping between the ground-state wave function and the ground-state electron density.
* **The second Hohenberg-Kohn theorem** states the electron density that minimizes the energy of the overall functional is the true electron density corresponding to the full solution of the Schrödinger equation.
* The energy functional can be written as
$$
\begin{split}
E[\left\{ \psi_i\ \right\}] 
 =& \frac{\hbar^2}{m}\sum\int\psi_i^* \nabla_i^2 \psi_i d^3 r \\
  &+ \int V({\bf r})n({\bf r})d^3 r\\
  &+ \frac{e^2}{2}\int\int\frac{n({\bf r})n({\bf r}')}{|{\bf r}-{\bf r}'|}\\
  &+ E_{ion}\\
  &+ E_{xc}[\left\{ \psi_i \right\}]
\end{split}
$$
These terms are electron kinetic energy, e-N Coulomb, e-e Coulomb, N-N Coulomb, and exchange-correlation energy, respectively. We can write down the term in a simple analytic form, except for the $E_{xc}$.

* Now we make use of the variational principle for the energy functional and minimize $E[n]$ with respect to the single particle states under the constraint of normalization. This procedure is entirely equivalent to the derivation of the Hartree and the Hartree-Fock equations. By this, we obtain the **Kohn-Sham equation**
$$
\left\{ -\frac{\hbar^2}{2m}\nabla^2 + v_{eff}({\bf r}) \right\} \psi_i({\bf r}) = \epsilon_i \psi_i({\bf r})
$$
* The effective one-electron potential acting on the electrons is given as
$$
v_{eff}({\bf r}) = v_{ext}({\bf r}) + v_H({\bf r}) + v_{xc}({\bf r})
$$
* $v_{ext}$ is the external potential, which describes the e-N Coulomb interaction. (CHECK)
* The Hartree potential
$$
v_{H} = e^2 \int d^3r' \frac{n({\bf r'})}{|{\bf r}-{\bf r'}|}
$$
describes the Coulomb repulsion between the electrons being considered in one of the Kohn-Sham equations and the total electron density defined by all the electrons in the problem.
* $v_H$ includes a so-called self interaction, which is a problematic issue present in the DFT (but not in wave-function-based theory). The reason for this is because the electron we are describing in the Kohn-Sham equation is also part of the total electron density, thus $v_H$ involves the Coulomb interaction between the electron and itself. The self-interaction is unphysical, and the correction for it is partially done with $v_{xc}$. Unfortunately, the perfect solution to the self-interaction problem is not obtained yet.
* The exchange-correlation potential $v_{xc}({\bf r})$ is the functional derivative of the **exchange-correlation (XC) functional** $E_{xc}[n]$ with respect to the electron density $n$.
$$
v_{xc}({\bf r}) = \frac{\delta E_{xc}[n]}{\delta n}
$$
* From this, $E_{xc}[n]$ can be written as an integral form like
$$
E_{xc}[n] = \int d^3 r n({\bf r}) \epsilon_{xc}[n]({\bf r})
$$
where $\epsilon_{xc}[n]({\bf r})$ is the exchange-correlation energy per particle at the point ${\bf r}$, but depends on the whole electron density distribution $n({\bf r})$.
* Historically, a lot of XC functionals has been proposed by many researchers. These functionals are ranked by the variables used in the XC functional, and sometimes expressed as **Jacob's ladder** (Fig.X) starting from the earth ("Hartree world") to the heaven (chemical accuracy). For example, the first level is the local density approximation (LDA), which uses only the electron density $n$ as variable. The second one is generalized gradient approximation (GGA), which uses $n$ and the gradient of the electron density ($\Delta n$). The third one is called meta-GGA, where the kinetic energy density $\tau$ is added to the variables. The fourth one is hybrid GGA, which uses the Hartree-Fock exchange (or exact exchange) information. Up to this level is commonly used, while the higher ladder is possible if the state-of-art methods like double hybrid, random phase approximation (RPA) etc. is used.
* In the following, we look into the XC functionals at each level of the Jacob's ladder.

##### LDA
* There is one case where the XC functional can be derived exactly: this situation is called the uniform electron gas. In this case, the electron density is constant at all points in the space; that is, $n({\bf r}) = {\rm constant}$. This approximation uses only the local density to define the XC functional, so it is called the **local density approximation (LDA)**, which is the simplest approximation to the true XC functional.
$$
E_{xc}^{LDA}[n] = \int d^3 r n({\bf r}) \epsilon_{xc}^{LDA}\left( n({\bf r}) \right)
$$
* For some regimes, these results were determined from quantum Monte Carlo calculations, as this gives the accurate solution of the Schrodinger equation, although it is quite computationally expensive.
* The LDA is exact in the uniform density limit, which says that in the limit of uniform density the XC energy should be the exact energy of the uniform electron gas at that density.
* In a wide range of bulk and surface problems, the LDA has been surprisingly successful. However, the LDA results are not sufficiently accurate for atoms and molecules. This is related to the constancy of the electron density, because the LDA is exact in the uniform density limit. Thus the LDA is more accurate when electron density is slowly varying, which is not at all true for atoms and molecules. Practically, the LDA often shows over-binding, i.e. adsorption energies or cohesive energies being too large compared with experimental values. This leads too short lattice constant in solid case, and too short bond length in molecules. The band gap is also inaccurate, as the deviation from experimental value amount to be ~50% in cases.

##### GGA
* The next-generation XC functional to the LDA is the **generalized gradient approximation (GGA)**. The physical idea behind the GGA is simple; real electron densities are not uniform, so including information on the spatial variation in the electron density would make a functional more flexible.
* In the GGA, the XC functional is expressed using both the local electron density $n({\bf r})$ and the gradient in the electron density $\nabla n({\bf r})$;
$$
E_{xc}^{GGA}[n] = \int d^3 r n({\bf r}) \epsilon_{xc}^{GGA}\left( n({\bf r}), |\nabla n ({\bf r})| \right)
$$
* Because there are many ways to include $\nabla n({\bf r})$ in a XC functional, a lot of GGA functional have been proposed. Two of the most widely used functionals in the solid-state physics calculations are the Perdew-Burke-Ernzerhof (PBE) and the Perdew-Wang (PW91) functionals. These two functionals satisfy the uniform density limit, and also the exact properties of the XC hole near nuclei.
* DFT calculations in the GGA achieve chemical accuracy (error <= 0.1 eV) for many chemical reactions.
* The application of the semiconductor needs to go beyond the GGA, mainly because of the well-known underestimation of the band gaps of semiconductors.
* Also, the GGA has the tendency to delocalize the electronic state.[Cohen,A.J.;Mori-Sańchez,P.;Yang,W.Science, 2008, 321, 792−794.]

##### Meta-GGA
* The third rung of the Jacob's ladder is defined by the **meta-GGA** functionals, which include information from $n({\bf r})$, $\nabla n({\bf r})$, and $\nabla^2 n({\bf r})$. In practice, the kinetic energy density of the Kohn-Sham orbitals,
$$
\tau({\bf r}) = \frac{1}{2}\sum_\text{occupied states}|\nabla \varphi_i({\bf r})|^2
$$
contains the same physical information as the Laplacian of the electron density, and using these quantities has a number of advantages; so $\tau({\bf r})$ may be used in the meta-GGA functionals instead of $\nabla^2 n({\bf r})$.
* Popular functionals in this level are SCAN, TPSS, MS etc., and their variants.

##### Hybrid functional
* The fourth rung of the ladder is important because it is the most common functional in quantum chemistry calculations with localized basis sets.
* The exact exchange energy can be derived from the exchange energy density, which can be written in terms of the Kohn-Sham orbitals as
$$
E^{ex}({\bf r})
 = -\frac{1}{2n({\bf r})}\int d^3 r' \frac{|\sum_{occ}\varphi_i^{*}({\bf r}')\varphi_i({\bf r})|^2}{|r-r'|}
$$
* A critical feature of this quantity is that it is nonlocal, that is, a functional based on this quantity cannot be evaluated at one particular spatial location unless the electron density is known for all spatial locations.
* Functionals that include contributions from the exact exchange energy within a GGA functional are called the **hybrid functional**.
* Hybrid functionals describe exchange using a mixture of the exact exchange and a GGA exchange functional. By far the most widely used of these functionals in for solids systems are PBE0 and HSE, and for molecular systems the B3LYP functional is most widely used. For example, the B3LYP functional has the following form;
$$
\begin{align*}
\epsilon_{xc}^{B3LYP} = \epsilon_{XC}^{LDA} &+ \alpha_1 (E^{ex} - \epsilon_{X}^{LDA}) \\
    &+ \alpha_2 (\epsilon_{X}^{GGA} - \epsilon_{X}^{LDA}) \\
    &+ \alpha_3 (\epsilon_{C}^{GGA} - \epsilon_{C}^{LDA})
\end{align*}
$$
* Here, $\epsilon_{X}^{GGA}$ is the Becke 88 exchange functional, $\epsilon_{C}^{GGA}$ is the Lee-Yang-Parr correlation functional, and $\alpha_i (i = 1, 2, 3)$ are three numerical parameters. The three parameters were chosen empirically to optimize the performance of the functional for a sizable test set of molecular properties.
* Because of the numerical details associated with solving the Kohn-Sham equations in a plane-wave basis set, introducing the non-locality of exact exchange greatly increases the numerical burden of solving these equations. This difficulty is not so severe when localized basis sets are used, thus hybrid functional is relatively common in quantum chemistry community than solid-state physics community.
* This numerical difficulty is partly resolved by using the screened hybrid functionals. In this approach, the exchange interaction is split into two components, a long-range and a short-range one. A degree of the exact exchange is applied only to the short-range portion. The HSE functional is based on this approach, and is also widely used in plane-wave basis set.

#### van der Waals interaction
* Dispersion interaction (also known as van der Waals interaction or London force) is a direct result of long-range electron correlation. This interaction plays an incredibly important role in our everyday lives. Consider the gasoline or diesel that all of us rely on as transportation fuels. The huge global infrastructure that exist to deliver and use these fuels heavily relies on the fact that they are liquids. The fact that nonpolar molecules such as hexane readily form liquids is a signature of the net attractive interactions that exist between hexane molecules. These attractive interaction arise directly as a result of dispersion interaction.
* The relationship between electron correlation and long-range forces between atoms was initially examined in the 1930s by London. He realized that although the time-averaged electron density around an atom or nonpolar molecules has no dipole moment, electron oscillations lead to deformations of the density resulting in a transient dipole moment.
* This instantaneous dipole moment can induce a temporary dipole moment (called induced dipole moment) on other atoms or molecules by distorting their electron density. London showed that the general form of the interaction between two spherically symmetric atoms at large distance was
$$
V^{\rm dispersion} = -\frac{C}{r^6}
$$
where $r$ is the distance between the atoms and $C$ is a collection of physical constants.
* Several studies illustrate the limitations of DFT calculations with respect to the dispersion interactions.
* One conceptually simple remedy for the shortcomings of DFT regarding dispersion forces is to simply add a dispersion-like contribution to the total energy between each pair of atoms in a material.
* This method is called the **DFT-D** method.
* In the DFT-D calculations, the total energy of a collection of atoms as calculated with DFT, $E_{\rm DFT}$, is augmented as follows:
$$
E_{\rm DFT-D} = E_{\rm DFT} - S\sum_{i \ne j}\frac{C_{ij}}{r_{ij}^6}f_{\rm damp}(r_{ij})
$$
* Here $r_{ij}$ is the distance between atoms i and j, $C_{ij}$ is a dispersion coefficient for atoms i and j, which can be calculated directly from tabulated properties of individual atoms, and $f_{\rm damp}(r_{ij})$ is a damping function to avoid unphysical behavior of the dispersion term for small distances. $S$ is the empirical scaling factor applied uniformly to all pairs of atoms.

#### DFT + U
* Other than the XC functionals introduced above, another popular DFT method is the DFT + U method. This should be considered as a variant (or branch) from the Jacobs' ladder, because the empirical parameter is used.
* DFT + U method is often used in the materials with the strong correlation, because the GGA often fails to properly describe the electronic property of these systems such as band-gap.
* Many semi-conductors are grouped in this category.
* The DFT + U method is usually used with the GGA functional, and the deficiency of the GGA is corrected via the $U$ parameter.

* The $U$ parameter is called on-site Coulomb parameter. The $U$ in GGA + U refers to the strong correlation repulsion energy of electrons with opposite spins in the Hubberd model.

#### Random phase approximation
* 

### 1-3. Other topics
#### Basis set
##### Choice of the basis function form
* In solving the Schrodinger equation, we need to express the wave function by some known mathematical functions. In other words, we need "expand" the wave function by some analytic functions. These analytic functions are mathematically called **basis set** or **basis functions**.
* In principle, any type of basis functions may be used; exponential, Gaussian, polynomial, wavelets, plane waves, etc.
* There are, however, two guidelines for choosing the basis functions
1. they should have a behavior that agrees with the physics of the problem. For this reason, the basis function is often the exact solution for the simple model systems; for the periodic systems, the plane waves ($\sin$ or $\cos$ functions) are the exact solution and these functions satisfies a necessary condition for the system i.e. periodicity of the bulk system. For atomic systems, the exponential functions is the exact solution of the hydrogen atom, and it satisfies a necessary condition for the electronic wave function of atoms; it should decay toward zero when the distance between the nucleus and the electron becomes large. Thus, these functions or similar ones are often used as the basis set, and we need to select an appropriate basis set for the problem in hand.
2. the second condition is a practical one: the chosen functions should make it easy to calculate all the required integrals. For this reason, the Gaussian type orbitals are preferred over the Slater type orbitals, as will be discussed later.

<!-- Tsuneda start -->
* For contracted Gaussian-type basis functions, various types of functions have been suggested. The following are several major Gaussian-type basis functions.

* *Minimal basis functions* (e.g. STO-$L$G) contain only minimal required contracted Gaussians for each atom. For example, since electrons occupy atomic orbitals up to the 2p orbital in the carbon atom, five contracted Guassian-type functions corresponding to 1$s$, 2$s$, 2$p_x$, 2$p_y$, and 2$p_z$ orbitals are necessary at the minimum. The minimal basis functions approximating Slater-type oribtals (STO) corresponding to atomic orbitals with $L$ primitive functions are called *STO-LG* basis functions.
* *Split valence basis functions* (e.g. 6-31G, 6-311G, DZ, TZ, and DVZ) use one type of contracted Gaussians for core orbitals, while multiple contracted Gaussians are used for valence orbitals. In most molecules, valence orbitals mainly contribute to chemical bonds, but core orbitals hardly participate in the bonds. It is threfore reasonable to use many basis functions for valence orbitals to calculate electronic states accurately while keeping the total number of basis functions compact. In this type of basis function, *Pople-type* basis functions such as 6-31G basis are included. "6-31" indicates the extent of the contraction and split, in which "6" means the use of contracted basis functions of 6 primitive functions for core orbitals, and "31" means the use of doubly-split basis functions combining contracted basis functions of 3 primitive functions with one uncontracted basis functions for valence orbitals. "6-311G" uses triply-split basis functions for valence orbitals, where one contrarcted ("6" part) and two uncontracted basis functions with different exponent values ("11" part) are included. Obviously, triply-split basis (denoted triple zeta) functions gives more accurate results thant the doubly-split one, but the computational cost increases as the number of basis functions increases. Among the split valence basis functions, widely used Dunning-Huzinaga type and Ahlrichs-type basis functions are also included: the formers are described "DZ", "TZ", and "QZ" for double-zeta, triple-zeta, and quadruple-zeta functions for valence orbitals, while the latter types are similarly written as "VDZ" and "VTZ" (or "TVZ").
* *Polarization-function-supplemented basis function* (e.g. 6-31G(d), DZP, pVDZ, and cc-pVXZ) add polarization functions to incorpolate the anisotropic nature of molecular orbitals originating from chemical bonds. Polarization functions usually have higher angular momenta than the highest angular momenta of the atomic oribtals that mainly makeup the molecular orbitals. In the Pople-type basis functions, the inclusion of polarization function is denoted as the function in the parenthesis, such as "6-31G(d)". The polarization function is also denoted by an asterisk "*", such as "6-31G*". When p orbital functions are also added to the hydrogen atoms, the basis is denoted as "6-31G(d,p)" or "6-31G**". In the case of the Dunning-Huzinaga basis set, polarization functions are represented by "P", for example "DZP" for a single polarization functions, while it is denoted like "pVDZ" for Ahlrichs basis set. Recently, highly accurate calculations often use the correlation consistent basis functions. Dunning et al. developed these basis functions so that they produce electron correlation equivalently among the added polarization functions. Tehse are represented as "cc-pVXZ (X = D, T, Q, 5, 6, ...) where X denotes the level of valence splitting.
* *Diffuse-function-augumented basis functions* (e.g. 6-31+G and aug-cc-pVXZ) employ diffuse functions to take weakly bound electrons into consideration. These diffuse functions are often required in the calculations of anions or excited states of small molecules as in the case where the Rydberg states are included. Adding diffuse functions is shown by a plus "+" ni Pople-type basis functions (e.g. "6-31+G"). Double-plus "++" means diffuse functions are also added to the hydrogen atoms. For correlatino consistent basis set, augumenting diffuse function is represented by "aug" at the head of the names, and one diffuse function is added to each angular momentum type of basis function. For example of the "aug-cc-pVDZ" basis for the carbon atom,the original "cc-pVDZ" contains up to d orbital functions, and one diffuse functions of s, p, and d orbital functions are added these.
* *Effective core potential (ECP) basis functions* (e.g. LanL2DZ and Stuttgart ECP) approximate core orbitals, which hardly affect chemical reactions or properties in most systems, as effective potentials are introduced to drastically reduce the number of basis functions. In particular, for the fourth-period atoms or later, ECP basis functions are used in most cases, except for the cases where core electrons explicitly play some roles. The most widely used ECP is the relativiestic ECP, incorpolating the relativistic effects of the core electrons. The ECP basis functions are the valence basis functions used together with the ECP, and they are specifically optimized to the combined ECPs. The best known ECP basis functions are LanL2DZ of Los Alamos National Labolatory, and the Stuttgart relativistic ECP basis functions. In particular, the Stuttgart basis functions include the correlation consistent ECP basis functions (cc-pVXZ-PP), and these are often used for accurate calculations including the transition metal elements.
<!-- Tsuneda end -->

##### Atomic orbitals and the Roothaan-Hall equation.
* Each MO $\phi$ is expanded by the basis functions $\chi$, conventionally called **atomic orbitals**
$$
\phi_i = \sum_{\alpha}^{N_{bas}}c_{\alpha i}\chi_{\alpha}
$$
* The Hartree-Fock equation may be written as
$$
{\bf F}_i\sum_{\alpha}c_{\alpha i}\chi_{\alpha} = \epsilon_i \sum_{\alpha}c_{\alpha i}\chi_{\alpha}
$$
* Multiplying from the left by a specific basis function and integrating yields the *Roothaan-Hall equations*. These are the Hartree-Fock equations in the atomic orbital basis, and all the $N_{bas}$ equations may be collected in a matrix notation.
$$
\begin{align*}
{\bf FC} &= {\bf SC \epsilon} \\
F_{\alpha \beta} &= \braket{\chi_{\alpha}|{\bf F}|\chi_{\beta}} \\
S_{\alpha \beta} &= \braket{\chi_{\alpha}|\chi_{\beta}}
\end{align*}
$$
* The ${\bf S}$ matrix contains the overlap elements between basis functions, and the ${\bf F}$ matrix contains the Fock matrix elements.

* There are two types of basis functions in atomic or molecular electronic structure calculations; Slater type orbitals (STO) and Gaussian type orbitals (GTO).
* STOs have the functional form as
$$
\chi_{\xi,n,l,m}(r, \theta, \varphi) = N Y_{l,m}(\theta,\varphi)r^{n-1}\exp(-\xi r)
$$
* Here $N$ is a normalization constant and $Y_{l,m}$ are spherical harmonic functions. The exponential dependence on the distance between the nucleus and electron comes from the nature of the exact wave function of the hydrogen atom.
* Although STOs have clear physical origin, the calculation of three- and four-electron two-electron integrals cannot be performed analytically.
* GTOs can be written in terms of polar or Cartesian coordinates as
$$
\begin{align*}
\chi_{\xi,n,l,m}(r,\theta,\varphi) &= N Y_{l,m}(\theta,\varphi)r^{2n-2-l}\exp(-\xi r^2) \\
\chi_{\xi,n,l,m}(x,y,z) &= N x^{l_x}y^{l_y}z^{l_z}\exp(-\xi r^2)
\end{align*}
$$
* The sum of $l_x$, $l_y$, and $l_z$ determines the type of the orbital; $l_x + l_y + l_z = 0, 1, 2, 3$ are s-, p-, d- and f-orbitals, respectively.
* Because of the $r^2$-dependence in the exponential of GTOs, they have poorer description than STOs for the "cusp" (discontinuous derivative) at $r = 0$. This is often alleviated by increasing the number of GTO functions with high exponent $\xi$ values (thus "sharp" function at small $r$).

* Having decided on the type of the functions (STO or GTO), the most important factor is the number of functions to be used. The smallest number of functions possible is a minimal basis set. Only enough functions are employed to contain all the electrons of the neutral atoms. For hydrogen and helium, this means a single s-function. For the first row in the periodic table, it means two s-functions (1s and 2s) and one set of p-functions (2px, 2py, and 2pz).
* The next improvement of the basis set is a doubling of all basis functions, producing ad Double Zeta (DZ) type basis. The chemical bonding occurs between valence orbitals. Doubling the 1s-functions in for example carbon allows for a better description of the 1s-electrons. However, the 1s-orbital is essentially independent of the chemical environment.
* A variation of the DZ type basis only doubles the number of valence orbitals, producing a split valence basis. In actual calculations, a doubling of the core orbitals would rarely be considered, and the term DZ basis is used also for split valence basis set or sometimes referred as VDZ, for valence double zeta.
* In most cases, higher angular momentum functions are also important, and these are denoted *polarization functions*.
* If methods including electron correlation are used, higher angular momentum functions are essential.
* Polarization functions can be added to the chosen basis set. Adding a single set of polarization functions (p-functions on hydrogens and d-functions on non-hydrogen atoms) to the DZ basis forms the *Double Zeta plus Polarization (DZP)* type basis set.
* There is a variation where the polarization functions are only added to non-hydrogen atoms. As hydrogen often accounts for a large number of atoms in the system, a saving of three (px, py, and pz) basis functions for each hydrogen is significant. If hydrogen plays an important role in the property of interest, it is of course not a good idea to neglect polarization functions on hydrogen.
* In some cases, extra basis functions called *diffuse functions* are added to the basis set. These are functions with small exponents, and are needed whenever loosely bound electrons are present (for example, anions or excited states).

##### Contracted basis set
* The inner-shell or core electrons are more energetically important than the valence electrons, as these are closer to the nucleus. However, chemistry is mainly dependent on the valence electrons because chemical bonds are always formed using the valence electrons.
* The fact that many basis functions focus on describing the energetically important, but chemically unimportant, core electrons is the foundation for contracted basis set.
* Consider for example a basis set consisting of ten s-functions (and some p-functions) for a carbon atom. Optimization of these ten exponents by variational calculation is presumably leads the basis set well describing the 1s-orbital.
* The important chemical regions is however the outer valence. Out of the ten functions, a few basis functions are actually used for describing the chemically interesting phenomena.
* As the core orbitals change very little depending on the chemical bond situation, the MO expansion coefficients in front of these inner basis functions also change very little. The majority of the computational effort is therefore spent in describing the chemically unimportant part of the wave function, which is a problematic.
* Considering now making the variational coefficient in front of the inner basis functions constant, i.e. they are no longer parameters to be determined by the variational principle. The 1s-orbital in this case is, thus, described by a fixed linear combination of, say, six basis functions.
* Similarly, the remaining four basis functions may be contracted into only twos functions, for example by fixing the coefficient in front of the inner three functions. In doing this, the number of basis functions to be handled by the variational principle has been reduced from ten to three.
* Combining the full set of basis functions, known as the *primitive GTOs*, into a small set of functions by forming *fixed* linear combinations is known as *basis set contraction*, and resulting functions are called the *contracted GTOs*.
$$
\chi({\rm contracted}) = \sum_i a_i \chi_i({\rm primitive})
$$

##### Plane wave basis functions
* The outer valence electrons in metals behave almost like free electrons, which leads to the idea of using solutions for the free electrons as basis functions. The solutions to the Schrodinger equation for a free electron in one dimension can be written either in terms of complex exponential or cosine and sine functions.
$$
\phi(x) = A\exp(ikx) + B\exp(-ikx) \\
\phi(x) = A\cos(kx) + B\sin(kx) \\
E = \frac{1}{2}k^2
$$
* For infinite systems, the molecular orbitals coalesce into bands, since the energy spacing between distinct level vanishes.
* The electrons in a band can be described by orbitals expanded in a basis set of plane waves, which in three dimensions can be written as a complex function
$$
\chi_k({\bf r}) = \exp(i{\bf k}\cdot{\bf r})
$$
* Here, $k = |{\bf k}|$. The wave vector ${\bf k}$ plays the same role as the exponent $\zeta$ in GTO or STO, and is related to the energy by means of Eq.X. ${\bf k}$ can be also thought of as a frequency factor, with high ${\bf k}$ values indicating a rapid oscillation.
* Plane wave basis functions are ideal for describing delocalized slowly varying electron densities, such as valence bands in a metal.
* Describing the core region adequately with plane waves requires a large number of rapidly oscillating functions, i.e. basis with very high ${\bf k}$ values.

* Quantum chemists have developed a particular nomenclature to describe the quality of a basis set.
* The simplest choice of just one atomic orbital per valence state is called "minimal basis set" or "single zeta". It is also possible to use multiple functions to the valence state; when two functions are used, it is called the "double zeta (DZ)" basis set. Three functions are used it becomes "triple zeta (TZ)". There also exists quadruple zeta (QZ) or quintuple zeta (5Z), but these are usually used for very accurate calculations.
* The polarization functions describe small displacements of the atomic orbitals from the nuclear centers. Then a "P" is added to the acronym of the basis set resulting in e.g., DZ2P.
* For rather delocalized states such as anionic or Rydberg excited states, diffuse functions are added.
* The plane wave basis set is used in the solid-physics type calculations.

#### Pseudopotentials 
* Systems involving elements from the lower part of the periodic table have a large number of core electrons. These core electrons are unimportant in a chemical sense, but they should be included to properly describe the valence orbitals.
* Relativistic effects further complicate matters in the lower part of the periodic table.
* These two problems may be "detoured" simultaneously by modelling the core electrons by a suitable function, and treating only the valence electrons explicitly.
* The function modeling the core electrons is usually called an *effective core potential (ECP)* in the chemistry community, while the physics community uses the term *pseudopotential (PP)*.
* There are four major steps in designing a pseudopotential;
  1. generate a good-quality all wave function for the atom. Typically, numerical Hartree-Fock, Dirac-Hartree-Fock, or DFT is used.
  2. replace the valence orbitals by a set of nodeless pseudo-orbitals. The pseudo-orbitals are designed such that they behave correctly in the outer part, but without the nodal structure in the core region.
  3. replace the core electrons by a potential, and expand them by a suitable set of analytical functions of the nuclear-electron distance. For example, a set of spherical Bessel or Gaussian function are used. The potential may be different for each angular momentum.
  4. fit the parameters of the above potential such that the pseudo-orbitals matches the all-electron valence orbitals.
* The pseudopotentials are typically characterized by a "core radius" $r_c$ (which may depend on the angular momentum of the valence orbitals).
* The potential for distances smaller than $r_c$ is described by the pseudo-wave function (some analytical function), and its first and second derivatives are required to match those of the reference wave function at $r_c$.
* It is clear that a "hard" (small $r_c$) pseudopotential require more plane wave basis functions for describing the region beyond $r_c$ than a "soft" (large $r_c$) pseudopotential, but a too large $r_c$ will deteriorate the quality of the calculated results and also make pseudopotential less transferrable.
* The *norm-conserving* pseudopotential by Hamman, Schluter, and Chaing require in addition to the above matching conditions at $r_c$ that the integral of the square of the reference and pseudo-wave from 0 to $r_c$ to be same, i.e. conservation of the "wave function norm".
* These pseudopotentials are rather "hard" and therefore require a relatively large energy cutoff for the plane waves.
* Vanderbilt proposed to relax the norm-conserving requirement to give the so-called *ultrasoft* pseudopotentials, thereby reducing the necessary number of plane waves for expanding valence orbitals by roughly a factor of two.
* The *projector augmented wave (PAW)* method is usually also considered as a pseudopotential method, although it formally retains all the core electrons.
* The PAW wave function is written as a valence term expanded in a plane wave basis plus a contribution from the region within the core radius.
* The contribution from a core region is expanded as a difference between two sets of densities, one arising from the all-electron atomic orbitals, and the other is from a set of nodeless pseudo-atomic orbitals. This terms allow the wave function within the core region to adjust for different environments. In all the applications so far, the all-electron atomic orbitals have been kept fixed from the isolated atoms.

#### Population analysis
##### Mulliken
* The electron density $\rho$ is a probability of finding an electron at a certain position ${\bf r}$, and we can calculate it for a single MO containing one electron as
$$
\rho_i({\bf r}) = \phi_i^2({\bf r})
$$
* Assuming that the MO is expanded in a set of normalized but non-orthogonal basis function $\chi$, this can be written as
$$
\phi_i   = \sum_{\alpha}^{M_{basis}}c_{\alpha i}\chi_\alpha \\
\phi_i^2 = \sum_{\alpha\beta}c_{\alpha i}c_{\beta i}\chi_\alpha \chi_\beta
$$
* Integrating and summing over all occupied MOs ($N_{occ}$) gives the total number of electrons $N_{elec}$, as
$$
\begin{align*}
\sum_i^{N_{occ}}\int\phi_i({\bf r})^2 d{\bf r} 
&= \sum_i^{N_{occ}}\sum_{\alpha\beta}^{M_{basis}}c_{\alpha i}c_{\beta i}\int\chi_\alpha\chi_\beta d{\bf r} \\
&= \sum_{\alpha\beta}^{M_{basis}}c_{\alpha i}c_{\beta i}\braket{\chi_\alpha | \chi_\beta} 
= \sum_i^{N_{occ}}\sum_{\alpha\beta}^{M_{basis}}c_{\alpha i}c_{\beta i}S_{\alpha\beta} 
= N_{elec}
\end{align*}
$$
* We may generalize this by introducing an *occupation number* $n$ for each MO. For a single-determinant wave function, $n$ is either 0, 1, or 2, and the above equation is rewritten as
$$
\sum_i^{N_{orb}}n_i\int\phi_i({\bf r})^2 d{\bf r} = \sum_{\alpha\beta}^{M_{basis}}\left( \sum_i^{N_{orb}}n_i c_{\alpha i}c_{\beta i} \right) S_{\alpha\beta} = \sum_{\alpha\beta}D_{\alpha\beta}S_{\alpha\beta} = N_{elec}
$$
* The sum of the product of MO coefficients and the occupation number is the *density matrix*, and the sum over the product of the density and overlap matrices elements is the number of electrons
* The **Mulliken population analysis** uses the ${\bf D}\cdot{\bf S}$ matrix for distributing the electrons into atomic contributions. A diagonal element $D_{\alpha\alpha}S_{\alpha\alpha}$ is the number of electrons in the $\alpha$ AO ($\chi_\alpha$), and the off-diagonal element $D_{\alpha\beta}S_{\alpha\beta}$ is (half) the number of electrons shared by AO $\alpha$ and $\beta$ (there is an equivalent $D_{\beta\alpha}S_{\beta\alpha}$ element).
* The contribution from all MOs located on a given atom A may be summed up to give the number of electrons associated with atom A.
* This requires a decision on how a contribution involving basis function on different atoms should be divided. The simplest, and the one used in the Mulliken scheme, is to partition the contribution equally between two atoms.
* The Mulliken electron population is thereby defined as
$$
\rho_A = \sum_{\alpha \in A}^{M_{basis}}\sum_\beta^{M_{basis}}D_{\alpha\beta}S_{\alpha\beta}
$$
* The net charge on atom A is then the sum of the nuclear contribution $Z_A$ and the electronic contribution $\rho_A$, noting that electrons have negative charge; $Q_A = Z_A - \rho_A$.

##### Bader (Atoms in Molecules)
* The population analysis shown in the previous section is based on the wave function. However, it is also possible to use the electron density for this purpose.
* However, the problem is the definition of an "atom" within a molecule, and we need to divide the electron density into the atomic basins, and several different schemes to this have been proposed.
* Perhaps the most rigorous way of divide a molecular volume into atomic subspaces is the **Atoms in Molecules (AIM)** method of R. Bader.
* The electron density is a function of three spatial coordinates, and it may be analyzed in terms of its topology (maxima, minima, and saddle points).
* In most cases, it is found that the only maxima in the electron density occur at the nuclei, as these are the only sources of positive charge. The nuclei tus act as *attractors* of the electron density.
* At each point in space, the gradient of the electron density points the direction of the strongest attractor nearby. This forms a way of dividing the physical space into atomic subspaces; starting from a given point in space, a series of small steps may be taken in the gradient direction until an attractor is encountered. The collection of all such points forms the atomic basis associated with the attractor (nucleus).
* Once the molecular volume has been divided up, the electron density may be integrated within each of the atomic basins to give atomic charges, and dipole moments, quadrupole moments, etc.

* For a point on a dividing surface between two atomic basins, the gradient path for such a point leads to a stationary point on the surface where the total derivative is zero.
* The basin attractor is also a stationary point on the electron density surface. The second derivative of the electron density, the Hessian, is a function of the three (Cartesian) coordinates, i.e. it is a 3x3 matrix.
* At stationary points, it may be diagonalized and the number of negative eigenvalues are determined. The basin attractor is an overall maximum, so it has three negative eigenvalues.
* Another type of stationary points i.e. saddle points are found and indeed they have chemically interesting character. Saddle points are usually found between a pair of nuclei that are "chemically bonded". Such points have a minimum in the electron density in the direction of the nuclei, while a maximum in the perpendicular directions; there is one positive and two negative eigenvalues in the Hessian. These points are known as *bond critical points*.
* A bond critical point thus correspond to a saddle point in the electron density, and thus is similar to the transition state in the potential energy surface (Section X.X).
* The second derivative of the electron density, the Laplacian $\nabla^2\rho$, provides information on where the electron density is depleted or increased. At a bond critical point, the sign of the Laplacian has been used for characterizing the nature of the bond, i.e. a negative value indicates a covalent bond, while a positive value indicates an ionic bond or a van der Waals interaction.

#### Geometry optimization
* 

##### Local optimization
* 

##### Global optimization
* 

#### Vibrational analysis
* The matrix of all the second derivatives is called **Hessian or Hesse matrix**, and it expresses a linear system of differential equations
$$
{\bf H}{\bf u}_k = \omega_k^2 {\bf M}{\bf u}_k
$$
for the vibrational eigenmode $u_k$ and their frequencies $\omega_k$ that will characterize the collective movement of the atoms.
* Using the finite difference method, the Hessian can be approximated as
$$
H_{ij} = \left. \frac{\partial E^2}{\partial u_i \partial u_j} \right|_0 = -\frac{\partial F_j}{\partial u_i}
$$
where $F_j$ is the forces.

#### Transition state search
* The transition state is connecting the two energy minima in the PES, and these PES are usually correspond to some stable molecular structure.
* A clear example of this is that the bond dissociation or formation process (Fig.X). The minimum at right corresponds to the state that has the O-H bond in this example, while the left minimum corresponds to the state with N-H bond. There is a saddle point between them, and the activation energy (Ea) can be measured by the energy difference between minima and the saddle point. This topic is also discussed in the later section.
* Since we need to evaluate the Ea for the calculation of the rate constant in the reaction rate, obtaining the transition state structure is critically important topic in the computational catalytic science and many numerical methods have been developed for this task.

##### Nudged elastic band
* The method that is most widely used for finding transition states in plane-wave DFT calculation is the nudged elastic band (NEB) method.
* This method was developed by Hannes Tonsson and co-workers as a refinement of earlier "chain-of-states" method. The aim of a chain-of-states calculation is to define the minimum energy path (MEP) between two local minima.
* If you remember that the forces on the atoms in any configuration are defined by $F = -\nabla E(r)$, where $r$ is the set of coordinates of the atom, then images can be separated into two groups; the starting and the final images are located at local minima, so $F = 0$. for all the other images, the forces on the atoms are nonzero.
* Now let us move to elastic band method. This method is based on the concept that images along the MEP should use the lowest among of energy to define a path between the two minima and that the images should be evenly spaced along the path. These two ideas can be expressed mathematically for a set of images $r_0, r_1, \cdots, r_p$ by defining the objective function
$$
M(r_1, r_2, \cdots, r_p) = \sum_i E(r_i) + \sum_i \frac{K}{2}(r_i - r_{i-1})^2
$$
* Here, $E(r_i)$ is the total energy of the i-th image, and $K$ is a constant that defines the stiffness of the harmonic springs (the "elastic band") connecting adjacent images.
* The objective function does not include the images 0 or P, because those images are held fixed at the energy minima.
* The minimization of this objective function moves all the images closer to the true MEP.
* The **nudged elastic band (NEB)** is defined to improve upon the elastic band method. From the current position of the images, we can estimate the direction of the path defined by the image.
* A useful estimate for this direction is to define th path direction for image i, $\hat{\tau_i}$ as unit vector pointing along the line defined by the two adjacent images, $r_{i+1}-r_i$.
* The images will satisfy the definition of a MEP, given above if the component of the force not pointing along the path direction is zero, that is, $F_i^{\perp} = F_i - (F_i\cdot\hat{\tau_i})\hat{\tau_i} = 0$.
* This description suggests a simple strategy for adjusting the images - move each of them "downhill" along the direction defined by $F_i^{\perp}$
* If we want to include harmonic springs between images, the nwe also need to include the spring forces in defining this downhill direction.
* A useful way to do this is to define
$$
F_{i, \rm{update}} = F_{i}^{\perp} + F_{i, \rm{spring}}
                   = F_{i}^{\perp} + K(|r_{i+1}-r_i| - |r_i - r_{i-1}|)
$$
* We want the spring forces to act only to keep the images evenly spread out along the path, and we do not want the spring forces to pull the images away from the MEP.
* To do this, we define
$$
F_{i, \rm{spring}}^{\parallel} = (F_{i, \rm{spring}}\cdot\hat{\tau_i})\hat{\tau_i}
$$
and then the update the positions using
$$
F_{i, \rm{update}} = F_i^{\perp} + F_{i, \rm{spring}}^{\parallel}
$$
* If all the images in the calculation lie on an MEP, then this update force is zero for every image and the calculation has converged.

<!--Jensen start-->

* The *Nudged Elastic Band (NEB)* method defines a target function ("elastic band") as a sum of energies of all images and adds a penalty term having the purpose of distributing the images along the path. A single spring constant $k$ will attempt to distribute the images evenly along the path, but it may also be taken to depend on the energy in order to provide a better sampling near the saddle point:
$$
T_{\rm NEB}({\bf R},{\bf x}_1,{\bf x}_2, \cdots, {\bf x}_M, {\bf P}) = \sum_{i=1}^{M}E({\bf x}_i) + \sum_{i=1}^{M-1}\frac{1}{2}k({\bf x}_{i+1}-{\bf x}_i)^2
$$
* A straightforward minimization of $T_{\rm NEB}$ gives a reaction path that has a tendency to cut courners if the spring constant $k$ is too large and a problem of image sliding donw towards the minima if the spring constant is too small.
* These problems can of course be solved by employing a large number of images, but that would render the optimization inefficient. The "corner-cutting" and "down-sliding" problemsm for a manageable number of images can be alleviated by "nudging" the elastic band, that is using only the component of the spring force parallel to the tangent of the path, and only the perpendicular component of the energy force in the optimization on $T_{\rm NEB}$.
* The parallel forces for each image are obtained by projecting the total force onto the reaction path tangent, and the perpendicular components by an orthogonal projection. Since the reaction path is presented by a discrete set of images, the tangent to the path at a give nimage must be estimated from the neighboring images.
* The magnituide of the spring constant influences the optimization efficiency; a small value causes an erratic coverage of the reaction path, while a large value forces the effort on distributing the image rather than on finding the reaction path, and consequently slows down the convergence.
* The NEB algorithm defines the tangent as the difference vector of the neighboring images.
* In the *Climbing Image NEB (CINEB)* version, one of the images is allowed to move along the elastic band, to become the exact saddle point.

<!--Jensen end-->

##### Dimer method
* 

### 1-X. Molecular dynamics
* Nuclei are heavy enough that they, to a good approximation, behaves as classical particles and their dynamics can thus be simulated by soloving the Newton's equation, ${\bf F} = m{\bf a}$, in differential form
$$
-\frac{dV}{d{\bf r}} = m\frac{d^2{\bf r}}{dt^2}
$$
* Here $V$ is the potential energy at potition ${\bf r}$. The vector ${\bf r}$ contains the coordinates for all the particles, that is in Cartesian coordinates it is a vector of length $3N_\text{atom}$. The left-hand side is the negative of the energy gradient, also called the force (${\bf F}$) on the particle(s).
* Given a set of particles with positions ${\bf r}_i$, the positions at small time step $\Delta t$ later are given by a Taylor expansion:
$$
\begin{align*}
{\bf r}_{i+1} = {\bf r}_i + \frac{\partial{\bf r}}{\partial t}\Delta t + \frac{1}{2}\frac{\partial^2{\bf r}}{\partial t^2}\Delta t^2 + \frac{1}{6}\frac{\partial^3{\bf r}}{\partial t^3}\Delta t^3 + \cdots \\
= {\bf r}_i + {\bf v}_i \Delta t + \frac{1}{2}{\bf a}_i\Delta t^2 + \frac{1}{6}{\bf b}_i\Delta t^3 + \cdots
\end{align*}
$$
* The velocities ${\bf v}_i$ are the first derivaties of the positions with respect to time at $t_i$, the accelerations ${\bf a}_i$ are the second derivatives at $t_i$, and the changes in accelerations (${\bf b}_i$) are the third derivates, etc.
* The positions at small time step $\Delta t$ earlier are derived from Eq.(15.3) by substituting $\Delta t$ with $-\Delta t$:
$$
{\bf r}_{i-1} = {\bf r}_i - {\bf v}_i\Delta t + \frac{1}{2}{\bf a}_i\Delta t^2 + \frac{1}{6}{\bf b}_i\Delta t^3 + \cdots
$$
* Addition of Eqs.(15.3) and (15.4) gives a recipe for predicting the position at time step $\Delta t$ later from the current and previons positions, and the current acceleration. The latter can be calculated from the force, or equivalently, the potneital $V$:
$$
{\bf r}_{i+1} = \left(2{\bf r}_i - {\bf r}_{i-1}\right) + {\bf a}_i\Delta t^2 + \cdots \\
{\bf a}_i = \frac{{\bf F}_i}{m_i} = -\frac{1}{m_i}\frac{dV}{d{\bf r}_i}
$$
* This is the *Verlet* algorithm for solving the Newton's equation numerically.
* Note that the term involving ${\bf b}$ disappears, that is the equation is correct to the thrid order in $\Delta t$.
* The Verlet algorithm has the numerical disadvantage that the new positions are obtained by adding a term proportional to $\Delta t^2$ to a difference in positions $(2{\bf r}_i - {\bf r}_{i-1})$. Since $\Delta t$ is a small number and $(2{\bf r}_i - {\bf r}_{i-1})$ is a difference between two large numbers, this may lead to truncation errors due to finite precision.
* The Verlet algorithm furthermore has the disadvantage that velocites do not appear explicitly, which is problem in generating ensembles with constant temperature, as discussed in next section.
* The numerical aspect and the lack of explicit velocities in the Verlet algorithm can be remediced by the *leap-frog* algorithm. Performing expansions analogous to Eqs.(15.3) and (15.4) with half a time step followed by subtraction gives
$$
{\bf r}_{i+1} = {\bf r}_i + {\bf v}_{i+1/2}\Delta t
$$
* The velocities is obtained by analogous expansions to give
$$
{\bf v}_{i+1/2} = {\bf v}_{i-1/2} + {\bf a}_i\Delta t
$$
* Eqs.(15.8) and (15.9) define the leap-frog algorithm.
* In terms of theoretical accuracy, the leap-frog is also of third-order, likewise the Verlet algorithm, but the numerical accuracy is better.
* The velocity now appears directly, which facilitates a coupling to an external heat bath.
* The disadvantage is that the positions and velocites are not known at the same time; they are always out of phase by half a time step. This problem can be removed by the *velocity Verlet* algorithm, where the atoms are propagated according to
$$
{\bf r}_{i+1} = {\bf r}_i + {\bf v}_i\Delta t + \frac{1}{2}{\bf a}_i\Delta t^2 \\
{\bf v}_{i+1} = {\bf v}_i + \frac{1}{2}\left({\bf a}_i + {\bf a}_{i+1}\right)\Delta t
$$
* The preference of Verlet and leap-frog-type algorithms over, for example, the Runge-Kutta methods for solving the differential equation in MD simulations is that they are time-reversible, which in general tend to improve the energy conservation over long simulation times.

#### Generating non-natural ensembles
* A standard MD simulation generates an NVE ensemble, that is the temperature and pressure will fractuate.
* The total energy is a sum of the kinetic and potential energies, and can be calcualted from the positions and velocities:
$$
E_{tot} = \sum_{i=1}^N \frac{1}{2}m_i{\bf v}_i^2 + V({\bf r}) \\
        = E_{kin} + E_{pot}
$$
* The temperature of the system is proportional to the average kinetic energy, as
$$
\braket{E_{kin}} = \frac{1}{2}\left(3N_\text{atom} - N_\text{constraint} \right) k_B T
$$
* The number of constraint $N_\text{constraint}$ is typically three, corresponding to conservation of linear momentum. Note that for 1 mole of particles, Eq.(15.19) reduces to the familiar expression $\braket{E_{kin}} = 3/2 RT$.
* Since the kinetic energy is the difference between the total energy (almost constant) and the potential energy (depends on the positions), the kinetic energy will vary significantly, thus the temperature will be calculated as an average value with an associated fluctuation.
* Similarly, if the volume of the system is fixed, the pressure will fluctuate.
* Although the NBE is the natural ensemble generated by an MD simulation, it is possible also to generate NVT or NPT ensemble by MD techniques by modifying the velocities or positions in each time step.
* The instant value of the temperature is given by the average of the kinetic energy, as indicated by Eq.(15.19).
* If this actual temperature $(T_\text{actual})$ is different from the desired temperature $T_\text{desired}$, all velocities may be scled by a factor of $\left(T_\text{desired}/T_\text{actual}\right)^{1/2}$ in each time step to achieve the desired temperature. Such an "instant" correction procedure actually alters the dynamics, such that the simulation no longer corresponds to a canonical (NVT) ensemble.
* Performing the scaling at larger intervals introduce some periodicity into the simulation, which is also undesirable.
* Altenatively, the system may be coupled to a "heat bath", which gradually adds or remove energy to/from the system with a suitable time constant, a procedure of ten called *thermostat*.
* The kinetic energy of the system is again modified by scaling the velocities, but the rate of heat transfer is controlled by a coupling parameter $\tau$:
$$
\frac{dT}{dt} = \frac{1}{\tau}\left(T_\text{desired} - T_\text{actural}\right) \\
\text{velocity scale factor} = \sqrt{1+\frac{\Delta T}{\tau}\left(\frac{T_\text{desired}}{T_\text{actual}}-1\right)}
$$
* Thermostat methods such as Eq.(15.20) are widely used but again do not produce a canonical ensemble; they do produce correct averages but give incorrect fluctuations of properties.
* In *Nos$\'{e}$-Hoover* methods, the heat bath is considered an integral part of the system and assigned fictive dynamic variables, which are evolved on an equal footing with other variable.s These methods are analogous to the extended Lagrange methods, and can be shown to produce true canonical ensemble.
* The pressure can similarly be held (approximately) constant by coupling to a "pressure bath". Instead of changing the velocities of the particles, the volume of the system is changed by scaling all coordinates according to
$$
\frac{dP}{dt} = \frac{1}{\tau}\left(P_\text{desired} - P_\text{actual}\right) \\
\text{coordinate scale factor} = \sqrt[3]{1+\kappa\frac{\Delta T}{\tau}\left( P_\text{actual} - P_\text{desired} \right)}
$$
* Here the constant $\kappa$ is the compressibility of the system. Such *barostat* methods are again widely used, but do not produce strictly correct emsembles. The pressure may alternatively be maintained by a Nos$\'{\rm{e}}$-Hoover type approach, in order to produce a correct ensemble.

#### Constrained and biased sampling methods
* A straightforward sampling of the reaction path is not possible since the dynamics at ordinary temperatures only very rarely visit the high-energy region near the TS(unless the activation energy is close to zero).
* In order to achieve a sampling of a specific region of the free energy surface with MD, the sampling must be giased toward a specific volume of phase space.
* A central component in biased sampling method is the selection of one or a few *collective variables (CV)* that can be used to describe the reaction path.
* A CV can be any function of atomic coordinates, but it must be selected and defined by the user based on the given application.
* A CV can be in the simplest case be the position of a point, a distance between two points, an angle between three points or torsional angle between four points, where the points can either be single atoms or center of mass for a large group of atoms.
* Once a set of CVs have been selected, the biasing can be done by two different methods, the penalty approach or the Lagrange type approach, analogously to the optimization of function with constraints.
* A number of closely related methods have been proposed for performing simulations with bias potentials, typically with the aim of calculating free energy profiles along a suitable reaction coordinate, and the following describe some of these approaches.
* The penalty approach corresponds to augumenting the energy surface with a biasing potential $\mathcal{U}$, for example a harmonic function centered at position $\xi_0$ of the CV with a suitable width $k_{\mathcal{U}}$:
$$
\mathcal{V}_{\text{umbrella}}(\xi) = \mathcal{V}(\xi) + \mathcal{U}(\xi) \\
\mathcal{U}(\xi) = k_{\mathcal{U}}\left(\xi - \xi_0\right)^2
$$
* By making the potential $\mathcal{U}$ sufficiently steep (large $k_{\mathcal{U}}$), the energy of the augmented energy surface far from $\xi_0$ will be come so high in energy that only the region near $\xi_0$ will be sample at ambient temperatures. This technique is called *umbrella sampling*.
* The ensemble calculated with the augmented potential $\mathcal{V}_\text{umbrella}$ will of course be non-Boltzmann, but this can be deconvoluted as shown for a property $P$.
$$
\braket{P} = \frac{\braket{P(\xi)\exp(\mathcal{U}/k_B T)}_{\mathcal{V}_\text{umbrella}}}{\braket{\exp(\mathcal{U}(\xi)/k_B T)}_{\mathcal{V}_\text{umbrella}}}
$$
* Here $\braket{\ }_{\mathcal{V}_\text{umbrella}}$ indicates an average over the ensemble generated by the augmented potential. By performing a series of simulations with biasing potential located at different positions along the reaction path, the free energy along the reaction path often called the *potential of mean force (PMF)*, can be calculated.
* The Lagrange approach constrains the sampling to the $N-1$ dimensional subspace corresponding to a specific value of the CV, where the constraint is fulfilled by means of an additional term in the Hamiltonian onvolving a Lagrange multiplier. This is related to the extended Lagrange techniques, and is usually referred to as *Blue moon sampling* in the literature (the term "Blue moon" denotes the relatively rare phenomenon).

#### Metadynamics
* A common feature of umbiased MD is the tendency of the simulation to revisit the same region of phase space by random thermal motion due to the presence of energy barriers substantially larger than the available thermal energy.
* **Metadynamics** attempts to prevent this situation by gradually adding repulsive Gaussian-shaped potential based on the sampling history and can be thought of as gradually filling up the energy minimum currently being sampled until the thermal energy is sufficient for escaping to another minimum.
* The bias potential is defiend in terms of CV and is for a single CV given by
$$
\mathcal{U}_\text{metadyn}(\xi) = \sum_{t_0=0, \Delta, 2\Delta, \cdots}^{t-\Delta}W\exp\left[-\frac{\left\{\xi(t)-\xi(t_0)\right\}^2}{2\sigma^2}\right]
$$
* Here $W$ and $\sigma$ are user-defined parameters controlling the height and width of the Gaussian, and $\Delta$ is the time interval between points on the trajectory where the bias potential are placed. Equation above is readily extended to include more than one C by simply including an additional term in the exponential function for each additional CV.

### 1-X. Force field
* In force field method (also referred as molecular mechanics method), calculation of the electronic energy for a given nuclear coordinate is bypassed by writing the electronic energy as a parametric function of the nuclear coordinates.
* Parameters included in this function is fitted to experimental data or higher level computational data, such as ab initio method.
* The fundamental objects in the force field method are atoms, so electrons are not considered as individual particles. This means that bonding information must be provided explicitly rather than being the result of solving the electronic Schrodinger equation.
* This also means that the quantum aspects of the nuclear motion are neglected, ans the dynamics of the atoms is treated by classical mechanics, i.e. the Newton's equation of motion.
* Molecules, bulk materials, and surfaces are all described by a "ball and spring" model in the force field method, with atom having different sizes and "softness" and bonds having different length and "stiffness".
* The idea of molecules being composed of atoms, which are structurally similar in different molecules, is implemented in force field as atom types. The atom type depends on the atomic number and the type of chemical bonding it is involved in. For example, sp2-hybridized carbon and sp3-hybridized carbon atoms are in different atom types.
* The force field energy is written as a sum of terms, each describing the energy required to distorting a molecule in a specific manner.
* $E_{str}$ is the energy function for stretching a bond between two atoms, $E_{bend}$ represents the energy required for bending an angle, $E_{tor}$ is the torsional energy for rotation around a bond. $E_{vdw}$ and $E_{el}$ describe the van der Waals and electrostatic atom-atom interactions, and finally $E_{cross}$ describe the coupling between the first three terms.

#### Advantage and limitation
* The main advantage of force field methods is the speed with which calculations can be performed, enabling large systems to be treated. Even with a laptop computers, molecules with several thousand atoms can be optimized. This puts the applications in the regions of modeling biomolecules such as proteins and DNA. Molecular modeling is now used by many pharmaceutical companies.
* For systems where good parameters are available, it is possible to make very good predictions of geometries and relative energies. One of ths main problems is of course the lack of parameters; if the molecule is slightly out of the common dataset, it is very like that only poor quality parameters exist, or not at all. This often happens for systems containing the transition metal elements.

#### Ensemble
* 

### 1-X. Kinetic Monte Carlo
* Nowadays kinetic Monte Carlo (KMC) method is a popular tool to describe a variety of phenomena related to e.g. transport (diffusion), structures and properties of materials (e.g., crystal growth) or equilibrium and non-equilibrium chemistry (catalysis).

* Many elementary processes involved at surfaces of solids exhibit high activation barriers. These barriers are usually much larger than $k_B T$ and the corresponding processes are thus classified as rare events, if only thermal energy is there to drive them.
* While the motion of of atoms occurs on picosecond time scales, the time between consecutive high-barrier events can therefore be many orders of magnitude longer.
* The "life" of our system in the long time span between these rare events is filled with vibrational motion around a single minimum on the PES.
* The relevant transitions to other (meta)stable states (or PES basins) occur only occasionally.
* On a mesoscopic time scale, the time evolution of our system therefore manifests itself as a series of consecutive jumps from state to state.
* The longer the time the system spends in one basin, the more it "forgets" how it actually got there. In other words, the state-to-state jumps of the system constitute a so-called *Markov chain*.

* The change of the probability $P_i(t)$ of the system to actually be in state $i$ at time $t$ depends only on the probabilities of hopping out of the current state $i$ into any other state $j$, which is denoted as $k_{ij}$.
* In the present context of chemical kinetics, these hopping probabilities are expressed as rate constants of the elementary processes with units ${\rm time}^{-1}$. 
* The overall change in $P_i(t)$ is thus governed by a simple balancing equation, called a **master equation**, that only contains these rate constants:
$$
$$

* As discussed above, the real trick of KMC is the KMC algorithm that generates stochastic trajectories in such a way that their appropriate averaging yields the time evolution of the probability $P_i(t)$ in the master equation.
* One of the most commonly used such KMC algorithm, initially developed by Bortz et al. in 1975 for Ising spin systems, is known as the *BKL algorithm*. This is also called as the *variable step size method* or *direct method*.
* Here for simplicity, let us consider a system in which only two states $A$ and $B$ connected by a barrier with associated rate constants $k_{AB}$ and $k_{BA}$ for the forward and backward transitions, respectively.
* In this system, only two elementary processes are possible; the system is sitting in state $A$, it can hop to $B$. As the rate constant for this is $k_{AB}$, one may naively think that the average time that will passed until such an event occurs is $\Delta t_{AB} = k_{AB}^{-1}$, as the unit of $k_{AB}$ is ${\rm time}^{-1}$. Obviously, for the hop back from state $B$ to $A$, the average time would be $\Delta t_{BA} = k_{BA}^{-1}$.
* Consequently, a KMC algorithm would generate a trajectory where after each hop the time is incremented by $\Delta t_{AB}$ or $\Delta t_{BA}$.
* Mathematically, this naive thinking is not entirely correct. In reality, while being in state $A$ for each short increment of time, the system will have the same probability of finding the escape path. This generates an exponentially decaying survival statistics, whose derivative represents the probability distribution $p_{AB}$ for the true time of first escape as
$$
$$
* The average escape time thus has to be appropriately weighted by this Poisson distribution. It can be shown that this is achieved by advancing the system clock by
$$
$$
where $\rho_2 \in [0,1]$ is a randomly drawn number.

#### Application of the KMC
* In this model, the RuO2(110) surface is considered to contain two times of active sites; bride (br) and coordinately unsaturated (cus) sites. Each site can be either empty or occupied by O or CO. A total of 26 processes are possible in this system, covering non-dissociative CO adsorption/desorption, dissociative O2 adsorption/desorption, diffusion of O and CO, and CO2 formation. The formed CO2 is assumed to desorb instantaneously and irreversibly due to its weak binding to the surface.
* The Figure X shows the temporal evolution of the system at 1 bar O2 and CO2 pressure and 450 K. The evolution starts from an empty RuO2 surface.
* In the very first KMC steps, the coverage builds up quickly; the O coverage builds up roughly double as fast as the CO coverage, as every O2 dissociative adsorption yields two O atoms.

* Once the steady-state solution has been reached, one can make use of the ergodicity of the KMC simulation to calculate the desired quantities as time averages instead of ensemble averages. Such quantities are surface composition, occurrence of various elementary steps, TOFs etc.
* The average reaction rate is, for example, can be calculated as follows;
$$
$$
* Here $t_{\rm KMC}$ is the total KMC simulation time, the first sum runs over all KMC steps $n$ (up to a total $N_{\rm KMC}$ steps), the second sum runs over all states $j$ that are accessible from the current state $i$, $k_{ij}^{\beta}$ is the rate constant for a process involving the production of the molecule $\beta$, and $\Delta t_{\rm escape, n}$ is the escape time for KMC step $n$.
* The total simulation time should be chosen long enough to reduce the statistical error on the sample quantities to a desired value.

### 1-X. Machine-Learning
* The effective use of ML not only facilitates the discovery of materials, but also helps to establish a deeper understanding of the relationships between the properties of materials and their functionalities.
* In this section, we overview the application of the ML technique to the catalytic chemistry.
* The basic tools for the ML is first presented, and the application to the catalytic chemistry is discussed later.

#### regression
* From the regression coefficient, one can identify the important descriptor for the target variables.

* Fig.X shows the heatmap of the regression coefficient on the adsorption energy(?) calculated by the DFT method.
* One can see that XXX is the most important descriptor.

#### Machine-learning force field
* To perform MD simulations, typically the Newtonian equations of motions are integrated numerically. This requires knowledge of the forces acting on individual atoms at each time step of the simulation.
* In principle, the most accurate way to obtain these forces is by solving the Schrodinger equation at each time step. However, this is prohibitively demanding in terms of computational cost and time.
* Instead, simple empirical functions are commonly used to model the relevant interactions. These are called force fields (FFs), and atomic forces can be readily derived analytically.
(Details of FF --> other part)
* Machine learning (ML) methods could help in closing the gap between the accuracy of ab initio methods and the efficiency of classical FFs.
* ML methods aim to learn the functional relationship between inputs (chemical descriptors) and outputs (properties) from patterns of structure in the data.
* Practically, ML models can take a shortcut by not having to solve any equations that follow from the physical laws governing the structure-property relation.
* For constructing ML-FFs, suitable reference data to learn the relevant structure-property relation include energy, forces, or a combination of both, obtained from ab initio calculations.

* By introducing a parametric dependency between energy and nuclei, the Born-Oppenheimer approximation implies the existence of a functional relation $f: \left\{ Z, r \right\}_{i=1}^N \rightarrow E$, which maps the nuclear charges $Z_i$ and positions ${\bf r}_i$ of $N$ atoms directly to their potential energy $E$.
* This function, called the potential energy surface (PES), governs the dynamics of a chemical system.
* At each time step of a dynamics simulation, the forces ${\bf F}_i$ acting on each atom $i$ must be known so that the equation of motion can be integrated numerically. They can be derived from the PES by using the relation ${\bf F}_i = -\nabla_{{\bf r}_i}E$, that is, the forces are the negative gradient of the potential energy $E$ with respect to the atomic positions ${\bf r}_i$.

##### NN
* In the simplest case, the fundamental building blocks of neural networks (NNs) are dense (or "fully-connected") layers; linear transformations from input vectors ${\bf x}\in \mathbb{R}^{n_{in}}$ to output vectors ${\bf y} \in \mathbb{R}^{n_{out}}$ according to
$$
{\bf y} = {\bf W x} + {\bf b}
$$
where both weights $W$ and biases $b$ are parameters, and $n_{in}$ and $n_{out}$ denote the number of dimensions of ${\bf x}$ and ${\bf y}$, respectively.
* Evidently, a single dense layer can only express linear functions.
* Nonlinear relations between inputs and outputs can only be modeled when at least two dense layers are stacked and combined with a nonlinear activation function $\sigma$:
$$
$$
* Provided that the number of dimensions of the "hidden layer" $h$ is large enough, this arrangement can approximate any mapping between inputs $x$ and outputs $y$ to arbitrary precision, i.e. it is a generation function approximator.
* In theory, shallow NNs as shown above are sufficient to approximate any functional relationship. However, deep NNs with multiple hidden layers are often superior and were shown to be more parameter-efficient.[153-156].
* To construct a deep NN, L hidden layers are combined sequentially
$$
$$
mapping the input ${\bf x}$ to several intermediate feature representations $h_l$, until the output $y$ is obtained by a linear regression on the features $h_L$ in the final layer.
* For PES construction, typically, the NN maps a representation of chemical structure ${\bf x}$ to a one-dimensional output representing the energy.
* The parameters ${\bf W}$, ${\bf b}$ of and NN cannot be fitted in closed form. Instead, they are initialized randomly and optimized (usually using a variant of stochastic gradient descent) to minimize a loss function that measures the discrepancy between the output of the NN and the reference data, such as the means squared error (MSE).

##### Neutral network potentials
* The first neural network potentials (NNPs) used a set of internal coordinates, for example, distances and angles, as structural representation to model the PES.[183-187]
* Behller and Parrinello were the first to propose so-called high-dimensional neutral network potentials (HDNNPs) where the total energy of a chemical system is expressed as a sum over atomic contributions.[114]
* The total energy $E$ of the system is obtained as the sum over all $N_{\rm atoms}$ atoms in the system, as
$$
E = \sum_{\mu=1}^{N_{\rm atoms}}E_{\mu}
$$
* For a multi-component system containing $N_{\rm elem}$ elements, this equation becomes a sum over all the atomic energy contributions of the elements.
$$
E = \sum_{\nu=1}^{N_{\rm elem}}\sum_{\mu=1}^{N_{\rm atoms}}E_{\mu}^{\nu}
$$
* For a given element, the atomic NNs are constrained to have the same architecture -- specifying the number of hidden layers and neurons -- and the same weight parameters.
* The input of each atomic NN is a vector of atom-centered symmetry functions, which describes the local chemical environments of the atoms. These are defined by a cutoff radius $R_c$.

* Several types of symmetry functions are available. These are all many-body functions that depend simultaneously on the positions of all the atoms inside $R_c$. Their numerical values are invariant with respect to rotation and translation as well as the order of the neighboring atoms, and, thus, possess all the required invariances.
* At the cutoff radius, they decay to zero in both the value and slope, according to the cutoff function
$$
f_c(R_{ij}) =
\left\{
\begin{align*}
&0.5\left[\cos\left(\frac{\pi R_{ij}}{R_c}\right)+1\right] & &(R_{ij} \leq R_c) \\
&0 & &(R_{ij} \gt R_c)
\end{align*}
\right.
$$
where $R_{ij}$ is the distance between the central atom $i$ and its neighbor $j$.
* The radial symmetry function is
$$
G_i^2 = \sum_j \exp\left[-\eta(R_{ij}-R_S)^2\right]\cdot f_c(R_{ij})
$$
* This is a sum of products of a Gaussian function of the interatomic distance and the cutoff function $f$. Radial functions with different Gaussian exponents $\eta$ provides a radial fingerprint of the neighboring atoms. The parameter $R_S$ can be used to shift he centers of the Gaussians to specific interatomic distances.
* Since the radial functions alone are unable to distinguish different angular arrangements of neighbors, a set of "angular functions" should be used, which is
$$
G_i^4 = 2^{1-\zeta}\sum_j\sum_{k\neq j}(1+\lambda\cos\theta_{ijk})^{\zeta}\cdot f_c(R_{ij})f_c(R_{ik})f_c(R_{jk})
$$
* This depends on the angles $\theta_{ijk}$ centered at atom $i$ and formed with neighbors $j$ and $k$, which both needed to be within $R_c$. The use of a set of functions with different exponents $\zeta$ allows a fingerprint of the angular distribution to be obtained, while $\lambda = \pm 1$ can be used to adjust the positions of the maxima and minima of these functions.

* The underlying assumption is that the energetic contribution $E_i$ of each atom depends mainly on its local chemical environment.
* The introduction of HDNNPs inspired many NN architectures that can be broadly categorized into two types; descriptor-based NNPs [116, 191-193] or the end-to-end NNPs.[160, 196-198].

###### steps
* The construction of any MLP consists of several steps.
* First, an electronic structure method is chosen as a reference, which is able to capture the correct physical behavior of the system. This is very important, because the goal of the subsequent construction of the MLP is to reacy a close numerical agreement of the energies (and forces) with this reference method.
* The construction of this reference set is often the computational bottleneck in the construction of MLPs.
* The structures must be carefully chosen to ensure that the important features of the PES are present in the data set.
* In a second step, the data are prepared for training the ML method. This is done by transforming the atomic positions into a special set of coordinates that has to meet several requirements. It is followed by the fitting process, in which the parameters of the ML method are adjusted to reproduce the reference data as accurately as possible. This is an essential step, because the functional forms of ML methods have no physical meaning, and the correct shape of the PES has to be learned from the available electronic structure data.
* Once the MLP is carefully validated, it is ready for applications.
* The training of the NN parameters can be carried out using total energies and, if desired, using forces, which provide valuable local information about the shape of the PES.
* As a consequence of the locality of the atomic interactions, it is possible to train HDNNPs by using small systems containing only a few hundred atoms, and once constructed, they can be applied to much larger systems.

##### Descriptor-based NNPs
* The first descriptor-based NNP introduced by Behller and Parrinello uses atom-centered symmetry functions (ACSFs) consisting of two-body terms $G_i^2$ and three-body terms $G_i^3$.
* When sufficiently many $G_i^2$ and $G_i^3$ with different parameters are combined and stored in a vector $x$, they form a "fingerprint" of the local environment of atom $i$. This environment descriptor is then used as input for a NN for predicting the energy contributions $E_i$ of atoms $i$ and the total energy $E$ which is obtained by summation.
* Since the ACSFs only use geometric information, they work best for systems containing only atoms of one element, for example crystalline silicon.[114] To describe multi-component systems, the symmetry functions are duplicated for each combination of elements and separate NNs are used to predict the energy contributions for atoms of the same type.[205]
* Most descriptor-based NNPs, such as ANI [194] and TensorMol,[195] use variations of $G_i^2$ and $G_i^3$ or Behller-Parrinello to construct the environment descriptors $x_i$.
* The common feature for all variations of the descriptor-based NNPs is that the functional form of the environment descriptor is predetermined and manually designed.

##### End-to-End NNs
* A potential drawback of the previously introduced ACSFs is that they must be chosen by an expert before training the NNs. If the choice of symmetry functions is poor, for example when the resulting descriptor is (nearly) identical for two very different structures, the expressive power of the NNs and the achievable accuracy are limited.
* Additionally, a glowing number of input dimensions can quickly become computationally expensive, both for calculating the descriptors and for evaluating NN. This is especially the case when modeling multi-component systems.
* In contrast to this, end-to-end NNPs directly take atomic types and positions as inputs to learn suitable representations from the reference data.
* Likewise the descriptor-based NNPs, many end-to-end NNPs obtain the total energy $E$ as a sum of atomic contributions $E_i$. However, those are predicted from learned features $x_i$ encoding information about the local chemical environment of each atom $i$. This allows them to adapt the features based on the size and distribution of the training set as well as the chemical property of interest during the training process.
* Within the deep tensor neural network framework,[160], this is achieved by iteratively refining the atomic features $x_i$ based on neighboring atoms.
* The atomic features $x_i$ can be written in general as
$$
$$
* Here the summation runs over all atoms within a distance $r_{\rm cut}$ and a cutoff function $f_{\rm cut}$ ensures smooth behavior when atoms cross the cutoff. The "atom-wise" function $G^t$ is used to refing the atomic features after they have been updated with information from neighboring atoms through the interaction function $F^t$.
* Both $F^t$ and $G^t$ functions are NNs, and varing between specific end-to-end NNPs.
* As only pairwise distances are used and the order of atoms is irrelevant due to the commutative property of summation, the features $x_i$ by Eq.X are automatically rototranslationally and permutationally invariant.
* Glimer et al. have cast graph networks of this structure as message-passing neutral networks and proposed a variant that uses a set2set decoder instead of a sum over energy contributions.[198, 208]
* SchNet takes an alternative view of the problems and models interactions between atoms with convolutions.[109]
* PhysNet modified the energy function to include explicit terms for electrostatic and dispersion interaction.[108]

###### Data collection
* A good starting point to assemble the reference data set is by sampling the PES using ab initio molecular dynamics (AIMD).
* Here, the temperature of the simulation determines which regions of the PES and what energy ranges are explored. Sampling at higher temperatures ensures that the model does not enter the "extrapolation regime", in which the NNPs are trained sufficiently. However, pure AIMD sampling is only advisable when the intended application of the final ML model involves the MD simulation for equilibrium or close to equilibrium properties, where rate events do not play a major role.

### 1-X. Electronic Structure Analysis and Catalysis
* Computational chemistry provides many types of information: the electronic state, energetics, structures, optical properties, magnetic properties, and so on.
* Finding the relationship among these properties is helpful for understanding and predicting the functions of materials.
* For example, the relationship between electronic structure and the energetics is useful for identifying the adsorption/desorption event or chemical reactivity.
* One of the most commonly used used way for this the d-band center theory.

#### Sabatier principle
* Optimal catalysts should bind the reactants strongly to break required chemical bonds, but weakly enough to allow the removal of intermediates or products. This physicochemical phenomenon is known as **Sabatier principle**.
* This principle is extremely useful for guiding the discovery of catalytic materials with improved performance, but also imposes serious constraints on design flexibilities and attainable outcomes.

#### Density of state (DOS)
* One of the primary quantities used to describe the electronic state of a material is the electronic **density of states (DOS)**;
$$
\rho(E) dE = \text{\{number of electronic states in } (E, E+dE) \text{\}}
$$
* Once a DFT calculation has been performed, the electronic DOS can be determined by integrating the resulting electronic density in k-space. As the details of the DOS are affected by the integration accuracy in k-space, it is generally recommended to use a large number of k-points to calculate the DOS. 
* The most straightforward definition of a metal is that metals are materials with a nonzero DOS at the Fermi level.
* One useful way to think about this is to consider what happens when an electric field is applied to the metal. An electric field will accelerate electrons to higher energies than the electrons have when there is no field. In a metal, there are electronic states just above the Fermi level that can be populated by these accelerated electrons, and as a result the material can readily conduct electricity.
* The DOS of the insulators have a different picture; its DOS can be divided into two separate regions, the valence band and the conduction band. The valence band is the collection of all occupied electronic states, while all states in the conduction band are unoccupied (at T = 0 K). The region of energy that separates the valence and conduction bands contains no electronic states at all; this is the band gap.
* Materials with a band gap are classified as either semiconductors if their band gap is "small", or insulators if their band gap is "wide". The distinction between these two types of materials is somewhat arbitrary, but band gaps larger than ~3 eV are typically considered as wide band gaps.

##### LDOS
* To interpret the electronic structure of a material, it is often useful to understand what states are important in the vicinity of specific atoms.
* One standard way to do this is to use the **local density of states (LDOS)**, defined as the number of electronic states at a specified energy weighted by the fraction of the total electron density for those states that appears in a specified volume around a nuclei.
* Typically, this volume is simply taken to be spherical; so to calculate the LDOS we must specify the effective radii of each atom of interest. This definition cannot be made unambiguously. If a radius that is too small is used, information on electronic states that are genuinely associated with the nuclei will be missed. If the radius is too large, the LDOS will include contributions from other atoms.

##### Bader charge
* It is often convenient to thing of atoms within bulk materials or molecules as having net charges. In an ionic material such as $\ce{NaCl}$, for example, it is convenient to associate charges of +1 and -1 (in units of electron charge) to $\ce{Na}$ and $\ce{Cl}$ atoms.
* The ambiguity in defining the volumes used for LDOS calculations illustrates why making this assignment from a calculated electron density is not necessarily a simple task.
* A widely used method within the plane-wave calculations is the **Bader decomposition**, which uses stationary points in the three-dimensional electron density to partition electrons among different atoms.

#### d-band model
* In the first-principle calculation, the analysis of the bulk or surface by the DOS is often done especially when the adsorption phenomena is important.
* To analyze the interaction between adsorbate and surface, the analysis based on **d-band model** is often carried out, which uses the electronic DOS of the surface.
* As the adsorbate approaches the surface, the adsorbate electronic state will begin to interact with the electronic state of the surface.
* It is useful to divide the electronic states of transition metal surfaces into two types: the sp-bands and the d-bands. The sp-bands originate from the metal valence s and p atomic orbitals that interact to form broad overlapping bands. The valence d orbitals of transition metals are more localized than the s and p orbitals, and they interaction more weakly and form narrower band close to the highest occupied state, the Fermi level.

* The metal-$\sigma_u^*$ interaction usually leads to an attraction because the $\sigma_u^*$-derived anti-bonding level remains unoccupied.
* It depends on the position of the Fermi level whether the overall interaction is purely attractive or not.
* For transition metals, the Fermi level lies within the d-band thus both the $\sigma_g$ and $\sigma_u^*$-derived bonding levels are occupied while the anti-bonding states are empty. This causes an attractive interaction.
* For a noble metal, both the bonding and anti-bonding state of the $\sigma_g-d$ interaction are occupied, making this interaction repulsive. This is the reason why nobel metals are noble i.e., less reactive than transition metals.
* Recently, several improvements of d-band center theory e.g. its spin-polarized version [Bhattacharjee, SciRep, 2016] or using upper d-band edge instead of the d-band center has been proposed.[Xin, PRB, 2014]

* The DOS projected onto the d-states that interact with the adsorbate state can be characterized by the moments of the d DOS. The first moment is the d-band center:
$$
\epsilon_d = \frac{\int_{-\infty}^{\infty}n_d(\epsilon)\epsilon d\epsilon}{\int_{-\infty}^{\infty}n_d(\epsilon) d\epsilon}
$$
and the higher moments (n > 1)
$$
\epsilon_d^{(n)} = \frac{\int_{-\infty}^{\infty}n_d(\epsilon)(\epsilon - \epsilon_d)^n d\epsilon}{\int_{-\infty}^{\infty}n_d(\epsilon) d\epsilon}
$$
describe the shape of the band in more detail; the width (n = 2), the skewness (n = 3), and the kurtosis (n = 4).

### 3-X. Chemical Process on Surface
#### adsorption and desorption
* There are two types of adsorptions. One is **physisorption**, mainly driven by van der Waals forces. This occurs due to an attraction between mutually induced dipoles of the electron clouds surrounding the atom or molecule and in the surface.
* Close to the surface, when the electron clouds of the adsorbate and surface begins overlapping, chemical bonds may for or broken. This stronger form of adsorption is called **chemisorption**.
* In principle, the adsorption energy can be measured; from the measured rate of desorption (**temperature-programmed desorption (TPD)**). Another way is to measure the temperature increase of a surface when it is covered by adsorbates (**calorimetry**)

* Many of the interesting things that happen with surfaces, of course, occur when other chemicals are present on them.
* As an example, let us perform some calculations to understand how H atoms bind on Cu(100) surface.
* Perhaps the simplest question we can ask is: where on the surface does H prefer to be? If we look at the Cu(100) surface, there are several sites with special symmetry that are intuitively appearing as the potential binding sites for H.
* We define a supercell containing a slab model of Cu(100), place an H atom relatively close to the top layer of the surface in some positions, and then minimize the systems's energy while allowing several of the surface layers and the H atoms positions to relax.

* How strongly do the atoms prefer to be on the surface instead of building in the gas phase? This question is usually answered by calculating the adsorption energy os species on the surface. One possible definition for this quantity is
$$
E^{atomic}_{ads} = E_{H/surf} - [E_{H(gas)} + E_{surf}]
$$
Here, the three terms on the right are the total energy of surface with H adsorbed on int, the total energy of a single H atom by itself in the gas phase, and the total energy of the bare surface. This definition is easy to use, but chemically unnatural because H atoms rarely exist by themselves for very long. A more physically meaningful quantity is the energy gained (or lost) by pulling two H atoms off the surface and forming an H2 molecule in the gas phase. This is 
$$
E_{ads} = E_{H/surf} - [1/2E_{H_2(gas)} + E_{surf}]
$$
* The surface coverage possible contributes the adsorption energy. Because we are using periodic boundary conditions, putting one adsorbate oin our supercell automatically means that each adsorbate "sees" a copy of itself in each of the neighboring supercell.
* The size of the supercell controls the distance between adsorbates. If the supercell is small, it defines a surface with a high density (or coverage) of adsorbates. If the supercell is large, surfaces with lower coverages can be defined. When there is one adsorbate per surface atom, the adsorbate layer is often said to have a coverage of 1 monolayer (ML). If these is one adsorbate per two surface atoms, the coverage is 0.5 ML and so on.

##### adsorbate-adsorbate interaction
* The adsorption energy will in general depend on the presence of other adsorbates on the surface.
* The adsorbate-adsorbate interaction for O atoms on the close-packed Pt surface is, for example, the repulsion at short distance.
* This is typically the case for strongly chemisorbed adsorbates on metal surfaces.
* There can also be attractive interactions between adsorbates.
* The attractive/repulsive interaction leads to the formation of ordered patterns.

#### diffusion
* It is generally so that diffusion barriers for simple adsorbates on metal surfaces are not much higher than 0.5 eV. As a result, diffusion is very fast at temperatures above 300 K where most industrial catalytic processes take place.
* If, however, there are high coverages of adsorbates, such that there are not free sites available to diffuse into, this picture may change.
* There is a liner relationship between adsorption energy and diffusion barrier.[Nilekar, AngewChemIntEd, 2006, 45, 7046]
* On metal oxide surface or other inorganic compound catalysts, the distance between adsorption sites are larger than that on the metal surfaces. Thus the diffusion barriers are often found to be larger.

#### surface reaction
* After the reactants are adsorbed and perhaps dissociated, they will diffuse and recombine to form new molecules before the product is desorbed in gas or liquid phase. This kind of chemical event is categorized in to the surface reaction.

## 3. Theoretical/Computational Approach for Reaction Engineering
### 3-1. Reactor and Reaction Engineering
* There are several questions that can be put forth about the operation of reactors and they can be used to form the basis of classifying and defining ideal conditions that are desireable for the proper measurements of reaction rates.
* The first question is whether the system exchange mass with its surroundings. If it does not, then the system is called a **batch reactor**. If it does, then the system is classified as a **flow reactor**.
* The second question involves the exchange of heat between the reactor and its surroundings. If there is no heat exchange, the reactor is **adiabatic**. If the reactor makes very good thermal contact with the surroundings, it can be held at a constant temperature and is thus **isothermal**.
* The third question is that the reactor at constant pressure or constant volume.
* The fourth question is whether time spent in the reactor by each volume element of fluid is the same. If it is not the same, there may exist a distribution of residence time, and the opposite extreme of a unique residence time is an exponential distribution.
* The fifth question focuses on a particular fixed volume element in the reactor and whether it changes as a function of time. If it does not, then the reactor is said to operate at a stationary state. If there are time variations, then the reactor is operating under transient conditions.

## 5. Modelling Catalytic Systems
### 5-1. Modeling Heterogeneous Systems
#### Basic
* Facets are usually denoted by their Miller indices. The most common facets are the most close packed, which for the fcc crystal structure is the (111) and the slightly more open (100) surface structure.
* Undercoordinated sites at edges and corners are often particularly important for catalysis. The same kind of sites is found at steps on the surface. For the fcc structure, a (211) step is often used to model undercoordinated sites.
* The PED for elementary surface reactions generally depends strongly on the surface structure. For example, the CO dissociation on Ni surface is strongly dependent on the facet.[Andersson, J.Catal., 2008, 255, 6]

* The number density of surface atoms in the stable planes of transition metal is the order of $\approx 10^{15} = (6 \times 10^{23})^{2/3}$ (CHECK), regardless of crystal face. This value of surface atom density is a good starting point for estimating the total number of adsorption sites or active sites on metal surfaces.
* Perhaps the most common method for measuring the number density of exposed metal atoms is selective chemisorption of a probe molecule like $\ce{H2}$, $\ce{CO}$, or $\ce{O2}$.

#### solid-state physics/chemistry
* There is no packing of hard spheres in space that creates higher density than the fcc and hcp structures.
* We defined the shape of the cell that is repeated periodically in space, the supercell, by lattice vectors ${\bf a}_1$, ${\bf a}_2$, and ${\bf a}_3$. If we solve the Schrödinger equation for this periodic system, the solution must satisfy a fundamental property known as a **Bloch's theorem**, which states that the solution can be expressed as a sum of terms with the form
$$
\phi_{\bf k}({\bf r}) = \exp(i{\bf k}\cdot{\bf r}) u_{\bf k}({\bf r})
$$
where $u_{\bf k}({\bf r})$ is periodic in space with the same periodicity as the supercell. That is, $u_{\bf k}({\bf r}) = u_{\bf k}({\bf r}+n_1{\bf a}_1+n_2{\bf a}_2+n_3{\bf a}_3)$ for anu integers $n_1$, $n_2$, and $n_3$.
* Because the functions $\exp(i{\bf k}\cdot{\bf r})$ are called plane waves, calculations based on this idea are frequently referred to as plane-wave calculations.
* The space of vectors ${\bf r}$ is called real space, and the space of vectors ${\bf k}$ is called reciprocal space (or k-space).
* A primitive cell is that it is a cell that is minimal in terms of volume but still contains all the information we need. We can define a primitive cell in reciprocal space, and it is called **Brillouin zone**.
* We need to evaluate
$$
\bar{g} = \frac{V_{\rm cell}}{(2\pi)^3}\int_{BZ}g({\bf k})d{\bf k}
$$
* The Monkhorst-Pack scheme is often used, in which all that is needed is to specify how many k-points are to be used in each direction in reciprocal space. If $M$ k-points are used in each direction, it is usual to label the calculation as using $M\times M\times M$ k-points.
* The symmetries of solid mean that the integrals in k-space do not need to be evaluated using the entire Brillouin zone. Instead, they can just be evaluated in a reduced portion of the zone that can ben be extended without approximation to fill the entire Brillouin zone using symmetry. This reduced region in k-space is called **irreducible Brillouin zone**. For very symmetricmaterials such as the perfect fcc crystal, using just the irreducible Brillouin zone reduced the numerical effort required to perform integrals in k-space.

* In a metal, the Brillouin zone can be divided into regions that are occupied and unoccupied by electrons. The surface in k-space that separates these two regions is called the **Fermi space**. This is significant complication numerically, because the functions that ar eintegrated change discontinuously from nonzero values to zero at the Fermi surface.
* One approach to this is the smearing method. The idea of these methods is to force the function being integrated to be continuous by "smearing" out the discontinuity. An example of a smearing function is the Fermi-Dirac function:
$$
f\left(\frac{k-k_0}{\sigma}\right) = \frac{1}{\exp\left[\left(\frac{k-k_0}{\sigma}\right)+1\right]}
$$
* As $\sigma \rightarrow 0$, the function approaches to a step function. Ideally, the result of the calculation should be obtained using some method that extrapolates the final result to the limit where the smearing is eliminated (i.e. $\sigma \rightarrow 0$).
* One widely used smearing method is **Methfessel-Paxton method**, in which more complicated smearing functions are used.

#### k-point tips
* Increasing the volume of supercell reduced the number of k-points needed to achieve convergence, because volume increases in real space correspond to volume decreases in reciprocal space.

#### Cutoff energy
* The solution of the Schrödinger equation for supercell have the form
$$
\phi_{\bf k}({\bf r}) = \exp(i{\bf k \cdot r})u_{\bf k}({\bf r})
$$
where $u_{\bf k}({\bf r})$ is periodic in shace with the same periodicity as the supercell. The periodicity of $u_{\bf k}({\bf r})$ means that it can be expanded in terms of a special set of plane waves;
$$
u_{\bf k}({\bf r}) = \sum_{\bf G}\exp(i{\bf G \cdot r})
$$
where the summation is all over vectors defined by ${\bf G} = m_1{\bf b}_1 + m_2{\bf b}_2 + m_3{\bf b}_3 $ with integer values for $m_i$. These set of vectors defined by ${\bf G}$ in reciprocal space are defined so that any real-space lattice vector ${\bf a}_i$ satisfies ${\bf G}{\bf a}_i = 2\pi m_i$.
* Combining two equation gives
$$
\phi_{\bf k}({\bf r}) = \sum_{\bf G}c_{\bf k+G}\exp\left[i({\bf k+G}){\bf r}\right]
$$
This corresponds to the solution of the Schrödinger equation with kinetic energy
$$
E = \frac{\hbar^2}{2m}|{\bf k+G}|^2
$$
It is reasonable to expect that solution with lower energies are more physically important than solutions with high energies. Thus, it is usual to truncate the infinite sum over $G$ to include only solutions with kinetic energy less than the **cutoff energy**
$$
E_{cut} = \frac{\hbar^2}{2m}G_{cut}^2
$$
The infinite sum then reduced to
$$
\phi_{\bf k}({\bf r}) = \sum_{{\bf |k+G|}<G_{cut}}c_{\bf k+G}\exp\left[i({\bf k+G}){\bf r}\right]
$$

#### Miller index
* Miller indices form a notation system in crystallography for directions and planes in crystal lattice. Lattice planes are determined by the three integers $(h, k, l)$ also called Miller indices.
* In a cubic lattice, these indices coincide with the inverse intercepts along the lattice vectors. Thus, $(hkl)$ simply denotes a plane that intercepts at the three lattice vectors at the points $a/h$, $b/k$, and $c/l$. If one of the indices is zero, the planes are parallel to that axis.

#### Miller index 2
* The orientation of the surface plane can be defined by stating the direction of a vector normal to the plane.
* To decide which of these many vectors to used, it is usual to specify the points at which the plane intersects the three axes of the material's primitive cell or the conventional cell. The reciprocals of these intercepts are then multiplied by a scaling factor that makes each reciprocal an integral and also makes each integer as small as possible. The resulting set of number is called the **Miller index** of the surface. For example, the plane intersects the z axis of the conventional cell at 1 (in units of lattice constant) and does not intersect the x and y axes, the reciprocals of the intersects are (1/$infty$, 1/$infty$, 1/1) and thus the surface is denoted (001).
    * If the vectors defining the surface normal requires a negative sign, that component of the Miller index is denoted with an overbar.
    * The plane is usually identified by three indices enclosed in parentheses (hkl), the vector that is normal to the plane (in cubic systems) is enclosed in square brackets: [hkl].
* When the plane intercepts the x, y, and z axes at 1, 1, and 1, the reciprocals of these intercepts are (1/1, 1/1/, 1/1) thus (111) surface. This surface is an important one because it has the highest possible density of atoms in the surface layer of any possible Miller index surface of an fcc material. Surfaces with the highest atom densities for a particular crystal structure are typically the most stable, and thus play an important role in real crystal at equilibrium.
* Note that, for bcc materials the surface with the highest density of atoms is the (110) surface.
* The low-index surfaces are important in the real world because of their stability. However, high-index surfaces have important properties as well. For example, (311) surface has regions with monoatomic steps. The atoms located at step edges have a lower coordination number than other atoms in the surface, which often leads to high reactivity. The step edge atoms play an important role in the catalytic synthesis of ammonia, for example.
* In the hexagonal close-packed (hcp) structure, spheres are arranged in a single close-packed layer to form a basal plane with each atoms surrounded by six other atoms. The next layer is added by placing spheres in alternating threefold hollows of the basal plane. If a third layer is added such that the sphere s are directly above the spheres in the basal plane, we obtain the hcp structure. If the third layer atoms are added in the hollows not directly above the basal plane, the fcc structure is obtained (**MOVE TO BULK PART?**)
* hcp materials have a sixfold symmetry axis normal to the basal plane. Using a three-axis system to defined the Miller indices for this structure is unsatisfactory, because some planes with different Miller index might be equivalent by symmetry.
* A solution to this is to use a four axis, four-index system for hcp solids. The Miller indices are formed as before by taking the reciprocals of the intercepts of the plane with the four axes. Because the four axes are not independent, the first three indices will always sum to zero.

#### slab model
* If our goal is to study a surface, our ideal model would be a slice of material that is infinite in two dimensions, but finite along the surface normal. In order to accomplish this, it may seem natural to take advantage of periodic boundary conditions in two dimensions, but not the third.
* This is achieved by a slab model. In this model, the empty space has been left above the atoms in the top portion of the supercell, and the supercell is repeated in all three dimensions. It defines a series of stacked slabs of solid material separated by empty spaces.
* The empty space separating periodic images of the slab along the z-direction is called the vacuum space. It is important when using such a model that there is enough vacuum space so that the electron density of the material tails off to zero in the vacuum and the top of one slab has essentially no effect on the bottom of the next.
* The number of atoms mimicking the surface layer is arbitrary. So, how many layers are enough? Typically, more layers are better, but using more layers also inevitably means using more computational time. The same thing is true of the vacuum layer thickness. There questions can be answered by carrying out calculations to see how some property (such as surface energy or adsorption energy) varies as teh number of layers increases. The number of layers and the vacuum layers thickness would be a compromise between computational cost and physical accuracy.
* The fact that the long dimension in the supercell includes a vacuum region leads some computational advantages: The "long" dimension in real space means the "short" dimension in the reciprocal space. Thus, if the vacuum region is large enough, the accurate results are possible using just one k-point in the surface-normal direction. As a result, it is usual in slab calculations to use an MxNx1 k-point mesh, where M and N are chosen to adequately sample k^space in the plane of the surface.

#### surface relaxation
* It is natural to expect that the spacings between layers near the surface might be somewhat different from those in the bulk. This phenomenon is called surface relaxation. Surface relaxation implies that the relaxed surface has a lower energy than the original surface. We imagine the bottom of the slab as representing the bulk part of the material, and constrain the atoms in the bottom layers in their ideal, bulk positions.
* The calculation then involves a minimization of the total energy of the supercell as a function of the positions of the atoms, with only the atoms in the top layers are allowed to move.
* Often the surface relaxation leads to a decrease in the distance between the first and second atomic layers.

#### surface energy
* The surface energy, $\sigma$, is the energy needed to cleave the bulk crystal. This quantity can be defined from a slab calculation using
$$
\sigma = \frac{1}{A}\left[ E_{slab} - n E_{bulk} \right]
$$
where $E_{slab}$ is the total energy of the slab model for the surface, $E_{bulk}$ is the energy of one atom or formula unit of the material in the bulk, $n$ is the number of atoms or formula units in the slab model, and $A$ is the total area of the surface in the slab model.
* Surface energy is typically expressed in units of Joules / meters squared ($J/m^2$) in macroscopically, while it is expressed as electron volt per angstrom ($eV/\r{A}$) in calculation (note that $ 1 J/m^2 = 16.02 eV/\r{A}^2 $)

#### Wulff construction
* The equilibrium shape of finite mesoscopic crystals can be directly derived from the surface energies by the so-called **Wulff construction** which is based on the concept that the crystal seeks to minimize its total surface energy subject to the constraint of fixed volume.

#### metal-support interaction

#### Volcano curve
* For the adsorption of reactant, the Bronsted-Evans-Polanyi (BEP) relationship is well-established for the dissociative adsorption of many diatomic molecules, such as CO, NO and N2.[Chem.Rev, 110, 2005, 2009; JCatal, 191, 301, 2000; JCP, 114, 8244, 2001; JCatal, 197, 229, 2001; JPC.C, 112, 1308 2008]
* This fact means that the barriers of dissociative adsorption are correlated to the enthalpy change across the periodic table. For instance, the barrier of CO dissociative adsorption is related to the stability of C and O; the more stable the C and O are, the lower the barrier is fo CO dissociative adsorption.
* Regarding the desorption step, the desorption steps of carbon, nitrogen and oxygen forming CH4, NH3 and H2O on many metal surfaces were calculated, and a linear relationship was found between the effective desorption barrier and the enthalpy change.
* In other words, the more stable the intermediate is, the higher the barrier is for the desorption step.
* Based on the result mention above, a heterogeneous catalytic reaction can be approximated to be two steps: adsorption and desorption, the barriers of which are both related to the stability of the key intermediate.
* Thus, the reaction rate can be represented as a function of enthalpy change. By solving the kinetic equation of the two-step model,[JPC.C, 112, 1308 2008] the plot of reaction rate in terms of turnover frequency (TOF) against enthalpy change is shown in Figure X.
* If the adsorption step is assumed to be rate-determining, the green curve in Fig.X will be found between the reaction rate and the enthalpy change, while the blue curve is obtained when the desorption step is rate-determining.
* By taking both curves into account, a volcano-shaped curve will be found between TOF and deltaH.
* For catalysts binding strongly with adsorbates, the enthalpy change will be too negative, resulting in a high barrier for the desorption step. Thus, the desorption will be rate-determining, and the activity of this type of catalysts can be found on the blue curve in Fig.X
* On the other hand, a low adsorption energy between the surface and adsorbates will lead to a high barrier for the adsorption step. Therefore, the adsorption will limit the overall activity, which is the case on the green curve.
* In summary, based on the BEP relation and the two-step model, a volcano-shaped curve can be obtained between the stability of key intermediate and the overall activity of the catalysts.
* A good catalyst should have the proper adsorption energy of key intermediate.
* Therefore, adsorption energy may be a good indicate for rational catalyst design.

* The above model assumed that only one reactant and one product involved in the catalytic reaction. This two-step model is a reasonable simplification for many reactions in which only one adsorption or desorption is dominant in the overall activity.
* For example in the NH3 synthesis, the dissociative adsorption on N2 and the desorption of NH3 are mainly affects the overall activity.
* However, in some complex surface reaction, more than one reactants or products affect the overall reaction rate. This leads to the three-dimensional volcano surface for the catalyst activity.

#### Scaling relationship
* the presence of linear correlations between adsorption energies of similar adsorbates, known as **scaling relations**, reduces the complexity of DFT-based catalytic models and facilitates the simultaneous analysis of numerous materials through Sabatier-type activity plots.

##### Scaling among adsorption energy
* It is known that adsorption energies of different surface intermediates that bind to the surface through the same atom scale with each other. Abild-Pederen et al. have found that, for example, the scaling relationship of CHn (n = 1, 2, 3) adsorption energy; there is a linear relationship between the C atom adsorption energy and the CHn adsorption energy over a number of metal surfaces (Fig.X).[Abild-Pedersen, PRL, 99, 016105, 2007]
* It is evident from the earlier discussion that scaling among adsorption energies should not be limited to transition metal surfaces. In fact, even for metal-terminated surfaces of more complex systems like transition metal compounds (oxides, nitrides, sulfides, and carbides), where there is mixed covalent, ionic bonding between the surface cations and anions, there is scaling between electronically similar adsorbates.[Fernandez, AngewChemIntEd, 47, 4683, 2008]

##### Scaling among $E_a$ and $\Delta E$
* It is not surprising that transition state energies also correlate with adsorption energies. As for the relationships between adsorption energies, scaling between adsorption energies and transition-state energies is extremely important in building an understanding of heterogeneous catalysis.
* They provide guidance in building kinetic models to understand trends in catalytic activity.
* Let $E^{TS}$ be a set of energies describing the energy needed to move between two minima on potential energy surface for a set of different catalysts. Furthermore, let $\Delta E_i$ be a set of adsorption energies, relevant ofr the process of moving between two minima. We can now define a functional form of $E^{TS}(\Delta E_i)$, which is a map from the space of adsorption energies to the space of transition state energies.
* To first order in $\Delta E_i$, $E^{TS}$ will be give as a linear combination of $\Delta E_i$:
$$
E^{TS} = \sum_i\gamma_i \Delta E_i + \xi
$$
* The set of functions defined above constitute a class of linear relations - the linear transition state energy scaling relations.
* The transition-state scaling relations imply the scaling relation for the activation energy $E_a$ of surface chemical reaction.
* Linear correlations between activation (free) energies and reaction (free) energies is a well-established approach in the understanding of trends in chemical relationship that dates back to Bronstead in 1928 and Evans and Polanyi 10 years later.
* We note that in principle there is a different line for every surface geometry, sone one should think of a family of transition state caling lines. For example in the NO dissociation, a large number of geometries have been investigated, and the lines for the close-packed and the stepped surface basically define the upper and lower bound to the scaling lines.

* It turns out that if one compares dissociation of a number of similar molecules, their transition state energy scale with the dissociative chemisorption energy in much the same way (Fig.X). This is a remarkable result indicating that the nature of the relationship between the transition state and the final state are quite similar among these molecules.
* The universal relationship for the close-packed surfaces is found to be
$$
E^{TS} = (0.90\pm0.04)\Delta E_{\rm diss} + (2.07\pm0.07)\ {\rm eV}
$$
and the relationship for the stepped surface is 
$$
E^{TS} = (0.87\pm0.05)\Delta E_{\rm diss} + (1.34\pm0.09)\ {\rm eV}
$$
* The slopes of the relations are similar to each other, indicating the transition states of the reactants considered are final-state like.
* The intercepts are different, and this difference identifies the structure dependence of the relationship; the stepped surfaces have barriers that are much smaller than on the closed-packed surfaces.[JCatal, 209, 275, 2002]
* Note that there are some exceptions to these scalings, especially when molecules with weak interatomic bonds are considered suc as in the dissociation of H2.

#### Ab initio thermodynamics
* Ab initio atomistic simulations are, under the Born-Oppenheimer approximation, useful tool to evaluate the total energy (i.e. electronic plus nuclear repulsion) from atomic configurations. This total energy is, however, a zero-temperature and zero-pressure quantity. Since catalytic reactions often take place at high-temperature and/or high-pressure condition, this temperature or pressure gap should be filled somehow. This is indeed possible, when the thermodynamic potentials like Gibbs free energies are evaluated from ab initio atomistic simulations.
* In this section, such approach called **ab initio (atomistic) thermodynamics (AITD)** is introduced.

* Under the condition of fixed temperature $T$ and pressure $p$, the appropriate thermodynamics potential to consider is the Gibbs free energy $G(T,p)$.
* The oxidation of metal $M$ surface is considered as an example here, because it is one of the most common application of the AITD technique.

* The stability of a particular surface is given by its surface free energy per unit area
$$
\gamma = \frac{1}{A}\left[ G(T, p_i, N_{\rm M}, N_{\rm O}) - N_{\rm M}g_{\rm M}(T, p) - N_{\rm O} \mu_{\rm O}(T, p) \right]
$$
where $\mu_{\rm O}$ is the chemical potential of O atom. The number of metal and O atoms are denoted as $N_{\rm M}$ and $N_{\rm O}$, respectively.$g_{\rm M}$ is the Gibbs free energy per metal atom, which can be replaced by the total energy of bulk or surface divided by $N_{\rm M}$. 

* When discussing the stability of phases that result from adsorbing species at surface, it is convenient to choose some reference surface. One can introduce the surface free energy of the clean surface as
$$
\gamma^{\rm clean}(T,p) = \frac{1}{A}\left[G(T,p,0,N_{M}^{\rm clean})-N_{\rm M}^{\rm clean}E_{\rm M}^{\rm total}\right]
$$
and the Gibbs free energy of adsorption can be measured from the clean surface as
$$
\begin{align*}
\Delta G^{\rm ad}(T,p) 
&= \gamma(T,p,N_{\rm O},N_{\rm M}) - \gamma(T,p,0,N_{\rm M}^{\rm clean}) \\
&= \frac{1}{A}[ G(T,p,N_{\rm O},N_{\rm M}) - G(T,p,0,N_{\rm M}^{\rm clean}) \\
&- N_{\rm O}\mu_{\rm O}(T,p) - (N_{\rm M}-N_{\rm M}^{\rm clean})E_{\rm M}^{\rm total} ]
\end{align*}
$$
* The last term arises only when clean and adsorbed surfaces have different number of metal atoms.
* We are interested in describing the thermodynamic stability of a metal, M, in equilibrium with gas-phase $\ce{O2}$ at pressure $P_{\ce{O2}}$ and temperature $T$. We assume that we already know a series of candidate crystal structures for the metal and its oxide. As an example, we will examine the copper and copper oxide surfaces, and thus crystal structures of $\ce{Cu}$ and $\ce{Cu2O}$ are assumed to be known from experiments.
* Thermodynamically, we would like to know which material minimizes the free energy of a system containing gaseous $\ce{O2}$ and a solid at the specified conditions. A useful way to do this is to define the grand (canonical?) potential associated with each crystal structure. The grand potential for a metal oxide containing $N_M$ metal atoms and $N_O$ oxygen atoms is defined by
$$
\Omega(T, \mu_O, \mu_M) = E(N_M, N_O) - TS - \mu_O N_O - \mu_M N_M
$$
where $E$ is the internal energy of the metal oxide, $S$ is the material's entropy, and $\mu_X$ is the chemical potential of atomic species $X$.
* If we have a series of different materials ($i = 1, 2, \cdots$), then we can use this expression to define the grand potential of each material, $\Omega_i$.
* We can interpret the internal energy in the grand potential as simply the total energy from a DFT calculation for the material.
* It is then sensible to compare the grand potentials of the different materials by normalizing the DFT energies so that every DFT calculation describes a material with the same number of metals atoms. If we do this, $\Omega_i$ can be rewritten as
$$
\Omega_i(T, \mu_O) = E_i - TS_i - \mu_O N_{O,i} - \Omega^M
$$
where $\Omega^M$ is an additive constant that is the same for every material.

* For ideal gases, the chemical potential can be rigorously derived from statistical mechanism. A useful definition of the chemical potential for $\ce{O2}$ is
$$
\mu_{\ce{O2}} = E_{\ce{O2}}^{\rm total} + E_{\ce{O2}}^{\rm ZPE} + \Delta \mu_{\ce{O2}}(T, P^o) + kT\ln(P/P^o)
$$
Here, $E_{\ce{O2}}^{\rm total}$ and $E_{\ce{O2}}^{\rm ZPE}$ are the total energy of an isolated $\ce{O2}$ molecule at T = 0 K and the zero-point energy (O2). We can obtain them from a simple DFT calculation.
* $\tilde{\mu_{\ce{O2}}}$ is the difference in the chemical potential of $\ce{O2}$ between T = 0 K and the temperature of interest (at the reference pressure). This chemical potential difference (and also ZPEs) for standard gases species are tabulated in the *NIST-JANAF Thermodynamical Tables*.

* Now we have derived the chemical potential of oxygen molecule. Since the chemical potential of molecular oxygen and atomic oxygen are related by
$$
\mu_{\ce{O}} = \frac{1}{2}\mu_{\ce{O2}}
$$
, we can calculated $\mu_{O}$ as
$$
\mu_O(T, p) = \frac{1}{2}E_{\ce{O2}}^{\rm total} + \frac{1}{2}E_{\ce{O2}}^{\rm ZPE} + \Delta \mu_O(T) + \frac{1}{2}k_B T \ln\left(\frac{p}{p^0}\right)
$$

##### Pd(100) case
* The Gibbs free energy change by adsorption is
$$
\begin{align*}
\Delta G^{\rm ad}(T, p) &= \frac{1}{A}[ E^{\rm total}(N_{\rm O}, N_{\rm M}) - E^{\rm total}(0, N_{\rm M}^{\rm clean}) \\
&-(N_{\rm M}-N_{\rm M}^{\rm clean})E_{\rm M}^{\rm total}
-\frac{N_{\rm O}}{2}E_{\ce{O2}}^{\rm total}
-N_{\rm O}\Delta\mu_{\rm O}(T, p) ]
\end{align*}
$$
* This above equation allows to directly plot the Gibbs free energy of adsorption for each surface model as a function of $\Delta \mu_{\rm O}$; see Figure X.
* This yields a straight line for each model, and at any given $\Delta \mu_{\rm O}$ the model with the lowest lying line is identified as the most stable as it has the lowest $\Delta G^{\rm ad}$.
* If one concentrate only on the most stable structures, it is possible to convert range of chemical potential into (T, p) ranges, and plot these stability ranges in surface phase diagram (Fig.X)

##### Ag
* A good application of the ab initio thermodynamics to catalysis is the ethylene epoxidation by Ag. Ag is a well-known commercial catalyst for epoxidation of ethylene, where the reaction takes place with an $\ce{O2}$ pressure of approximately atmospheric pressure and temperature of 200-300 $^oC$.
* If the reaction is attempted using clean Ag under ultra-high-vacuum conditions, the surface is found to be effectively inactive as a catalyst. Therefore, the surface of active Ag catalyst under industrial conditions is not mettalics Ag.
* An example of theoretical phase diagram for oxygen interacting with the Ag(111) surface was given by Li, Stampfk, and Scheffler.[Ref]
* It can be seen from Fig.X that, at most temperatures, there is a range of pressures spanning several orders of magnitude for which the surface oxide structure is more stable than either the bulk oxide or the clean metal surface.
* This phase diagram strongly suggests that the working catalyst under industrial conditions is a surface oxide rather than bare Ag.

##### two-dimensional
* As a matter of fact, the AITD technique is not limited to the oxidation of metal surface. Moreover, one can extend above formulation when gas phase has two components.
* One good example for this is the CO oxidation reaction on RuO2(110) surface studies by Reuter and Scheffler.[PRL, 90, 046103, 2003] In this case, CO and O2 are present in gas phase and CO oxidation to CO2 take place with the surface reaction since gas-phase CO oxidation is quite slow. CO molecule and O atom (as a result of O2 dissociative adsorption) may take several adsorption positions on the RuO2(110) surface, while in this case two Ru sites (Ru bridge, ${\rm Ru_{br}}$ or Ru coordinatively unsaturated site ${\rm Ru_{cus}}$) are possible.
* To identify thermodynamically stable surfaces, AITD diagram would be helpful. Since the surface energy of RuO2-CO-O system is calculated likewise the Eq.X as follows;
$$
\begin{align*}
\gamma(T, p_{\ce{CO}}, p_{\ce{O2}}) =
& \frac{1}{A}[G(T, p_{\ce{CO}}, p_{\ce{O2}})-G_{\rm clean}(T) \\
&-N_{\rm CO}\mu_{\rm CO}(T, p_{\ce{CO}})-N_{\rm O}\mu_{\rm O}(T, p_{\ce{O2}})]
\end{align*}
$$
* The AITD diagram in this case should be three-dimensional ($\mu_{\rm CO}, \mu_{\rm O}, \gamma$) but when only the lowest-energy surface structures are plotted it becomes a two-dimensional plot; see Fig.X for the result.
* In the figure, four stable surface structures (${\rm O^{br}/-}$, ${\rm O^{br}/O^{cus}}$, ${\rm O^{br}/CO^{cus}}$, ${\rm CO^{br}/CO^{cus}}$) are shown. Among them, ${\rm O^{br}/CO^{cus}}$ surface is directly related to the function for the catalytic CO oxidation on RuO2 surface, because both CO and O are co-exist in the proximity on the surface. On the other hand, other three cases the surface is covered with only one reactant or adsorption of CO is weak (${\rm O^{br}/-}$). Therefore, we can expect the higher catalytic activity for CO oxidation in $(\mu_{\rm CO}, \mu_{\rm O})$ regions corresponding to ${\rm O^{br}/CO^{cus}}$.

### 5-2. Modeling Homogeneous Systems
#### Solvation models
* An important aspect of computational chemistry is to evaluate the effect of the environment, such as a solvent.
* Methods for evaluating the solvent effect may broadly be divided into two types: those describing the individual solvent molecules and those that treat the solvent as a continuous medium.
* Combinations are also possible, for example by explicitly considering the first solvation shell and treating the rest by a continuum model.
* The effect of solvation can be partitioned into the non-specific solvation and the specific solvation.
* The non-specific effets are primarily solvent polarization and orientation of the solvent electric multipole moments by the solute, where the dipole interaction is usually the most important.
* These effects causes a screening of charge interactions, leading to the (macroscopic) dielectric constant being larger than 1.
* The microscopic interactions are primarily located in the first solvation shell, although the second solvation shell may also be important for mutiple-charged ions. The microscopic interactions depend on the specific nature of the solvent molecule, such as the shape and ability to form hydrogen bonds.

##### Continuum solvation models
* Continuum models consider the solvent as a uniform polarizable medium with a dielectric constant of $\epsilon$ and with the solute $M$ placed in a suitable shaped hole (or cavity) in the medium (Figure 15.9).
* Creation of a hole in the medium costs energy i.e. a destabilization, while dispersion interaction between the solvent and solute add a stabilization (this is rougly the van der Waals energy between the solvent and solute).
* The electric charge distribution of $M$ will polarize the medium (induce charge moments), which in turn acts back on the molecule, thereby producing an electrostatic stabilization.
* One may also consider an exchange-repulsive component.
* The solvation (free) energy may thus be written as a sum of components as shown below:
$$
\Delta G_{\rm solvation} = \Delta G_{\rm elec} + \Delta G_{\rm cavity} + \Delta G_{\rm dispersion} + \Delta G_{\rm exchange} + \cdots
$$
* The dielectric medium is normally taken to have a constant value of $\epsilon$, but may for some purposes also be taken to depend, for example, on the distance from $M$.
* For dynamical phenomena, it can also be allowed to be frequency dependent, that is the response of the solvent is different for a "fast" process, such as an electronic stransition, and a "slow" process, such as the reorientation of solvent molecules.

* The simplest shape for the cavity is a sphere or an ellipsoid. This has the advantage that the electrostatic interaction between $M$ and the dielectric medium may be calculated analytically.
* More realistic models employ molecular-shaped cavities, generated by e.g. interlocking shperes located on each nucleus.
* Taking the atomic radius mutiplied a suitable factor (typically ~1.2) forms a van der Waals surface.
* Such a surface may have small "pockets" where no solvent molecules can enter and a more appropriate descriptor may be defined as the surface traced out by a spherical particle of a given radius (which depends on a solvent molecule) rolling on the van der Waals surface. This is denoted as *Solvent Accessible Surface (SAS)* as is illustrated in Figure 5.10.
* Since an SAS is computationally more expensive to generate than a van der Waals surface, and since the difference is small, a van der Waals surface is often used in practice.
* An additional disadvantage in using a SAS is that a very small displacement of an atom may alter the SAS in a discontinuous fashion, as a "pocket" in the van der Waals surface suddenly becomes too small to allow the solvent probe to enter.
* Whichever the SAS or the van der Waals surface is used, it is generally found that the shape of the cavity is important and these molecular-shaped cavities are necessary to be able to obtain good agreement with experimental data, such as solvation energies.
* The energy required in creating the cavity (entropy factors and loss of solvent-solvent interactions), the stabilization due to dispersion between the solute and solvent, as well as the exchange repulsion energy is usually assumed to be proportional to the surface area.
* The corresponding energy terms are often parametrized as being proportional to the total cavity area (a single proportionality constant) or parametrized by having a constant $\xi$ specific for each atom type (analogous to van der Waals parameters in force field methods), with the $\xi$ parameters being determined by fitting to experimental solvation data:
$$
\Delta G_{\rm cavity} + \Delta G_{\rm dispersion} + \Delta G_{\rm exchange} = \alpha {\rm SAS} + \beta \\
\Delta G_{\rm cavity} + \Delta G_{\rm dispersion} + \Delta G_{\rm exchange} = \sum_i^{N_{atoms}}\xi_i S_i
$$
* The electrostatic component of Eq.(15.65) can be described at several different levels of approximation, as discussed in the following sections.

##### Possion equation
* The *Poisson equation* is a second-order differential equation describing the connection between the electrostatic potential $\phi$, the charge distribution $\rho$ and the dielectric constant $\epsilon$:
$$
\nabla\cdot\left[ \epsilon({\bf r})\nabla \phi({\bf r}) \right] = -4\pi\rho({\bf r})
$$
* Note that the dielectric "constant" may depend on the position, but when its dependence is neglected the above equation becomes
$$
\nabla^2 \phi({\bf r}) = -\frac{4\pi}{\epsilon}\rho({\bf r})
$$
* If the charge distribution is point charge, the solution of Eq.(15.69) reduces to the Coulomb interaction.
* Eq.(15.68) can be used to describe, for example, the solvation of protein in water, where the protein region is taken to have a low dielectric constant ($2 < \epsilon < 5$) while the solvent has a high dielectric consant ($\epsilon = 78$ if water).
* The Poisson equation is a differential equation that must be solved numerically, typically by a grid representation, and the results give information about the electrostatic potential at any point in space.
* It can be mapped on to the surface of the solute, where it may suggest regions for inteaction with other polar molecules. It can also be used for generating the *reaction field*, defined as the difference between the potential in the presence of a solvent and in vacuum; $\phi_{\rm reac} = \phi_{\rm solv} - \phi_{\rm vac}$.Multiplication of the reaction field with the solute charges in either a continuous charge ($\rho$) or partial charge ($Q$) description gives the electrostatic component of the free energy:
$$
\Delta G_{\rm elec} = 
\left\{
\begin{align*}
&\frac{1}{2}\int \rho({\bf r})\phi_{\rm reac}({\bf r}) d{\bf r} & &(\text{continuous charge}) \\
&\frac{1}{2}\sum_i Q_i({\bf r}_i)\phi_{\rm reac}({\bf r}_i) & &(\text{partial charge})
\end{align*}
\right.
$$

##### Born/Onsager/Kirkwood models
* The numerical aspects of solving the Poisson or Poisson-Boltzmann equations make them too demanding for use in connection with, for example, geometry optimizations or simulations of macromolecules.
* For certain special cases, however, the Poisson equation (Eq.(15.68)) can be solved analytically, and this forms the basis for many approximate models for estimating the electronic component in Eq.(15.65).
* The simplest reaction field model is a spherical cavity, where only the lowest-order electric moment of the molecule is taken into account.
* For a net charge $q$ in a cavity with radius $a$, the difference in energy between a vacuum and a medium with a dielectric constant $\epsilon$ is given by the *Born model*.
$$
\Delta G_{\rm elec}(q) = -\left(1-\frac{1}{\epsilon}\right)\frac{q^2}{2a}
$$
* It can be noted that the Born model predicts equal solvation energies for positive and negative ions of the same size, which is not the observed behavior in solvent such as water.
* The reciprocal dependence on the dielectric constant furthermore means that the calculated solvent effect is sensitive to the variation of $\epsilon$ in the low dielectric limit but is virtually unaffected in the high dielectric limit.

* Using partial atomic charges in Eq.(15.76) is iften caleld the *Generalized Born model*. In this method, the Coulomb interaction between the partial charges $Q$ is combined with the Born formula by means of a function $f_{ij}$ depending on the internuclear distance and Born radii for each of the two atoms, $a_i$ and $a_j$.
$$
\Delta G_{\rm elec}(Q_i, Q_j) = -\left(1-\frac{1}{\epsilon}\right)\frac{Q_i Q_j}{f_{ij}} \\
f_{ij} = \sqrt{r_{ij}^2 + a_{ij}^2\exp(-D)} \\
a_{ij}^2 = a_i a_j, \ \  D=\frac{r_{ij}^2}{4a_{ij}^2}
$$
* The effective Born radius for a given atom depends on the nature and position of all the atoms. The dependence on the other atoms is in practive relatively weak, and updates of the $a_i$ parameters can be done at suitable intervals, for example when updating the non-bonded list in an optimization or simulation.
* The boundary between the solute and solvent is usually taken as a modified van der Waals surface generated from the unification of atomic van der Waals radii scaled by a suitable factor.
* The cavity/dispersion terms are parametrized according to the SAS, as in Eq.(15.67).

* The dipole in a spherical cavity is known as the *Onsager model*, which for a dipole moment $\mu$ leads to an energy stabilization given by
$$
\Delta G_{\rm elec}(\mu) = -\frac{\epsilon-1}{2\epsilon+1}\frac{\mu^2}{a^3}
$$
* The *Kirkwood model* refers to a general multipole expansion in a spherical cavity.

* The charge distribution of the molecule can be represented either as atom-centerd partial harges or as a multipole expansion. The lowest-order approximation for a neutral molecule considers only the dipole moment.
* This may be a quite poor approximation and fails completely for symmetric molcules that do not have a dipole moment.
* It is often necessary to extend the expansion up to order six or more, in order to obtained converged results, that is including dipole, quadrupole, octapole, etc., moments.
* Furthermore, only for small and symmetric molecules can be approximation of a spherical or ellipsoidal cavity be considered realistic. The use of the Born/Onsager/Kirkwood models should therefore only be considered as a rough estimate of the solvent effects, and quantitative results can rarely be obtained.

##### Self-consistent reaction field models
* A classical description of the molecule $M$ in Figure 15.9 can be a force field with (partial) atomic charges, while a quantum chemical description is also possible and it involves calculation of the electronic wave function.
* When a quantum description of $M$ is employed, the calculated electric moments induce charges in the dielectric medium, which in turn acts back on the molecule, causing the wave function to respond and thereby changing the electric moments, etc.
* The interaction with the solvent model must thus be calcualted by an iterative procedure, leading to various *self-consistent reaction field (SCRF)* models.
* For spherical or ellipsoidal cavities, the Poisson equation can be solved analytically, but for molecular-shaped surfaces it must be done numerically.
* This is typically done by reformulating it in terms of a surface integral over surface charges and solving this numerically by dividing the surface into small fractions called *tesserae*, wach having an associated charge $\sigma({\bf r}_s)$.
* The surface charges are related to the electric field ${\bf F}$ (the derivative of the potential) perpendicular to the surface by
$$
4\pi\epsilon\sigma({\bf r}_s) = (\epsilon-1){\bf F}({\bf r}_s)
$$
* Once $\sigma({\bf r}_s)$ is determined, the associated potential is added as an extra term to the Hamiltonian operator:
$$
\phi_\sigma({\bf r}) = \int \frac{\sigma({\bf r}_s)}{|{\bf r}-{\bf r}_s|}d{\bf r}_s \\
H = H_0 + \phi_\sigma
$$
* The potential $\phi_\sigma$ from the surface charge is given by the molecular charge distribution (Eq.(15.81)), but also enters the Hamiltonain and this influences the molecular wave function. The procedure is therefore iterative.
* FOr the case of the Onsager model (spherical cavity, dipole moment only), the term added to the Hamiltonian is given by
$$
\phi_\sigma = -{\bf r}\cdot{\bf R}
$$
* Here, ${\bf r}$ is the dipole moment operator (i.e. the position vector) and ${\bf R}$ is proportional to the molecular dipole moment, with the proportionality constant depending on the radius of the cavity and the dielectric constant:
$$
{\bf R} = g \boldsymbol\mu \\
g = \frac{2(\epsilon-1)}{(2\epsilon+1)a^3}
$$
* The $\phi_\sigma$ operator at the Hatree-Fock level of theory corresponds to the addtion of an extra term to the Fock matrix element
$$
F_{\alpha\beta} = \braket{\chi_\alpha|{\bf F}|\chi_\beta} - g\boldsymbol\mu\braket{\chi_\alpha|{\bf r}|\chi_\beta}
$$
* The addtional integrals are just expectation values of $x$, $y$, and $z$ coordinates, and their inclusion requires very little additional computationa effort. Generalization to higher-order multipole is straightforward.
* The cavity size in Born/Onsager/Kirkwood models strongly influences the calculated stabilization, but there is no consensus on how to choose the cavity radius.
* More sophisticated models employ the molecular-shaped cavities, but there is again no consensus on the exact procedure.
* The cavity is often defined based on the van der Waals radii of the atoms in the molecule, multiplied by an empirical factor e.g. 1.2. The molecular volume may alternatively be calculated directly from the electronic wave function, for example by using an isodensity surface corresponding to a value of $10^{-3} - 10^{-4}$.
* The *Polarizable Continuum Model (PCM)* employs a van der Waals cavity formed by interlocking atomic van der Waals radii scaled by an empirical factor, a detailed description of the electrostatic potential, and parametrizes the cavity/dispersion contributions based on the surface area. Several different implementations have been published, of which the integral equation formalism PCM (IEFPCM) is the most general.
* The *Conductor-like Screening Model (COSMO)* also employs molecular-shaped cavities and represents the electrostatic potential by partial atomic charges. COSMO may be considered as a limiting case of the PCM model, where the dielectric constant is set to infinity. The *COSMO-RS (real solvent)* includes additional terms in order to model, for example, hydrogen bonding in terms of the surface charges.
* The *Solvent Models x (SMx, x being a version number)* developed by Cramer and Truhlar are generalized Born-type model, where the partial atomic charges are calculated from a wave function and the cavity/dispersion terms in Eq.(15.65) are parametrized based on the solvent exposed surface area (Eq.(15.67)). The version number of these models reflects increasingly sophisticated parametrizations.

* It should be noted that the parametrization of continuum solvent models, such as, for example the cavity size defined by the atomic radii, is against experimental free energies, which implicitly include entropy and finite temperature effects.
* These effects should therefore not be added from calculated frequency, as this effectively would be a double counting.
* Furthermore, the use of the ideal-gas rigid-rotor harmonic oscillator approximation for calculated finite temperature effects is unlikely to be a good approximation for condensed phases.

* The "mixed" solvent models, where the first solvation shell is accounted for by including a number of solvent molecules, implicitly include the solute-solvent cavity/dispersion terms, although the corresponding terms between the solvent molecules and the continuum are usually neglected.
* Once discrete solvent molecules are included, however, the problem of configuration sampling arises. Furthermore, a parametrization of the continuum model against experimental data must be done by explicitly taking the first solvation shell into account.
* Nevertheless, the first solvation shell is in many cases by far the most important, and mixed model may yield substantially better results than pure continuum model, at the price of an increased computational effort.

### 5-3. Comparing Experiment and Computational Results
##### TPD
* Temperature-programmed desorption (TPD)
* TPD spectrum can be simulated with microkinetic analysis, by using the first-order Polanyi-Wigner equation
$$
r(t) = \frac{d\theta}{dt} = \nu\theta\exp\left(-\frac{G_a}{RT}\right)
$$
where $r(t), \theta, \nu, G_a$ are the reaction (desorption) rate at time $t$, surface coverage, pre-exponential factor, and desorption free energy.

## 6. Computational Tools for Catalytic Chemistry
### 6-1. Tools in Electronic Structure Theory
#### VASP
#### Gaussian

### 6-2. Tools in Chemical Kinetics
#### Chemkin

### 6-3. Tools in Chemical Engineering
#### COMSOL

#### OpenFOAM

#### ANSYS
