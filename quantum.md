# Introduction
* Most of the chemical reactions employed in the chemical industry are performed in the presence of a catalyst. This is because catalysts are useful in promoting the chemical reaction, and consequently higher product yield can be expected.
* Another popular usage of catalysts is to convert hazardous waste into less harmful products. This is often related to the automobiles, as catalysts are used to convert the harmful nitrogen oxides (NOx), carbon monooxides (CO), and unburnt hydrocarbons in car exhaust.
* In addition, there is a growing consensus that the world's increased demand for fuels and base chemicals will need to be met by so-called "carbon-neutral" technologies. For example, researchers are trying to convert the carbon dioxides (CO2) into the useful hydrocarbon compounds. Catalysts are of great importance for this purpose, meaning the catalytic science is definitely a key technology for the global environmental issues.

* Traditionally, the field for catalytic science is divided into three areas: heterogeneous, homogeneous, and enzyme catalysis.
* Heterogeneous catalysis are present in a phase that is different from that of the reactants: typically, the reactants are in gas or liquid phase but the catalyst is a solid material. Thus the reactions mainly take place at the surface of the solid material.
* Homogeneous catalysts operate in the same phase as the reactants. This category involves many important reactions in organic chemistry, where the molecular catalyst such as transition metal complexes work as catalyst and convert the reactant species into product species. As a result the homogeneous catalysis plays a significant role in drug synthesis or polymer science.
* Enzyme catalysts are specialized proteins. The chemically active part of enzymes is often a tiny part of the protein, and enzyme catalysis can be viewed as a special kind of heterogeneous catalysis. Obviously, this field is important in biochemistry or pharmacology.

# I. Basics in Theoretical and Computational Methods
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

#### meta-GGA
* The logical extension of GGA method is to allow the EX functional to include the variables of higher order derivatives of the electron density, for example the Laplacian ($\nabla^2 \rho$) as the second-order derivative term.
* The same information is actually gained by introducing the orbital kinetic energy density $\tau$, which is
$$
\tau({\bf r}) = \frac{1}{2}\sum_i^{occ}|\nabla \phi_i({\bf r})|^2
$$
* and this approach is more often used.
* Including of either the Laplacian or the orbital kinetic energy density $\tau$ as a variable leads to so-called *meta-GGA functionals*.

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

#### geometry optimization
* 

##### local optimization
* 

##### global optimization
* 

#### vibrational analysis
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

#### transition state search
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

##### Dimer method
* 

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
\begin{cases}
0.5\left[\cos\left(\frac{\pi R_{ij}}{R_c}\right)+1\right] & (R_{ij} \leq R_c) \\
0 & (R_{ij} \gt R_c)
\end{cases}
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

## 2. Chemical Kinetics and Catalytic Activity

<!-- Jensen begin -->
### Transition state theory
* Consider a chemical reaction of the type $A + B \rightarrow C + D$. The rate of reaction may be written as
$$
\frac{d[C]}{dt} = \frac{d[D]}{dt} = -\frac{d[A]}{dt} = -\frac{d[B]}{dt} = k[A][B]
$$
* Here $k$ being the rate constant.
* If $k$ is known, the concentration of the various species can be calculated at any given time from the initial concentrations. At the microscopic level, the rate constant is a function of the quantum states of $A$, $B$, $C$ and $D$, that is the translational, rotational, vibrational, and electronic quantum numbers.
* The macroscopic rate constant is an average over such "microscopic" rate constants, weighted by the probability of finding a molecule with a given set of quantum numbers.
* For systems in equilibrium, the probability of finding a molecule in a certain state depends on its energy by means of the Boltzmann distribution and the macroscopic rate constant thereby becomes a function of temperature.
* Stable molecules correspond to minima on the potential energy surface within the Born-Oppenheimer approximation and a chemical reaction can be described as nuclei moving from one minimum to another.
* In the lowest level of approximation, the motion is assumed to occur along the path of least energy and this path forms the basis for the **transition state theory (TST)**.
* The *transition state (TS)* is the configuration that divides the reactant and product parts of the potential energy surface.
* In the multidimensional case, the TS is a first-order saddle point on the potential energy surface, a maximum in the reaction coordinate direction and a minimum along all other directions.
* TST is a semi-classical theory where the motion along the reaction coordinate is treated clasically, while the perpendicular directions take into account the quantization of, for example, the vibrational energy.
* The probability of finding a molecule in a given quantum state is proportional to $\exp(-\Delta E/k_B T)$, which is a Boltzmann distribution.
* Assuming that the molecules at the TS are in equilibrium with the reactant, the macroscopic rate constant can be expressed as
$$
$$
* $\Delta G^{\ddagger}$ is the Gibbs free energy difference between the TS and reactant, and $k_B$ is the Boltzmann's constant.
* The TST expression only holds if all molecules that pass from the reactant over the TS go on to the product.
* The dividing surface separating the reactant from the product is a hyperplane perpendicular to the reaction coordinate at the TS.
* The TST assumption is that no-crossings occur, that is all molecules passing through the dividing surface will go on to form the product. These assumptions mean that the rate constant calculated from Eq.(14.2) will be an upper limit to the true rate constant.

### Statistical mechanics
* Most experiments are performed on macroscopic samples, which usually contains ~$10^20$ particles. Calculations, on the other hand, are performed on relatively few particles, typically $1-10^3$ particles.
* The macroscopic result of an experimental measurement can be connected with properties of the microscopic system. The temperature, for example, is related to the average kinetic energy of the particles, as
$$
\braket{E_{kin}} = \frac{3}{2}RT
$$
* The connection between properties of a microscopic system and a macroscopic sample is provided by statistical mechanics.
* At a temperature of 0 K, all molecules are in their energetic ground state but at a finite temperature there is a distribution of molecules in all possible quantum energy states.
* The relative probability $P$ of a molecule being in a state with an energy $\mathcal{E}$ at a temperature $T$ is given by a Boltzmann factor, as
$$
P \propto \exp\left(-\frac{\mathcal{E}}{k_B T}\right)
$$
* The exponential dependence on the energy means that there is a low (but non-zero) probability for finding a molecule in a high-energy state.
* The key feature in statistical mechanics is the **partition function**. Just as the wave function is the cornerstone in quantum mechanics, the partition function allows calculation of all macroscopic functions in statistical mechanics.
* The partition function for a single molecule is usually denoted $q$, and is defined as a sum of exponential terms involving all possible quantum energy states.
$$
q = \sum_{i={\rm states}}^\infty \exp\left(-\frac{\mathcal{E}_i}{k_B T}\right)
$$
* The partition function can also be written as a sum over all distinct energy levels, multiplied with a degeneracy factor $g_i$ that indicates how many states are with the same energy $\mathcal{E}_i$.
$$
q = \sum_{i={\rm levels}}^\infty g_i \exp\left(-\frac{\mathcal{E}_i}{k_B T}\right)
$$
* The partition function may also viewed as the normalization factor for the Boltzmann probability distribution, as
$$
P(\mathcal{E}_i) = \frac{\exp\left(-\frac{\mathcal{E}_i}{k_B T}\right)}{q}
$$
* $q$ is the partition function for a single particle, and corresponding quantity $Q$ for a collection on $N$ non-interacting particles (ideal gas) is given as
$$
Q = q^N \ \ (\text{different particles}) \\
Q = \frac{q^N}{N!} \ \ (\text{identical particles})
$$
* If the particles are interacting (liquid or solid states), the partition function $Q$ must be calculated by summing over all energy states $\mathcal{E}_i$ for the whole system.
* The significance of the partition function $Q$ is that thermodynamic functions, such as the internal energy $U$ and Helmholtz free energy $A \ (A=U-TS)$ can be calculated from it.
$$
U = k_B T^2 \left(\frac{\partial \ln Q}{\partial T}\right)_V \\
A = -k_B T \ln Q
$$
* Macroscopic observables, such as pressure $P$ and heat capacity (at constant volume) $C_V$, may be calculated as derivatives of thermodynamic functions.
$$
P = \left(\frac{\partial A}{\partial V}\right)_T = k_B T \left(\frac{\partial \ln Q}{\partial V}\right)_T \\
C_V = \left(\frac{\partial U}{\partial T}\right)_V = 
2 K_B T \left(\frac{\partial \ln Q}{\partial T}\right)_V + k_B T ^2 \left(\frac{\partial^2 \ln Q}{\partial T^2}\right)_V
$$
* Other thermodynamic functions, such as the enthalpy $H$, the entropy $S$, and Gibbs free energy $G$ may be constructed from these relations
$$
H = U + PV = k_B T^2 \left(\frac{\partial \ln Q}{\partial T}\right)_V + k_B T V \left(\frac{\partial \ln Q}{\partial V}\right)_T \\
S = \frac{U-A}{T} = k_B T \left(\frac{\partial \ln Q}{\partial T}\right)_V + k_B \ln Q \\
G = H - TS = k_B T V \left(\frac{\partial \ln Q}{\partial V}\right)_T - k_B T \ln Q
$$
* In order to calculate the partition function, one needs to know all possible quantum states for the system. In principle, these can be calculated by solving the nuclear Schrodinger equation.
* For an isolated polyatomic molecule, the energy levels for a single conformation can be calculated within the *rigid-rotor harmonic oscillator approximation*, where the electronic, vibrational, rotational, and translational degrees of freedom are assumed to be separable.

### The ideal gas, rigid-rotor harmonic-oscillator approximation
* For an isolated molecule, the total energy can be approximated as a sum of terms involving translational, rotational, vibrational, and electronic states, and this is a good approximation for the large majority of the systems.
* The assumption that the energy can be written as a sum of terms implies that the partition function can be written as a product of terms.
* As the enthalpy and entropy contributions involve taking the logarithm of $q$, the product of $q$'s thus transforms into sums of enthalpy and entropy contributions.
$$
\mathcal{E}_{tot} = \mathcal{E}_{trans} + \mathcal{E}_{rot} + \mathcal{E}_{vib} + \mathcal{E}_{elec} \\
H_{tot} = H_{trans} + H_{rot} + H_{vib} + H_{elec} \\
S_{tot} = S_{trans} + S_{rot} + S_{vib} + S_{elec} \\
q_{tot} = q_{trans}q_{rot}q_{vib}q_{elec}
$$
* For each of the partition functions, the sum over allowed quantum states runs to infinity.
* However, since the energies becomes larger, the partition functions are finite. Let us examine each of the $q$ factors in a little more detail.

#### Translational degrees of freedom
* The translational degrees of freedom can be exactly separated from the other $3N-3$ coordinates.
* The allowed quantum states for the translational energy are determined by placing the molecule in a "box", i.e. the potential is zero inside the box but infinite outside.
* The solutions to the Schrodinger equation for such a "particle in a box" are standing waves, cosine and sine functions.
* The energy levels are associated with a quantum number $n$, and depend only on the total molecular mass $M$,
$$
\mathcal{E}_n = \frac{n^2h^2}{8\pi^2 M}
$$
* Although the energy levels are quantized, the energy difference between levels is so small that the distribution can be treated as continuous.
* The summation involved in the partition function can therefore be replaced by an integral (an integral is just a sum in the limit of infinitely small contributions).
$$
q_{trans} = \sum_{n=0}^\infty\exp\left(-\frac{\mathcal{E}_n}{k_B T}\right)
    \approx \int_{n=0}^\infty dn \exp\left(-\frac{\mathcal{E}_n}{k_B T}\right)
$$
* Inserting the energy expression and performing the integration gives
$$
q_{trans} = \left(\frac{2\pi M k_B T}{h^2}\right)V
$$
* The only molecular parameter that enters is the total molecular mass $M$.
* The volume $V$ depends on the number of particles, and it is customary to work on a molar scale, in which $V$ is the volume of 1 mol of (ideal) gas.

#### Rotational degrees of freedom
* Within the rigid-rotor approximation, the rotation of the molecule is assumed to occur with a fixed geometry.
* The energy levels calculated from the rotational Schrodinger equation for a diatomic "rigid rotor" are given in terms of a quantum number $J$ and the moment of inertia $I$.
$$
\mathcal{E}_J = J(J+1)\frac{h^2}{8\pi^2 I}
$$
* $J$ runs from zero to infinity.
* The moment of inertia is calculated from the atomic masses $m_1$ and $m_2$ and the distances $r_1$ and $r_2$ of the nuclei relative to the center of mass.
$$
I = m_1 r_1^2 + m_2 r_2^2
$$
* For all molecules except very light species such as H2 or LiH, the moment of inertia is so large that the spacing between the rotational energy levels is much smaller than $k_B T$ at ambient temperatures.
* As for $q_{trans}$, this means that the summation in Eq.(13.10) can be replaced by an integral.
$$
q_{rot} = \sum_{J=0}^\infty \exp\left(-\frac{\mathcal{E}_J}{k_B T}\right) \approx \int_0^\infty dJ \exp\left(-\frac{\mathcal{E}_J}{k_B T}\right)
$$
* This yields
$$
q_{rot} = \frac{8\pi^2 I k_B T}{h^2 \sigma^2}
$$
* The symmetry number $\sigma$ is 2 for homonuclear system and 1 for a heteronuclear diatomic molecule.
* For a polyatomic molecule, the equivalent of Eq.(13.24) is a 3x3 matrix,
$$
{\bf I} = 
\begin{pmatrix}
\sum_i m_i (y_i^2 + z_i^2) & -\sum_i m_i x_i y_i & -\sum_i m_i x_i z_i \\
-\sum_i m_i x_i y_i & \sum_i m_i (x_i^2 + z_i^2) & -\sum_i m_i y_i z_i \\
-\sum_i m_i x_i z_i & -\sum_i m_i y_i z_i & \sum_i m_i (x_i^2 + y_i^2)
\end{pmatrix}
$$
* Here the coordinates are again relative to the center of mass.
* By choosing a suitable coordinate transformation, this matrix may be diagonalized, with the eigenvalues being the *moments of inertia* and the eigenvectors called *principle axes of inertia*.
* For a general polyatomic molecule, the rotational energy levels cannot be written in a simple form. A good approximation, however, can be obtained from classical mechanics, resulting in the following partition function
$$
q_{rot} = \frac{\sqrt{\pi}}{\sigma}\left(\frac{8\pi^2k_B T}{h^2}\right)^{3/2}\sqrt{I_1 I_2 I_3}
$$
* Here $I_i$ are the three moments of inertia.
* The symmetry number $\sigma$ is the order of the rotational subgroup in the molecular point group.
* The rotational partition function requires only information about the atomic masses and positions, i.e. the molecular geometry.

#### Vibrational degrees of freedom
* In the lowest approximation level, the molecular vibrations may be described by a harmonic oscillator. This can be derived by expanding the energy as a function of the nuclear coordinates $R$ in a Talyer series around the equilibrium geometry $R_0$.
$$
E(R) = E(R_0) + \frac{dE}{dR}(R-R_0) + \frac{1}{2}\frac{d^2E}{dR^2}(R-R_0)^2 + \cdots
$$
* The first term may be taken as zero, since this is just the zero point for the energy. The second term vanishes since the expansion is around the equilibrium geometry.
* Keeping only the lowest non-zero term results in the harmonic approximation, where $k$ is the force constant.
$$
E(\Delta R) \approx \frac{1}{2}\frac{d^2E}{dR^2}\Delta R^2 = \frac{1}{2}k\Delta R^2
$$
* Thus the energy levels obtained from the Schrodinger equation for a one-dimensional harmonic oscillator (diatomic system) are given as
$$
\begin{align*}
\epsilon_n &= \left( n + \frac{1}{2} \right) h\nu \\
\nu &= \frac{1}{2\pi}\sqrt{\frac{k}{\mu}} \\
\mu &= \frac{m_1 m_2}{m_1 + m_2}
\end{align*}
$$
* Here, $n$ is the vibrational quantum number running from zero to infinity and $nu$ is the vibrational frequency given in terms of the force constant $k$ ($=\frac{d^2 E}{dR^2}$) and the reduced mass $\mu$.
* In contrast to the translational and rotational energy levels, the spacing between vibrational energy levels is comparable to $k_B T$ for $T$ around 300 K. Thus the summation for the vibrational partition function $q_{vib}$ cannot be replaced by an integral. However, due to the regular spacing, the infinite summation can be written in a closed form as
$$
\begin{align*}
q_{vib} 
&= \sum_{n=0}^\infty \exp(-\epsilon_n/k_B T) = \exp(-h\nu/2k_B T) + \exp(-3h\nu/2k_B T) + \exp(-5h\nu/2k_B T) + \cdots \\
&= \exp(-h\nu/2k_B T)\left\{ 1 + \exp(-h\nu/k_B T) + \exp(-2h\nu/k_B T) + \cdots \right\} \\
&= \frac{\exp(-h\nu/2k_B T)}{1 - \exp(-h\nu/k_B T)}
\end{align*}
$$
* In the infinite sum, each successive term is smaller than the previous one by a constant factor $\exp(-h\nu/k_B T)$, which is smaller than one. Therefore, the infinite sum can be expressed in a closed form.
* For a polyatomic molecule, the force constant $k$ is replaced by $3N_{atom} \times 3N_{atom}$ matrix (the Hessian matrix) containing the second derivative of energy with respect to the coordinate as its matrix element.
* By mass-weighting and diagonalizing this matrix, a new coordinate system called as the *vibrational normal modes* is obtained.
* In the vibrational normal modes, the 3N-dimensional Schrodinger equation can be separated into 3N one-dimensional equations, each having the form of a harmonic oscillator.
* Of these 3N modes, three modes describe the translational motion, and three (non-linear molecule) or two (linear molecule) modes describe the rotational motion. Thus, $3N-6$ or $3N-5$ modes describe the vibration.
* In summary, within the harmonic approximation, the vibrational degrees of freedom are decoupled in the normal modes. Since the energy of the $3N-6$ vibrations can be written as a sum the partition function can be written as a product over $3N-6$ vibration partition functions, as
$$
E_{vib} = \sum_{i=0}^{3N-6}\left( n_i + \frac{1}{2} \right)h\nu_i \\
q_{vib} = \sum_{i=0}^{3N-6}\frac{\exp(-h\nu_i/2k_B T)}{1-\exp(-h\nu_i/k_B T)}
$$
* The vibrational frequencies $\nu_i$ are needed for calculation $q_{vib}$. These are obtained from the force constant matrix and atomic masses.
* If the stationary point is a minimum on the potential energy surface, the eigenvalues of the Hessian are all positive. If, however, the stationary point is a TS (a first-order saddle point), one of the eigenvalues becomes negative thus the frequency is imaginary number. In this case, the number of vibrational modes reduced to 3N-7, and the summation up to this number should be taken.

#### Electronic degrees of freedom
* The electronic partition function $q_{elec}$ involves the sum over electronic quantum states.
* These are the solutions to the electronic Schrodinger equation, i.e. the ground state and an possible excited states.
* In almost all molecules, the energy difference between the ground and excited states is large compared with $k_B T$, which means that only the first term (the ground state energy) in the partition function summation is important.
$$
g_{elec} = \sum_{i=0}^\infty g_i \exp(-\epsilon_i/k_B T) \approx g_0 \exp(-\epsilon_0/k_B T)
$$
* $g_i$ is the degeneracy of the $i$-th electronic wave function, and it may be either in the spin part (1 for singlet, 2 for doublet, etc.) or in the spatial part (1 for $A$, $B$, $\Sigma$ representation in the point group, or 2 for $E$, etc.). The large majority of stable molecule have $g_0=1$.

#### Enthalpy and entropy contributions
* Given the partition function, the enthalpy and entropy terms may be calculated by carrying out the required differentiations in Eq.(13.18).
* For one mole of molecules, the results for a non-linear system are
$$
H_{trans} = \frac{3}{2}RT \\
H_{rot} = \frac{3}{2}RT \\
H_{vib} = R\sum_i^{3N-6}\left\{ \frac{h \nu_i}{2k} + \frac{h \nu_i}{k} \frac{1}{\exp(h \nu_i/k_B T)-1} \right\} \\
H_{elec}^{TS} = \Delta E^\ddagger \\
S_{trans} = \frac{5}{2}R + R\ln\left\{ \frac{V}{N_A}\left( \frac{2\pi M k_B T}{h^2} \right)^{3/2} \right\} \\
S_{rot} = R \left[ \frac{3}{2} + \ln\left\{ \frac{\sqrt{\pi}}{\sigma}\left(\frac{8\pi^2 k_B T}{h^2}^{3/2}\sqrt{I_1I_2I_3} \right) \right\} \right] \\
S_{vib} = R\sum_i^{3N-6} \left[ \frac{h\nu_i}{k_B T}\frac{1}{\exp(h \nu_i/k_B T)-1} -\ln\left\{ 1-\exp(-h\nu_i/k_B T) \right\} \right]
$$
* The rotational terms are slightly different for a linear molecule, and the vibrational terms will contain one vibrational contribution more.
$$
H_{rot}({\rm linear}) = RT \\
S_{rot}({\rm linear}) = R\left[1+\ln\left(\frac{8 \pi^2 I k_B T}{\sigma h^2}\right)\right]
$$
* In summary, to calculate rate and equilibrium constants we need to calculate $\Delta G^\ddagger$ and $\Delta G_0$. This can be done within the RRHO approximation if the geometry, energy, and force constants are known for the reactant, TS, and product.
* The translational and rotational contributions are trivial to calculate, while the vibrational frequencies require the full force constant matrix (i.e. all energy second derivatives), which may be a significant computational effort.
* The activation enthalpies and entropies in principle depend on temperature (Eq.X), but only weakly. Thus, for a limited temperature range they may be treated as constants. Obtaining these quantities from experiments is possible by measuring the reaction rate as a function of temperature, and plotting $\ln(k/T)$ against $T^{-1}$.
$$
k = \frac{k_B T}{h}\exp \left( -\frac{\Delta G^\ddagger}{RT} \right) \\
\ln\left(\frac{k_{rate}}{T}\right) = \ln\left(\frac{k_B}{h}\right) + \frac{\Delta S^\ddagger}{R} - \frac{\Delta H^\ddagger}{RT}
$$
* Such plots should produce a straight line with the slope being equal to $-\Delta H^\ddagger / R$ and the intercept equal to $\ln(k_B/h) + \Delta S^\ddagger / R$.
* Experimentalists often analyze their data in terms of an Arrhenius expression instead of the TST expression of Eq.(13.39), by plotting $\ln(k_{rate})$ against $T^{-1}$.
$$
k = A\exp\left(-\frac{\Delta E^\ddagger}{RT}\right) \\
\ln k = \ln A - \frac{\Delta E^\ddagger}{RT}
$$

<!-- Jensen end -->

### 2-1. Thermodynamics
#### Potential energy
* The strength of the interaction is measured by the change in potential energy of the system as a function of the distance $z$, of the adsorbate above the surface:
$$
\Delta E(z) = E_{pot}(z) - E_{pot}(\infty)
$$
* This is called the **potential energy curve**, and its two-dimensional version is called the **potential energy surface (PES)**.
* The lowest energy pathway from one potential minimum on as PES to another is called the **minimum energy path (MEP)**. The MEP is formally defined as the path of least action: at any point along the path, the gradient of the potential has no component perpendicular to the path.
* MEP is often considered as the reaction coordinate, and the highest energy along the MEP is the transition state. This will be discussed later.

#### Zero-point energy
(need some simple introduction)
* The quantized energy solutions of the harmonic oscillator are given by 
$$
E(n) = E_{pot} + h\nu_i \left( n+\frac{1}{2} \right)
$$
where is $\nu_i$ is the vibrational frequency, $h$ is the Planck's constant, and $n$ is the quantum number. The ground state, which is the lowest energy level, $n=0$, thus has an energy
$$
E_0 = E_{pot} + \frac{1}{2}h\nu_i
$$
* The contributions from different vibrational modes are additive. Thus if there are $M$ vibrational modes in the potential energy minimum, the total **zero-point energy (ZPE)** correction would be
$$
ZPE = \sum_i^M \frac{1}{2} g_i h\nu_i
$$
Note that one should add degeneracy factor ($g_i$) for degenerated modes; $g_i = 2 (E)$, and $g_i = 3 (T)$ (CHECK).
* Frequencies larger than 1000 ${\rm cm^{-1}}$ have contributions to the ZPE at ~0.1 eV order. The H2 frequency is, for example, as high as 4395 ${\rm cm^{-1}}$ thus has significant contribution from ZPE.
* The reaction energy for ammonia synthesis, for example, comes out as much as 0.83 eV too exothermic if one does not account for the ZPE corrections.
* Often experimentally measured frequencies are used for the ZPE calculation. NIST webbook is often used.[NIST webbook]

#### Internal energy
* The amount of energy that needs to be transferred to warm up the system (per temperature increase at a given temperature) is called the **heat capacity**.
* The sum of the ZPE-corrected ground-state energy $E_0$ and the thermal energy of the system is called the **internal energy**, and if we know the heat capacity (at constant pressure) this is expressed as 
$$
U(T) = E_0 + \int_{T=0}^{T} C_P(T')dT'
$$

#### Gibbs free energy
* The key thermodynamic concept for describing equilibria in chemical process is the *Gibbs free energy**, which is defined as
$$
G = H - TS
$$
where $H$ is the enthalpy and $S$ is the entropy. Enthalpy is defined from the internal energy $U$ as
$$
H = U - pV
$$
where $p$ is the pressure and $V$ is the volume.
* The potential energy describes a mechanical systems's potential for carrying out mechanical work, while the Gibbs energy describes a closed chemical systems' potential for carrying out non-expansion work.
* The maximum amount of work carried out by a mechanical system is only possible if all motion is frictionless; likewise, the maximum amount of work by a chemical system is achieved if all reactions are reversible.
* The equilibrium of the system requires that any possible chemical reaction that the system could undertake must satisfy the relation
$$
\Delta G = G_{final} - G_{initial} = 0
$$
where $G_{initial}$ and $G_{final}$ are the Gibbs free energies of a given set of atoms before and after, respectively, they undergo some chemical reaction.
* The absolute Gibbs free energy incorporates a certain level of arbitrariness through the dependence of the internal energy on the ground-state potential energy. The commonly accepted standardized choice is to define a standard state for every substance at a set of standard conditions.
* The standard condition for gas is chosen to be p = 1 bar, while for solution the molarity of C = 1 mol solute/kg solvent at the standard pressure of p = 1 bar.
* For historical reasons, $G^0$ is defined to be zero for the pure elements in their reference states. For other substances, $G^0$ is the reaction Gibbs energy for making it in their standard state from the constituting elements in their reference state.
* One typically always need the Gibbs energies of reaction, not the absolute Gibbs energies themselves. The following relation is therefore often used
$$
\Delta G^0 = \Delta H^0 - T\Delta S^0
$$
* The $\Delta$ here indicate a change in the given quantity as a reaction is undertaken.
* For a pure substance X, we can write the Gibbs free energy as 
$$
G(p) = G^0(p) + k_B T \ln \left( \frac{p}{p^0} \right)
$$
where $p^0$ is the standard pressure.
* For a reaction $\ce{A -> B}$, the change in Gibbs free energy is given by
$$
\begin{split}
\Delta G 
 &= \Delta G^0 + k_B T \left( \ln \frac{p_B}{p^0} - \ln \frac{p_A}{p^0} \right) \\
 &= \Delta G^0 + k_B T \ln \frac{p_B}{p^A}
 \end{split}
$$
* At equilibrium $\Delta G = 0$ and we thus obtain
$$
\left. \frac{p_B}{p_A} \right|_{eq} = \exp \left( -\frac{\Delta G^0}{k_B T} \right) = K_{eq}
$$
where we have introduced the equilibrium constant of the reaction, $K_{eq}$, which is unitless.
* Consider as an example of the ammonia synthesis reaction: $\ce{N2 + 3H2 -> 2NH3}$.
* Using the experimental value at 700 K of $\Delta H^0 = -0.95\ eV$ and $\Delta S^0 = -2.05\ meV/K$, $\Delta G^0 = 0.49\ eV$ then $K_{eq} \approx 0.0002$.
* Thus, at this temperature, the equilibrium is strongly shifted toward the reactants (N2 and H2). 

#### Enthalpy
$$
C_P = k_B + C_{V,\rm trans} + C_{V,\rm rot} + C_{V,\rm vib} + C_{V, \rm elec}
$$
* The translational heat capacity is $3/2 k_B$. The rotational heat capacity is $k_B$ and $3/2 k_B T$ for linear and non-linear molecule, respectively. The electronic heat capacity is zero.
* The integrated form of the vibrational heat capacity is
$$
\int C_{V,\rm vib} dT = \sum_{i}^{N_{\rm mode}}\frac{\epsilon_i}{e^{\epsilon_i/k_B T}-1}
$$
* Here, the degree-of-freedom is $3N-5$ and $3N-6$ for linear and non-linear molecules, respectively.

#### Entropy
* The entropy can be thought of as a measure of the number of accessible quantum states of a system. The entropy is proportional to the logarithm of the number of accessible states $\Omega$, as
$$
S = k_B \ln \Omega
$$

* The ideal gas entropy can be calculated as a function of temperature and pressure, as
$$
\begin{split}
S(T,P) &= S(T,P^0) - k_B \ln\frac{P}{P^0} \\
&= S_{\rm trans} + S_{\rm rot} + S_{\rm vib} + S_{\rm elec} - k_B \ln\frac{P}{P^0}
\end{split}
$$
where $S_{\rm trans}, S_{\rm rot}, S_{\rm vib}, S_{\rm elec}$ are translational, rotational, vibrational, and electronic components of entropy. These are calculated as below.
* The translation entropy is
$$
S_{\rm trans} = k_B\left\{ \ln\left[ \left(\frac{2\pi M k_B T}{h^2}\right)^{3/2}\frac{k_B T}{P^0} \right] + \frac{5}{2} \right\}
$$
* The rotational entropy is
$$
S_{\rm rot} = \left\{
\begin{split}
&k_B \left[ \ln\left(\frac{8\pi^2 I k_B T}{\sigma h^2}\right)+1 \right] \\
&k_B \left\{ \ln\left[ \frac{\sqrt{\pi I_A I_B I_C}}{\sigma}\left(\frac{8\pi^2 k_B T}{h^2}\right)^{3/2} \right] + \frac{3}{2} \right\}
\end{split}
\right.
$$
where the upper and lower is linear and non-linear molecules. Here, $I_x$ is the principle moments of inertia. $I$ is the degenerate moments of inertia, and $\sigma$ is the symmetry number of the molecule.
* The vibrational entropy is
$$
S_{\rm vib} = k_B \sum_i \left[ \frac{\epsilon_i}{k_B T(e^{\epsilon_i/k_B T}-1)} - \ln(1-e^{-\epsilon_i/k_B T}) \right]
$$
* The electronic entropy is
$$
S_{\rm elec} = k_B (2S+1)
$$
where $S$ is the spin multiplicity.

### 2-2. Chemical Kinetics
#### Stoichiometry and reaction rate
* We assume that a chemical process can be described by a single **stoichiometric equation** or **overall reaction equation**. Real chemical systems corresponding to such a single chemical reaction, that is, when the reactants react with each other forming products immediately, are in fact very rare.
* In most cases, the reaction of the reactants produces intermediates, these intermediates react with each other, and the final products are formed at the end of many coupled reaction steps. Each of the individual steps is called an **elementary reaction** or **elementary step**.
* Within elementary reactions, these is no macroscopically observable intermediates between the reactants and the products.

* The overall reaction equation of the production of water from hydrogen and oxygen is very simple:
$$
\ce{2H2 + O2 = 2H2O}
$$
* From a stoichiometric point of view, a chemical equation can be rearranged, similarly to mathematical equation. For example, all terms can be shifted to the right-hand side as
$$
\ce{0 = -2H2 -O2 + 2H2O}
$$
* Let us denote the formula of the chemical species by the vector ${\bf A}=(A_1, A_2, A_2)=("\ce{H2}", "\ce{O2}", "\ce{H2O}")$ and $\nu_1 = -2, \nu_2 = -1, \nu_3 = +2$. The corresponding general stoichiometric equation is
$$
0 = \sum_{j=1}^{N_s}\nu_j A_j
$$
where $N_s$ is the number of species.
* The general stoichiometric equation of any chemical process can be defined in the same way, where $\nu_j$ is the stoichiometric coefficient of th $j$-th species and $A_j$ is the chemical formula of the $j$-th species in the overall reaction equation.
* We show here the stoichiometric coefficients for an overall reaction equation, but the same approach is taken for the elementary reactions. In general, for elementary reactions within a chemical mechanism, the stoichiometric coefficients are integers.

* Let us think about the time-dependent behavior of a chemical system and how we might describe it using information from the kinetic reaction system.
* By measuring the molar concentration $Y_j$ of $j$-th species at several consecutive time points, then by applying a finite-difference approach, the production rate of the $j$-th species $dY_j/dt$ can be calculated.
* The rate of chemical reaction defined by stoichiometric equation is the following:
$$
r = \frac{1}{\nu_j}\frac{dY_j}{dt}
$$
* Reaction rate $r$ is independent of index $j$. This means that the same reaction rate is obtained when the production rate of any of the species is measured.
* Within a narrower range of concentrations, the reaction rate $r$ can always be approximated by the following equation
$$
r = k\prod_{j=1}^{N_s}Y_j^{\alpha_j}
$$
where the positive scalar $k$ is called the **rate coefficient**. The exponents $\alpha_j$ are positive real numbers or zero.
* When the reaction rate is calculated by Eq.X, molar concentrations (such as $mol cm^{-3}$) should always be used.
* The rate coefficient $k$ is independent on the concentration, but may depend on temperatures, pressures, etc. Thus, it is better to call this quantity as *rate coefficient* and not *rate constant*.
* The exponent $\alpha_j$ is called the *reaction order of species $A_j$* and the sum of exponents ($\alpha = \sum_j A_j$) is called the *overall order of the reaction*.
* The physical dimension of the rate coefficient $k$ depends on the overall order of the reaction step. When the order of the reaction step is 0, 1, 2 or 3, the dimension of $k$ is $\text{(concentration)}\times\text{(time)}^{-1}$, $\text{(time)}^{-1}$, $\text{(concentration)}^{-1}\times\text{(time)}^{-1}$, or $\text{(concentration)}^{-2}\times\text{(time)}^{-1}$, respectively.

* As stated above, intermediates are formed within most reaction systems. Hence, in order to define the time-dependent dynamics of a system appropriately, a reaction model should include elementary reactions were such intermediates are formed from reactants and go on to form products.
* For example, we can propose a simple reaction mechanism for the hydrogen combustion reaction (Eq.X) as
$$
\begin{align*}
&R1\hspace{5mm} \ce{H2 + O2      = H + HO2} \\
&R2\hspace{5mm} \ce{O2 + H       = OH + O} \\
&R3\hspace{5mm} \ce{H2 + OH      = H + H2O} \\
&R4\hspace{5mm} \ce{H2 + O       = H + OH} \\
&R5\hspace{5mm} \ce{O2 + H (+ M) = HO2 (+ M)} \\
&R6\hspace{5mm} \ce{H2O + OH     = H2O + OH} \\
\end{align*}
$$
where $M$ represents any species in the mixture.
* Each elementary reaction step $i$ can be characterized by the following stoichiometric equation:
$$
\sum_j\nu_{ij}^L A_j = \sum_j\nu_{ij}^R A_j
$$
where the stoichiometric coefficient on the left-hand side ($\nu_{ij}^L$) and the right-hand side ($\nu_{ij}^R$) of an elementary reaction $i$ should be distinguished. The net stoichiometric coefficient of species $j$ can be obtained as $\nu_{ij} = \nu_{ij}^R - \nu_{ij}^L$.
* $\nu_{ij}^L$ and $\nu_{ij}^R$ should be positive integers, so $\nu_{ij}$ can be positive or negative.
* Elements $\nu_{ij}^L$, $\nu_{ij}^R$, and $\nu_{ij}$ constitute the left-hand side, the right-hand side, and the overall *stoichiometric matrix*, respectively.
* To emphasize the analogy with mathematical equations, so far the equality sign ($=$) was always used for chemical reactions. From now on, arrows will be used for one-way (or irreversible) chemical reactions (such as $\ce{A -> B}$). Reversible reactions will be denoted by double arrows ($\ce{A <=> B}$).

* The stoichiometric equation of the elementary reaction should corresponds to the real molecular changes.
* Because of this, most of the elementary reactions are *unimolecular* or *bimolecular* reactions, where one and two particles (atoms, molecules, ions) are involved, respectively.
* Unimolecular reactions include the photochemical reaction (e.g. $\ce{NO2 + h\nu -> NO + O}$), the isomerization reaction, etc.
* Most elementary reactions are bimolecular reaction, in which two particles meet and both particles change chemically.
* In some elementary reactions, a *third-body* is involved. This can be a molecule of the bath gas (in most experiments argon or nitrogen).

* The rates of elementary reactions can be calculated by assuming the **law of mass action**, as
$$
r_i = k_k \prod_j^{N_S}Y_j^{\nu_{ij}^L}
$$
where $r_i$ and $k_i$ are the rate and the rate coefficient, respectively, of reaction step $i$ and $Y_j$ is the molar concentration of species $j$.

* The *kinetic system of ordinary differential equations (ODEs)* defines the relationship between the production rates of the species and rates of the reaction steps $r_i$:
$$
\frac{dY_j}{dt} = \sum_i^{N_R}\nu_{ij}r_i \hspace{5mm} j=1,2,\cdots,N_S
$$
* The above equation can also be written in a simpler form using the vector of concentrations ${\bf Y}$, the stoichiometric matrix ${\boldsymbol \nu}$ and the vector of the rates of reaction steps ${\bf r}$, as
$$
\frac{d{\bf Y}}{dt} = {\boldsymbol \nu}{\bf r}
$$
* This means the number of equations in the kinetic system ODEs is equal to the number of species in the reaction mechanism.
* The rates of the reaction steps can be very different and may span (even 10-25) orders of magnitude. Such differential equations are called *stiff ODEs*.

##### Examples
* The example for the creation of the kinetic system of ODEs will be based on a hydrogen combustion mechanism.
* Using the law of mass action, the rates $r_1$ to $r_6$ of the reaction steps can be calculated from the species concentrations and rate coefficients
$$
\begin{align*}
&R1\hspace{5mm} \ce{H2 + O2      = H + HO2}   & &r_1 = k_1[\ce{H2}][\ce{O2}]\\
&R2\hspace{5mm} \ce{O2 + H       = OH + O}    & &r_2 = k_2[\ce{O2}][\ce{H}]\\
&R3\hspace{5mm} \ce{H2 + OH      = H + H2O}   & &r_3 = k_3[\ce{H2}][\ce{OH}]\\
&R4\hspace{5mm} \ce{H2 + O       = H + OH}    & &r_4 = k_4[\ce{H2}][\ce{O}]\\
&R5\hspace{5mm} \ce{O2 + H (+ M) = HO2 (+ M)} & &r_5 = k_5[\ce{O2}][\ce{H}][\ce{M}]\\
&R6\hspace{5mm} \ce{H2O + OH     = H2O + OH}  & &r_6 = k_6[\ce{HO2}][\ce{OH}]\\
\end{align*}
$$
Here $[\ce{M}]$ is the sum of the concentration of all species present.
* The calculation of the production rates is based on Eq.X. For example, the hydrogen atom H is produced in R1, R3, and R4, and, it is consumed in R2 and R5, and it is not present in R6.
* The line of the kinetic system of ODEs, corresponding to the production of H is the following:
$$
\frac{d[\ce{H}]}{dt} = k_1[\ce{H2}][\ce{O2}] - k_2[\ce{O2}][\ce{H}] + k_3[\ce{H2}][\ce{OH}] + k_4[\ce{H2}][\ce{O}] - k_5[\ce{O}][\ce{H}][\ce{M}]
$$
* In a similar way, the production of water can be described by the following equations:
$$
\frac{d[\ce{H2O}]}{dt} = k_3[\ce{H2}][\ce{OH}] + k_6[\ce{HO2}][\ce{OH}]
$$

##### Parametrising rate coefficients
###### Temperature dependence of rate coefficients
* The temperature dependence of rate coefficient $k$ is usually described by the **Arrhenius equation**:
$$
k = A\exp(-E/RT)
$$
where $A$ is the pre-exponential factor or $A$-factor, $E$ is the activation energy, $R$ is the gas constant. The dimension of $E/R$ is temperature, thus $E/R$ is called the *activation temperature*.
* If the temperature dependence of the rate coefficient can be described by the Arrhenius equation, named after a Swedish chemist Svante Arrhenius who proposed this equation in 1889. Plotting $\ln(k)$ as a function $1/T$, or *Arrhenius plot*, gives a straight line. The slope of this line is $-E/R$, and the intercept is $\ln(A)$.
* In high-temperature gas-phase kinetic systems, such as combustion and pyrolytic systems, the temperature dependence of the rate coefficient is usually described by the *modified Arrhenius equation (or extended Arrhenius equation)*:
$$
k = AT^n \exp(-E/RT)
$$

#### Reversible reaction steps
* In theory, all thermal elementary reactions are *reversible*, which means that the reaction products may react with each other to reform the reactants.
* Within the terminology used for reaction kinetics, a reaction step is called *irreversible*, either if the backward reaction is not taken into account or the reverse reaction is represented by a pair of opposing irreversible reaction steps.
* The irreversible and reversible reactions are denoted by a single arrow $\ce{->}$ and the two-way arrow $\ce{<=>}$, respectively.
* In such cases, a forward rate expression may be given such as Arrhenius forms, and the reverse rate is calculated from the thermodynamics properties of the species through the equilibrium constants. Hence, if the forward rate coefficient $k_{f_i}$ is known, the reverse rate coefficient can be calculated from
$$
k_{r_i} = \frac{k_{f_i}}{K_{c_i}}
$$
where $K_{c_i}$ is the equilibrium constant expressed in molar concentrations. $K_{c_i}$ is obtained from the thermodynamic property of species.

(overlapping with CFD section)
* In combustion systems, thermodynamics properties are often calculated from 14 fitted polynomial coefficients called NASA polynomicals for each species. Seven are used for the low-temperature range $T_{\rm low}$ to $T_{\rm mid}$ and seven for the high-temperature range $T_{\rm mid}$ to $T_{\rm high}$. Typical values are $T_{\rm low} = 300\ K$, $T_{\rm mid} = 1000\ K$, and $T_{\rm high} = 5000\ K$.
* The polynomial coefficients are determined by fitting to tables of thermochemical or thermodynamics properties, which are either measured values or calculated using theoretical methods and statistical thermodynamics.
* The polynomial coefficients can then be used to evaluate various properties at a given temperature ($T$), such as standard molar heat capacity ($C_p^{\rm o}$), enthalpy ($H^{\rm o}$) and entropy ($S^{\rm o}$) as follows:
$$
\begin{align*}
\frac{C_p^{\rm o}}{R} &= a_1 + a_2 T + a_3 T^2 + a_4 T^3 + a_5 T^4 \\
\frac{H^{\rm o}}{RT}  &= a_1 + \frac{a_2}{2}T + \frac{a_3}{3}T^2 + \frac{a_4}{4}T^3 + \frac{a_5}{5}T^4 + \frac{a_6}{T} \\
\frac{S^{\rm o}}{R}   &= a_1\ln(T) + a_2 T + \frac{a_3}{2}T^2 + \frac{a_4}{3}T^3 + \frac{a_5}{4}T^4 + a_7 \\
\end{align*}
$$
where $a_n$ parameters are the NASA polynomial coefficients.
(end overlapping)
* The standard molar reaction enthalpy of species $j$ ($\Delta_r H_j^{\rm o}$) and entropy ($\Delta_r S_j^{\rm o}$) can be calculated from the following equations:
$$
\frac{\Delta_r S_j^{\rm o}}{R} = \sum_{i=1}^I \nu_{ij}\frac{S_i^{\rm o}}{R} \\
\frac{\Delta_r H_j^{\rm o}}{R} = \sum_{i=1}^I \nu_{ij}\frac{H_i^{\rm o}}{RT}
$$
* The equilibrium constant $K$ in terms of normalized pressure $p/p^{\rm o}$ is then obtained from
$$
\Delta_r G^{\rm o} = -RT\ln K = \Delta_r H^{\rm o} - T \Delta_r S^{\rm o} \\
K = \exp\left(-\frac{\Delta_r H^{\rm o}}{RT}+\frac{\Delta_r S^{\rm o}}{R}\right)
$$
where $p^{\rm o}$ is the standard pressure (such as $10^5$ Pa, 1 atm, or 1 bar).
* The equilibrium constant in concentration units $K_c$ is related to the equilibrium constant in normalized pressure units $K$ by the following:
$$
K_c = K\left(\frac{p^{\rm o}}{RT}\right)^{\Delta \nu}
$$
where $\Delta \nu = \sum_i \nu_i$ is the sum of stoichiometric coefficients.
* In this way, the reverse rate coefficient for a thermal reaction can be defined by its forward rate coefficient and the appropriate NASA polynomials for the component species within the reaction.

#### Sensitivity analysis
##### Local sensitivity analysis
* For a spatially homogeneous, dynamical system, the change of the concentrations in time can be calculated by solving the following initial value problem:
$$
\frac{d {\bf Y}}{dt} = {\bf f}({\bf Y}, {\bf x}), \hspace{5mm} {\bf Y}(t_0) = {\bf Y}_0
$$
where the parameter vector ${\bf x}$ having $m$ elements may include rate coefficients, Arrhenius parameters, thermodynamics data, etc.
* Let us now look at the effect of changing parameter values on the solution over time. Assume that the solution of the system of differential equation above is calculated from time $t = 0$ to time $t = t_1$. Here the value of parameter $j$ is changed by $\Delta x_j$ and the solution of the ODE is continued until time $t_2$. We denote $Y_i(t_2)$ as the original solution and $\tilde{Y}_i(t_2)$ as the modified solution at time $t_2$. The sensitivity coefficient can be approximately calculated by a finite-difference approximation:
$$
\frac{\partial Y_i}{\partial x_j}(t_1, t_2) \approx \frac{\Delta Y_i(t_2)}{\Delta x_j} = \frac{\tilde{Y}_i(t_2)-Y_i(t_2)}{\Delta x_j}
$$
* The effect of changes in parameter set ${\bf x}$ on the concentrations at a given time can be characterized by the following Taylor expansion:
$$
\begin{align*}
Y_i(t, x+\Delta x) =& Y_i(t,x) + \sum_{j=1}^m\frac{\partial Y_i}{\partial x_j}\Delta x_j \\
&+\sum_{k=1}^m\sum_{j=1}^m\frac{\partial^2 Y_i}{\partial x_k \partial x_j}\Delta x_k \Delta x_j + \cdots
\end{align*}
$$
* Here, the partial derivative $\partial Y_i/\partial x_j$ is called the first-order local *sensitivity coefficient*. The second-order partial derivative $\partial^2 Y_i/\partial x_k \partial x_j$ is called the second-order local sensitivity coefficient, etc.
* Usually, only the first-order sensitivity coefficients $\partial Y_j/\partial x_j$ are calculated and interpreted.
* The local sensitivity coefficient shows how the model solution $Y_i$ changes as a consequence of a change in value of parameter $x_j$ by a small amount, assuming that all other parameters are fixed at their nominal values.
* The *sensitivity matrix* $S$ has $S_{ij} = \partial Y_i/\partial x_j$ as its element, thus this matrix is the linear approximation of the effects of parameter changes on the model solutions.
* The value of the sensitivity coefficient depends on the units used, such as the rate coefficients. However, the unit of the rate coefficients differs for the first-order reaction, second-order reaction, etc. Thus, to make the sensitivity coefficients comparable, *normalized sensitivity coefficients* $x_j/Y_i(\partial Y_i/\partial x_j)$ is often used. Normalized sensitivity coefficients are dimensionless and their values are independent of the units of the model solution and the model parameters.
* If the model solution $Y_i$ can be zero, usually *semi-normalized sensitivity coefficients* $x_j(\partial Y_i/\partial x_j)$ are calculated.
* The sensitivity coefficients can be used to assess the effect of changing one of the parameters by $\Delta x_j$ at time $t_1$ as follows:
$$
\tilde{Y}_i(t_2) \approx Y_i(t_2) + \frac{\partial Y_i}{\partial x_j}\Delta x_j
$$
* When several parameters are changed simultaneously at time $t_1$ according to the parameter vector $\Delta {\bf x}(t_1)$, then the "perturbed" solution $\tilde{\bf Y}$ at time $t_2$ can be calculated knowing the original solution ${\bf Y}(t_2)$ and the sensitivity matrix, as
$$
\tilde{\bf Y}(t_2) = {\bf Y}(t_2) + {\bf S}(t_1, t_2) \Delta {\bf x}(t_1)
$$
* The sensitivity matrix depends on the time $t_1$ of the parameter perturbation and the time $t_2$ of inspection of the effect of the parameter perturbation. Usually, time $t_1$ is set to be identical to the initial time of the kinetic system of differential equations, which is usually $t = 0$.

##### The brute force method
* In general, the local sensitivity matrix can only be determined numerically. If the original system of kinetic differential equations can be solved numerically, then the local sensitivity matrix can also be calculated using the finite-difference approximations.
* To calculate the sensitivity matrix in this way, we have to known the original solution and the $m$ solution obtained by perturbing each parameter one by one. All in all, the kinetic system of OEDs has to be solved $m+1$ times. This procedure is called the *brute force method*.

##### Principal component analysis of the sensitivity matrix
* Elements of the local concentration sensitivity matrix show the effect of changing a single parameter on the calculated concentration of a species. However, we are frequently interested in the effect of parameter changes on the concentrations of a group of species. This effect is indicated by the *overall sensitivity*,
$$
B_j = \sum_i\left(\frac{x_j}{Y_i}\frac{\partial Y_i}{\partial x_j}(t)\right)^2
$$
* Quantity $B_j$ shows the effect of changing parameter $x_j$ on the concentration of all species present in the simulation at time $t$.
* The utilization of such overall sensitivity measures for identifying unimportant species or reactions, as part of the model reduction process which will be discussed further in Sec. X.
* If a group of species is important for us in the time interval $[t_1, t_2]$, we may be interested in which parameters are highly influential on the measured concentrations. To answer this question, the following scalar valued function, called the objective function, will be used.
$$
Q(x) = \int_{t_1}^{t_2}\sum_i\left(\frac{\tilde{Y}_i(t)-Y_i(t)}{Y_i(t)}\right)^2 dt
$$
where $Y_i(t)$ is the value of variable $i$ at time $t$ calculated by the original parameter set and $\tilde{Y}_i(t)$ is the corresponding value calculated with the altered parameter set. The objective function shows the relative deviation of the two values, integrated over the time interval $[t_1, t_2]$.
* The *principal component analysis* of matrix ${\bf S}$ investigates the effect of the change in parameters on the value of the objective function. The objective function $Q$ can be approximated using the local sensitivity matrix ${\bf S}$:
$$
Q({\boldsymbol \alpha}) = (\Delta{\boldsymbol \alpha}^T)\tilde{{\bf S}}^{\rm T}\tilde{{\bf S}}(\Delta{\boldsymbol \alpha})
$$
where $\Delta {\boldsymbol \alpha} = \Delta \ln{\bf x}$, index $T$ indicates the transpose and matrix $\tilde{\bf S}$ is defined in the following way:
$$
\tilde{\bf S} = [\tilde{\bf S}_1, \tilde{\bf S}_2, \cdots, \tilde{\bf S}_l]^T
$$
* Sensitivity matrices $\tilde{\bf S}_m = {\partial \ln Y_i(t_m)/\partial \ln x_k}$ belong to a series of $l$ time points within time interval $[t_1, t_2]$, and the rows of these matrices belong to the variables present in the summation of the objective function. From the calculation of matrix $\tilde{{\bf S}}^{\rm T}\tilde{{\bf S}}$, the sum of elements of matrices $\tilde{\bf S}_m$ is obtained, and therefore the integral present in Eq.X is replaced by a summation in Eq.X.
* Let us denote ${\boldsymbol \lambda}$ as the vector of eigenvalues of $\tilde{{\bf S}}^{\rm T}\tilde{{\bf S}}$ and ${\bf U}$ the matrix of eigenvectors. Another form of the objective function is the following:
$$
Q({\boldsymbol \alpha}) = \sum_i\lambda_i(\Delta \Psi_i)^2
$$
where the transformed parameters $\Delta {\boldsymbol \Psi} = {\bf U}^{\rm T}\Delta{\boldsymbol \alpha}$ are called principal components.
* The eigenvector belonging to the largest eigenvalue defined a parameter groups which have the largest influence on the objective function.
* On the other hand, the eigenvectors belonging to the small eigenvalues has little effect on the solution of the model.
* Thus, the principal component analysis provides a useful way to interpret complex sensitivity information and to identify important parameters that influence target outputs.

#### Reduction of reaction mechanisms
##### Directed relation graph method
* Methods for species and reaction removal based on directed relation graphs (DRGs) with specified accuracy requirements have been introduced by Lu and Law[Ref].
* In their development of the method, Lu and Law suggested that graph-based methods are highly suited to exploring couplings between species.
* This means that such methods may be applied to remove groups of species that may be internally coupled, through, for example, fast reactions, but are not strongly coupled to important processes with the mechanism.
* Each node in the DRG represents a species from the mechanism, and an edge from vertex A to vertex B exists if and only if the removal of species B would directly induce significant error to the production rate of species A.
* This means that an edge from A to B means that B has to be kep in the mechanism to correctly evaluate the production rate of species A.
* DRG methods start from the selection of important species, called "target species" in the DRG terminology. Using a DRG method, all species closely connected to the target species are identified.
* The various DRG-based reduction methods all state a *connection weight* between paris of species.
* These weights define the directed relation graph structure.
* Starting from the target species, an importance coefficient is calculated for all other species, which quantifies how strongly a given species is connected to the target species.
* Then, all species are eliminated from the mechanism (with their reactions) whose importance coefficient is below a user-defined threshold.
* The original DRG method of Lu and Law defines the connection weight from species $i$ to species $j$ in the following way:
$$
R_{i \rightarrow j} = \frac{\sum_{\alpha \in C(i,j)}|\nu_{i \alpha}r_{\alpha}|}{\sum_{\alpha \in R(i)}|\nu_{i \alpha}r_{\alpha}|}
$$
where $R(i)$ is the set of reactions that are related to species $i$. $C(i,j)$ is the set of reactions in which both species $i$ and $j$ participate. $\nu_{i \alpha}$ is the stoichiometric coefficient of species $i$ in reaction $\alpha$ and $r_{\alpha}$ is the net reaction rate i.e. the difference of the forward and backward rates.
* A variant of the DRG method was suggested by Luo et al for the reduction of reaction mechanism containing many isomers. Luo et al. recommended the application of the maximum norm instead of the summation:
$$
R_{i \rightarrow j} = \frac{\max_{\alpha \in C(i,j)}|\nu_{i \alpha}r_{\alpha}|}{\max_{\alpha \in R(i)}|\nu_{i \alpha}r_{\alpha}|}
$$
* The original DRG method of Lu and Law defines the importance coefficient of species $i$ as
$$
I_i =
\left\{
\begin{align*}
& 1 \hspace{5mm} (\text{if species i is a target species}) \\
& \max_{j \in S}(\min(R_{j \rightarrow i},\ I_j)) \hspace{5mm} (\text{otherwise})
\end{align*}
\right.
$$
* Here $S$ is the full set of chemical species and $R_{j\rightarrow i}$ is a connection weight, and it is simply assumed that if two species are not connected, then $R_{j\rightarrow i}$ = 0.
* This approach defines the importance coefficient for species $i$ as the smallest connection on any path towards a target species.
* $I_i$ is calculated iteratively.
* A small threshold value $\epsilon$ can be defined, and if $I_i < \epsilon$, then species $i$ is considered to be redundant for the simulation of the target species.
* Hence, in Fig.X, if A is an target species, then D must be retained within the scheme since although it is not directly coupled to A, it is part of the dependent set of A by meing directly coupled to B, where B is coupled to A.

* The DRG method was first applied to model system of ethylene combustion (Lu and Law 2005, Lue 2011) with a full scheme of 70 species.
* A value of $\epsilon$ of 0.16 gave a skelton of 33 species, i.e. quite a substantial degree of reduction.
* In application to $n$-heptane and $iso$-octane combustion using full schemes of 561 and 857 species (Lu and Law 2006), $\epsilon$ values of 0.19 and 0.17 resulted in reduced schemes of 188 and 233 species, respectively.

##### from CHEMKIN
* The directed relation graph (DRG) method identifies unimportant species in a reaction mechanism by resolving species coupling without any a priori knowledge of the system.
* It has a distinct advantage over the principle component analysis in that it does not require the sensitivity analysis results.
* Direct species coupling can be defined by the immediate error to the production rate of a species $A$, introduced by the removal of another species $B$ from the mechanism. Such immediate error, noted as $r_{AB}$, can be expressed as
$$
r_{AB} = \frac{\sum|\nu_{A,i}\omega_i\delta_{B,i}|}{\sum|\nu_{A,i}\omega_i|}
$$
where
$$
\delta_{B,i} =
\left\{
\begin{align*}
& 1 \hspace{5mm} (\text{if the i-th reaction involves species B}) \\
& 0 \hspace{5mm} (\text{otherwise})
\end{align*}
\right.
$$
Here, $A$ and $B$ indicate the species, $i$ indicate the $i$-th reaction of the mechanism, $\nu_{A_i}$ indicates the stoichiometric coefficient of species $A$ in teh $i$-th reaction, and $\omega_i$ is the reaction rate of the $i$-th reaction.
* When the immediate error $r_{AB}$ is larger than a user-specified error-tolerance level, it means that the removal of species $B$ from the mechanism will introduce an error to the production rate of species $A$ that is beyond that acceptable to the user.

##### ?
* One can combine a series of elementary reactions into a complete catalytic process. For example in the ammonia synthesis, the overall reaction
$$
\ce{N2 + 3H2 -> 2NH3}
$$
can be divided into following six elementary reactions
$$
\begin{split}
\ce{& N2 + 2$*$     -> 2N$*$} \\
\ce{& H2 + 2$*$     -> 2H$*$} \\
\ce{& N$*$ + H$*$   -> NH$*$ + $*$} \\
\ce{& NH$*$ + H$*$  -> NH2$*$ + $*$} \\
\ce{& NH2$*$ + H$*$ -> NH3$*$ + $*$} \\
\ce{& NH3$*$        -> NH3 + $*$} \\
\end{split}
$$
* Each elementary reactions have their own PED. Thus, one can combine these PEDs into one, which consists of the PED for the full ammonia synthesis reaction. The potential energy difference between the initial and final states is the reaction energy of NH3 synthesis. 

#### TOF
* For heterogeneous reactions, the areal rate is a good choice. However, the catalyst can have the same surface area but different concentrations of active sites. Thus, a definition of the rate based on the number of active sites appears to be the best choice. The **turnover frequency (TOF)** is the number of times the catalytic cycle is completed (or turned-over) per active site per unit time, and expressed as
$$
r_t = \frac{1}{S}\frac{dn}{dt}
$$
where $S$ is the number of active sites.
* Note that $r_t$ is a rate and not a rate constant, so it is always necessary to specify all conditions such as temperature or pressure.

#### Transition state theory
* In this section we take up the question of how DFT calculations can be used to calculated the rates of chemical process.
* Everything is based on the potential energy surface (PES).
* In situations where the energy required for a chemical process is considerably larger than the typical thermal energy, the net rate of the process can be calculated using the transition state theory (TST)
* The fundamental idea of the TST is that the rate constant of the process is (MORE)
$$
\begin{align*}
k = (\text{average thermal velocity of crossing the border}) \\
\times(\text{probability of finding atom at TS})
\end{align*}
$$

* TST can be applied under the following assumptions:
  * When developing a rate theory based on the existence of a PES, we are therefore implicitly involving the Born-Oppenheimer approximation. This approximation assumes that motion of the electrons is instantaneous compared with the motion of the nuclei, such that for whatever motion, we shall never move on an electronically excited state.
  * We shall in addition assume that the rate of quantum tunneling through the potential barriers is negligible compared with the rate obtained from the classical treatment.
  * The last assumption is that once the system attains a TS configuration with a velocity toward the product region, it will necessarily "react" and become the product. For such a configuration to read to a reactive event, the system should not re-enter the reaction region again shortly after moving into the product region.

* The rate at which the system traverses the infinitesimal region of thickness around the TS, $\delta x$, is given by the average velocity orthogonal to the configuration space TS:
$$
r = \frac{|v|}{\delta x}
$$
The average orthogonal velocity can be evaluated directly from the Boltzmann distribution
$$
\begin{split}
&|v| = |{\bf v} \cdot {\bf n}({\bf x})| \\
&= \frac{\int_{TS}dx \int_{-\infty}^{\infty}dv |v\cdot n(x)|
					e^{-V(x)/k_B T}
					e^{-\sum_i\frac{1}{2}m_i v_i^2/k_B T}}
       {\int_{TS}dx\int_{-\infty}^{\infty}dv
			 		e^{-V(x)/k_B T}
					e^{-\sum_i\frac{1}{2}m_i v_i^2/k_B T}} \\
&= \sqrt{\frac{k_B T}{2\pi\mu}}
\end{split}
$$
* The probability of being in the $\delta x$ region compared to the reactant region is given
$$
P = \frac{q^{TS\pm \frac{\delta x}{2}}}{q^R}
  = \frac{\sum_{i} \exp(-\epsilon_i/k_B T)}{\sum_{k} \exp(-\epsilon_k/k_B T)}
$$

$$
E_n = \frac{h^2 n^2}{8\mu \delta x^2}
$$

$$
q^{\delta x} = \sum \exp(-\epsilon_n/k_B T) = \sum \exp\left( -\frac{h^2n^2}{8\mu\delta x^2}/k_B T \right)
$$

$$
q^{\delta x} \approx \int_{n=0}^{\infty} \exp\left( -\frac{h^2n^2}{8\mu\delta x^2} /k_B T \right) dn
= \sqrt{\frac{2\pi\mu k_B T}{h^2}}\cdot \delta x
$$

$$
P = \frac{q^{TS\pm\frac{\delta x}{2}}}{q^R}
  = \frac{q^{\delta x}\cdot q^{TS}}{q^R}
	= \sqrt{\frac{2\pi\mu k_B T}{h^2}} \delta x \frac{q^{TS}}{q^R}
$$

* This gives us the TST rate constant as the product of the rate at which the TS is traversed multiplied by the probability of the system to be found in the TS:
$$
k_{TST} = r\cdot P =
\left( \sqrt{\frac{k_B T}{2\pi\mu}} \frac{1}{\delta x} \right) 
\left( \sqrt{\frac{2\pi\mu k_B T}{h^2}} \delta x \frac{q^{TS}}{q^R} \right)
$$
which reduces to
$$
k_{TST} = \frac{k_B T}{h}\cdot\frac{q^{TS}}{q^R}
$$
* The partition function, $q^{TS}$, is now the partition function for being "in" the TS and not just the $\pm \delta/2$ vicinity of the TS, thus now excluding the mode perpendicular to the TS.
* The relative partition function expression, $q^{TS}/q^R$, is as mentioned earlier an expression for the relative probability of finding the system in the TS compared to finding it in the reactant state.
* Since one of the fundamental assumptions underlying the TST is that the distribution of states in the reactant region and the TS is thermally equilibrated. Thus it is reasonable to set the relative partition functions equal to an equilibrium constant, which is defined through the differences in standard Gibbs energy between the TS and the reactant state:
$$
K_{TST} = \frac{q^{TS}}{q^R} = \exp\left( -\frac{G_{TS}^0 - G_R^0}{k_B T} \right) 
= \exp\left(-\frac{\Delta G_{TS}^0}{k_B T}\right)
$$
This means the overall rate constant has a very simple form:
$$
k_{TST} = \frac{k_B T}{h}\exp\left(-\frac{\Delta G_{TS}^0}{k_B T}\right)
$$
* It is expected that the TS for a unimolecular reaction may have a structure similar to that of the reactant except for bond elongation prior to breakage. If this is the case,
$$
\frac{k_B T}{h} \approx 10^{13} s^{-1}
$$
* It has been experimentally verified that numerous unimolecular reactions have rate constants with pre-exponential factors on the order of $10^{13} s^{-1}$. However, the pre-exponential factor can be either larger or smaller than $10^{13} s^{-1}$ depending on the details of the TS.

#### pre-exponential
* Based on the collision theory, pre-exponential factor for the adsorption is calculated with following equation
$$
A_{ads} = \frac{\sigma^0(T,\theta)P_A A_{site}}{\sqrt{2\pi m_A k_B T}}
$$
where $A_{site}$ is the area of a single active site for adsorption, $m_A$ is the molecular mass. $\sigma^0(T, \theta)$ is the sticking probability.

#### microkinetics
* A catalytic reaction often involves several species and elementary reactions. 
* Microkinetic analysis (or detailed reaction mechanism) treats a set of elementary reactions.
* Therefore, if is often more accurate than the kinetic analysis based on the global rate expression.
* If the thermal fluctuations are small enough compared with the configuration? of the potential energy surface, but large enough that the adsorbates can slowly diffuse around, the adsorbate structuring on the surface will occur.
* This complexity is often difficult to treat exactly, and the solution of such problems with many interacting bodies is a key research area of statistical mechanics. In statistical mechanics, an often-utilized approximation to a solution is the so-called **mean field model**.
* In a mean field model, one replaces all the detailed interactions between any one body and the rest of the system with an average or "effective" interaction.
* This replacement turns a many-body problem into a set of one-body problems.
* In catalysis, the term "mean field microkinetic modeling" commonly refers to the specific microkinetic model in which all repulsive (or attractive) interactions between adsorbates have been removed.

##### Microkinetics of elementary surface processes (CAN BE OMITTED?)
* Let us assume we have determined a rate constant, $k_{-}$, for the desorption of species A from a surface. The rate at which desorption occurs from a given site on the surface is proportional to the product of the rate constant and the probability that a given site is occupied by A. This probability is equal to the fractional coverage of the adsorbate on the surface, $\theta_A$. So,
$$
r_\text{desorption} = k_{-}\theta_{\rm A}
$$
* Likewise, we would expect $r_\text{adsorption}$ to be proportional to the probability that a site is unoccupied, $\theta_{*}$. At equilibrium, this process is described by
$$
r_\text{adsorption} = k_{+}p_{\rm A}\theta_*
$$
* For the elementary reaction
$$
\ce{A + $*$ <=> A$*$}
$$
we can define the reaction fractions by
$$
\left. \frac{\theta_{\rm A}}{p_{\rm A}\theta_*}\right|_{\rm eq} = K_{\rm ads} = \exp\left(-\frac{\Delta G^o}{k_B T}\right)
$$
and at equilibrium,
$$
\gamma = \frac{\frac{\theta_{\rm A}}{p_{\rm A}\theta_*}}{\left.\frac{\theta_{\rm A}}{p_{\rm A}\theta_*}\right|_{\rm eq}}
$$
* It is convenient to define the **approach to equilibrium** (also sometimes called the **reversibility**) for an elementary reaction or for an overall reaction as the rate of the backward reaction divided by the rate of the forward reaction.

##### Microkinetics of several coupled elementary surface processes
* Let us first focus on the simple catalytic cycle, whereby two reactant gases, A2 and B reaction to form the gas AB. We assume that the process proceeds in the two elementary steps,
$$
\begin{align*}
\ce{A2 + 2$*$ &-> 2A$*$}&    \ &(\text{R1})\\
\ce{A$*$ + B  &-> AB + $*$}& \ &(\text{R2})
\end{align*}
$$
and that the rate constants for the elementary reactions have already been determined.
* We can now write up the reaction rate expressions
$$
\begin{align*}
& R_1 = r_1 - r_{-1} = k_1 p_{\rm A_2}\theta_*^2 - k_{-1}\theta_A^2 \\
& R_2 = r_2 - r_{-2} = k_2 p_{\rm B}\theta_{\rm A} - k_{-2}p_{\rm AB}\theta_*
\end{align*}
$$
* In these expressions, $R_1$ describes the net rate of $\ce{A_2}$ removal, and $R_2$ describes the net rate of AB formation. Since both of the two elementary reactions need to occur in order for the product to be formed, and also since step1 changes the coverage that goes into step 2, the two rate expressions must be solved simultaneously in order to obtain the overall reaction rate.
* If we think of the net rates in terms of what surface species they create and remove, we see that reaction 1 creates 2A*, while reaction 2 removes 1A*. This establishes a differential equation for the time development of the coverages of A:
$$
\frac{\partial \theta_{\rm A}}{\partial t} = 2 R_1 - R_2 = 2k_1 p_{\rm A_2}\theta_*^2 - 2 k_{-1}\theta_{\rm A}^2 - k_2 p_{\rm B}\theta_{\rm A} + k_{-2}p_{\rm AB}\theta_*
$$
* Though there are two coverages i.e. $\theta_A$ and $\theta_{*}$, one differential equation is enough to completely specify the reaction, because the fractional coverage should sum to one (*site conservation* or *site balance*), as
$$
\sum_i \theta_i = 1
$$
* Significant reduction in the complexity of solving the microkinetic equations can be achieved if we employ the approximation that the rate of change of all the coverages is zero, that is, the *steady-state approximation*:
$$
\frac{\partial\theta_i}{\partial t} = 0 \ \text{for all i}
$$
* The steady-state approximation effectively turns the microkinetic model from a set of coupled non-linear differential equations in time into a time-independent algebraic root-finding problem, which is simpler to solve.
* For the reactions R1 and R2, the steady-state approximation amounts to
$$
2k_1p_{\rm A_2}(1-\theta_{\rm A}^2) - 2k_{-1}\theta_{\rm A}^2 - k_2p_{\rm B}\theta_{\rm A} + k_{-2}p_{\rm AB}(1-\theta_{\rm A}) = 0
$$
* For many heterogeneous catalytic reactions, the overall reaction rate is typically determined by one specific elementary reaction being particularly "slow".
* Here, a useful concept is that of a *strongly rate-determining reaction step*. By that, we shall mean a reaction step that is so difficult that all the other reaction steps are equilibrated. (Langmuir-Hishelwood-Hougen-Watson?)
* To solve the reaction set (R1 and R2) in the strongly rate-determining approximation, we write the rates in terms of the following reversibilities:
$$
\begin{align*}
R_1 = k_1 p_{\rm A_2} \theta_*^2 (1-\gamma_1),\hspace{5mm}  \gamma_1 = \frac{\theta_{\rm A}^2}{p_{\rm A_2}\theta_*^2}\frac{1}{K_1} \\
R_2 = k_2 p_{\rm B} \theta_{\rm A}(1-\gamma_2),\hspace{5mm} \gamma_2 = \frac{p_{\rm AB}\theta_*}{p_{\rm B}\theta_{\rm A}}\frac{1}{K_2}
\end{align*}
$$
* We note that a general feature for a serial catalytic cycle, which relies on reaction step $i$ being carried out $n_i$ times, the overall equilibrium constant is given by $K_{\rm eq} = \prod_i K_i^{n_i}$ and the reversibility of the overall reaction is $\gamma = \prod_i\gamma_i^{n_i}$. For the reaction in question, we have $n_1 = 1$ and $n_2 = 2$, so the overall equilibrium constant is $K_{\rm eq} = K_1 \cdot K_2^2$, and the overall reversibility is $\gamma_{\rm eq} = \gamma_1\cdot\gamma_2^2$.
* Let us now assume that R1 is the rate-determining step. This means that $\gamma_2 = 1$ and thus $\gamma_1 = \gamma_{\rm eq}$.
* From Eq.X, we get the coverage of free sites as
$$
\theta_{*} = \frac{1}{1+p_{\rm AB}p_{\rm B}^{-1}K_2^{-1}}
$$
and the coverage of A as
$$
\theta_{A} = \frac{p_{\rm AB}p_{\rm B}^{-1}K_2^{-1}}{1+p_{\rm AB}p_{\rm B}^{-1}K_2^{-1}}
$$
* The rate of reaction then found by taking the rate of the rate-determining step (R1) as the rate of R2 is zero. By inserting the appropriate coverages and using the reversibility of the overall reaction ($\ce{A2 + 2B -> 2AB}$)
$$
\gamma = \frac{p_{\rm AB}^2}{p_{\rm A_2}p_{\rm B}^2}\frac{1}{K_{\rm eq}} = \frac{p_{\rm AB}^2}{p_{\rm A_2}p_{\rm B}^2 K_1K_2^2}
$$
one finally has the rate of the overall reaction as
$$
R = k_1p_{\rm A_2}\theta_{*}^2(1-\gamma) = k_1p_{\rm A_2}\frac{1}{(1+p_{\rm AB}p_{\rm B}^{-1}K_2^{-1})^2}\left(1-\frac{p_{\rm AB}^2}{p_{\rm A_2}p_{\rm B}^2 K_1K_2^2}\right)
$$

#### Activity map
##### one-dimensional
* From the scaling relation between $\Delta E$ and $E_a$, there is only one independent variable, which we choose to be $\Delta E$.
* Based on this feature, the reaction rate R or turn-over-frequency (TOF) can be evaluated as the function of $\Delta E$, i.e. $R = f(\Delta E)$.
* Here, following assumptions are made
  * the reaction energy of the overall reaction is set to $-0.3\ eV$
  * the entropies of all the gas-phase species (A2, AB, and B) are set to $0.001\ eV/K$
  * the entropy of A adsorbed on the surface is assumed to be negligible
  * the pre-exponential factor for the dissociative adsorption of A2 (R1) is set to $k_B T/h$
  * dissociation of A2 follows a transition state scaling relationship of $\alpha = 0.87$ and $\beta = 1.34\ eV$.
* The dependence of $k_1$ and $\theta_*$ on $\Delta E$ are shown in Fig.X.
* Since the reaction rate or TOF depends on the product of $k_1$ and $\theta_*$, it also depends on $\Delta E$, as shown in the Figure.
* By decreasing the $\Delta E$, the number of free sites $\theta_*$ decreases dur to the poisoning of the surface by species $A*$, while $k_1$ increases because it leads to the smaller Ea value.
* As a result, the maximum TOF is obtained at $\theta_* ~ 0.5?$.
* This represents the **Sabatier principle**, which states that the catalytic activity for a given reaction follows a volcano-shape curve, because only an intermediate binding of intermediates on the surface of a catalyst will give a reasonably active catalyst.

##### two-dimensional
* Scaling relations always provide a way to reduce the number of independent variables characterizing the catalyst, but in most cases, more than one descriptor is necessary.
* An example where two descriptors are necessary is the CO methanation reaction (or Sabatier reaction)
$$
\ce{CO + 3H2 -> CH4 + H2O}
$$
* The elementary reactions should be
$$
\begin{align*}
\ce{should} \\
\ce{see} \\
\ce{the paper}
\end{align*}
$$
* The reaction rate is given by
$$
R = k_2 \theta_\ce{CO}\theta_*(1-\gamma)
$$
* The activity map for the CO hydrogenation is shown in Figure X. It shows a single maximum for $(\Delta E_C, \Delta E_O) \approx (0.5\ eV, -3.0 \ eV)$, thus Ru and Co is the transition metal elements closest to the maximum. This agrees with the experimental trends.
* Knowing the optimum value for ($\Delta E_C, \Delta E_O$) allows a search for other catalysts for this process.
* SUMMARY: to draw the activity map, one should (i) express the overall reaction rate by some descriptors, (ii) draw the rate by the function of selected descriptors. The number of descriptors is equal to the dimension of the activity map.

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

## 4. Theoretical Background of Electrocatalysis
#### basic
* The reversible potential, also called equilibrium potential, is a voltage at which there is no net ion flux.

* The potential difference between the anode and the cathode gives rise to a variation in the electrostatic potential through the cell. Since the electrolyte is conducting, there is no electrical field there and the potential is constant. The potential variation happens in the so-called dipole layer. They are formed near the two electrodes and sets up strong electrical fields there. The field is set up by the electrons in the electrode (or holses for a positive electrode) and the counterions in the electrolyte.
* We will discuss three? important ways in which elementary surface electrochemical reactions may differ from their gas-phase counterparts:
	1. The chemical potential of the electrons entering the reaction is controlled by the potential of the electrode. This is by far the largest effect. Changing the potential by 1 V changes the reaction free energy by 1 eV, when one $e^{-}$ is involved, for example.
	2. The surface species will be solvated by the electrolyte. This effect is notably larger for water due to the hydrogen bonding network.
	3. The electronic field at the solid-electrolyte interface will change the adsorption energy.

* Since the energy of the electron entering is now potential-dependent, it means that the potential energy diagram becomes potential-dependent.
* The energy-scaling relationship also holds in the electrochemical reaction. This means that one can tune the activation energy by shifting the reaction energy, as
$$
E_a = \gamma \Delta E + \xi
$$
* Since the reaction energy $\Delta E$ depends on the potential as
$$
\Delta E = \Delta E (U=0) + eU = \Delta E_0 + eU
$$
the activation energy depends on the potential as
$$
E_a(U) = \gamma(\Delta E_0 + eU) + \xi = \gamma eU + E_{a0}
$$
If we neglect variations in the coverage of intermediates with potential, the rate of elementary step will vary with potential as
$$
\begin{split}
r(U) &= \nu \exp(-E_a(U)/k_B T)  \\
     &= \nu \exp(-E_{a0}/k_BT) \exp(-\gamma eU / k_B T) \\
     &= A \exp(-\gamma eU / k_B T)
\end{split}
$$
which is the Tafel equation. A plot of the logarithm of the rate versus U should give a linear plot, known as the **Tafel plot**, with slope $\gamma$.
* One can calculated the **polarization curve** (current density vs. potential) from the theoretical PED.[Hansen, JPC.C, 2014, 118, 6706]

#### example: OER/ORR
* Since energy barriers are quite small for proton transfer reactions and since the barrier scales with the reaction energy, it is often useful to consider simplified free energy diagrams where only the energies of intermediates are included.
* The potential dependence of the reaction is easily seen by showing the variation for the free energy diagram with potential.
* For all the electrocatalysts, the OER does not start to have an applicable current density $j$ until a substantial overpotential
$$
\eta(j) = |U(j) - U_{eq}|
$$
is applied.
* Typically, the best OER catalysts ($\ce{RuO2}$ and $\ce{IrO2}$) and ORR catalyst (Pt) only give current densities above $|j| = 5 mA/cm^2$ at potential that are $|\eta| \approx 0.3 V$ above or below the equilibrium potential, $U_{eq} = 1.23 V$.

* The OER/ORR is considered possible if and only if all the involved reaction steps are neutral or downhill in the Gibbs free energy. For a given reaction we can determine the lowest potential at which this is the case. This is the limiting potential.
* If one neglect the barrier, this would be a lower bound to the overpotential.

* We can write the free energy change for any elementary step $i$ iin an electrochemical reaction with a transfer of one electron and a proton as 
$$
\Delta G_i(U) = \Delta G_{0,i} \pm eU
$$
where the sign of the last term depends on whether the electron transfer is from or to the surface.
* We can now define the limiting potential for elementary reaction step $i$ as the potential where the free energy difference for the reaction is zero:
$$
U_{L,i} = \frac{\mp\Delta G_{0,i}}{e}
$$
* For the ORR, the electrons are transferred from the surface to the reactant (thus plus sign). The minimum of the U for the elementary steps defines the potential where all steps are exergonic, and this potential is termed the **limiting potential** for the reaction
$$
U_{L,red} = \min(U_{L,i})
$$
* We define the theoretical overpotential as
$$
\eta_{red} = U_{eq} - U_{L,red}
$$
* For ORR,
$$
U_{L,ox} = \max(U_{L,i})
$$
and the theoretical overpotential is
$$
\eta_{ox} = U_{L,ox} - U_{eq}
$$

#### Nernst equation
* The relationship between the electrode potential and the activity is expressed by the **Nernst equation**. For the reaction $\ce{aA + bB <=> cC + dD}$,
$$
\begin{split}
E &= E^0 - \frac{k_B T}{ne}\ln K \\
  &= E^0 - \frac{k_B T}{ne}\ln \frac{a_C^c a_D^d}{a_A^a a_B^b}
\end{split}
$$
at standard temperature,
$$
E = E^0 - \frac{0.059}{n}\log\frac{a_C^c a_D^d}{a_A^a a_B^b}
$$
* when the reaction involves proton, $E$ can be expressed as the function of $pH$ by using 
$$
pH = -\log_{10}a_{H^+}
$$
* If the electrode potential is raised above the equilibrium potential, oxidation is preferred, while reduction is preferred if below the equilibrium potential.
* use $\ln(10) = 2.303$ to convert $\ln$ to $\log$

#### Pourbaix diagrams
* As shown in Section X, ab initio thermodynamics is able to evaluate the stability of the surface structure under reaction conditions as a function of the temperature and the reactant partial pressure. In electrocatalysis, rather than temperature and pressure, the pH value and the electrostatic potential ($U$) are extra variables, and these are actually more important variables to control the electrochemical reaction. Thus, the surface phase diagrams depending on these variables are desirable. In this section, we will see how to construct such diagram by ab initio method.
* Potential-pH diagrams are usually called **Pourbaix diagram**, after the name of their originator Marcel Pourbiax (1904-1998), a Russian-born Belgian chemist. In such diagrams, the redox potential is plotted on a vertical axis and the pH on a horizontal axis. Pourbaix diagrams show the most stable bulk phase of an element as function of pH and potential.

* Water splitting potentials in acid and alkaline electrolytes are known from cyclic voltammetry. However, the molecular structure of the surface is not known. It is obvious that the molecular structures on the surface is dependent on the metal, potential, and pH of the electrolyte.
* These diagrams constructed from the calculation based on Nernst equation.
* A horizontal line represents a reaction that does not involve pH: that is, neither $\ce{H+}$ or $\ce{OH-}$ is involved.
* A vertical line involves $\ce{H+}$ or $\ce{OH-}$ but not electrons
* A sloping line involve $\ce{H+}$, $\ce{OH-}$, and electrons.

##### e.g. water electrolysis
* anode reaction: $\ce{H2 <=> 2H+ + 2e-}$
* cathode reaction: $\ce{1/2O2 + 2H+ + 2e- <=> H2O}$
$$
\begin{split}
E_1 &= E^0 - 2.303\frac{RT}{zF}\log\frac{1}{[H^+]^2} = - 0.059\times pH \\
E_2 &= E^0 - 2.303\frac{RT}{zF}\log\frac{1}{[H^+]^2} = 1.229 - 0.059\times pH
\end{split}
$$
* Below the E_1 line, the reduction is preferred in the anode reaction thus H2 is stable. Above it, H+ is formed.
* Similarly in E_2 line, below it the H2O is stable but above it O2 is generated.

* **Surface Pourbaix diagrams** predict the thermodynamically most stable adsorbed species and thier coverage as a function of pH and the potential.
* The Pourbaix diagram can be also made for the molecular systems, which is interpreted as the phase diagram that shows the pH and potential at which certain molecular species are energetically stable.

##### ORR example
* Hansen et al. have shown that the ab initio surface Pourbaix diagram is quite useful in analyzing the oxygen reduction reaction (ORR) by metals such as Ag, Pt, or Ni.[PCCP]
$$
\ce{O2 + 4H+ + 4e- -> 2H2O}
$$
* Here it is assumed that the surface is in equilibrium with protons and liquid water at 298 K, so that O and OH may be exchanged between the surface and a reference electrolyte through the following steps
$$
\ce{H2O(l) + $*$ <=> OH$*$ + H+(aq) + e-} \\
\ce{OH$*$ <=> O$*$ + H+(aq) + e-}
$$
* When we assume the computational hydrogen electrode (CHE), the H2 evolution reaction becomes equilibrium at zero pH and zero potential vs. SHE ($U_{\rm SHE}$);
$$
\ce{H+(aq) + e- <=> 1/2H2(g)}
$$
* Using this, the above equations become
$$
\ce{H2O(l) + $*$ <=> OH$*$ + 1/2H2(g)} \\
\ce{OH$*$ <=> O$*$ + 1/2H2(g)}
$$
* The Gibbs energy change along the H2 evolution (Eq.X) depends on pH and potential as
$$
\Delta G = e U_{\rm SHE} + k_B T \ln10 {\rm pH}
$$
* In this case, adsorption Gibbs energy of OH* is, for example, 
$$
G(\ce{OH$*$}) = \Delta G_0(\ce{OH$*$}) - e U_{\rm SHE} - k_B T \ln10 {\rm pH}
$$
where $\Delta G_0$ is the reaction Gibbs energy of Eq.X at zero pH and zero $U_{\rm SHE}$. $\Delta G_{\rm field}$ is the change in the adsorption energy due to the electric field in the electrochemical double layer. This term can be neglected.

* Among the terms in Eq.X, one should evaluate $\Delta G_0$ i.e. the Gibbs energy of adsorption at standard condition. This term can be calculated with ab initio method, as
$$
\Delta G_0 = \Delta E + \Delta ZPE - T\Delta S
$$
where $\Delta E$ is the total energy $\Delta ZPE$ is the change in zero-point energy of the adsorbate. The entropy change along the adsorption $\Delta S$ is approximated from the loss of entropy of the gas phase molecules upon binding them to the surface.
* Thus, in summary, the Pourbaix diagram for the ORR can be constructed by ab initio method likewise the ab initio thermodynamics diagram (Section X). In the Pourbaix diagram case, the variables in the Gibbs energy of adsorption is $U_{\rm SHE}$ and pH, as seen in Eq.X. We can then construct the diagram along these two variables.

* We can start from the one-dimensional version of Pourbaix diagram by fixing pH. See Fig.X, in which the Gibbs energy of O* and OH* on Ag(111) at pH = 0 are plotted as a function of potential $U_{\rm SHE}$ and pH. The Gibbs energy dependence on these variables are evaluated by Eq.X. In the plot, the lowest line determines the surface with the lowest Gibbs energy at a given potential.
* The usual (or two-dimensional) Pourbaix diagram is readily plotted by extending this plot by evaluating the Gibbs energy by $U_{\rm SHE}$ and pH; see Fig.X.
* The potential for reversible water oxidation will show up as lines with a slope of $-k_B T/e\ln 10$ (or $-0.059*pH (\text{in V, at 298 K})$) in the Pourbaix diagram.
* From this Pourbaix diagram, following features are found;
  * at acidic condition, dissolution of the Ag electrode occurs at $U_{\rm SHE} > 0.8\ V$.
  * at moderately acidic condition, water oxidation occurs before the Ag dissolution.
  * at basic condition, the substitution of Ag atoms by O atoms are observed. This can be seen as the onset of the Ag electrode oxidation.
  * at pH = 14, water oxidation occurs at $U_{\rm SHE} = 0.10\ V$ (or $U_{\rm RHE} = 0.93\ V$). This value can be compared with cyclic voltammogram experiment, as confirmed in their work.

* The difference between the potential vs. SHE ($U_{\rm SHE}$)or the potential vs. reversible hydrogen electrode (RHE) ($U_{\rm RHE}$) is explained in Section X. The RHE easier to use in experiments thus more often used, and their values are related as follows;
$$
U_{\rm RHE} = U_{\rm SHE} + \frac{k_B T}{e} \ln10 
$$

## 5. Modelling Catalytic Systems
### 5-1. Modeling Heterogeneous Systems
#### basic
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
#### solution

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
