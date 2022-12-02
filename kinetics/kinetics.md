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
