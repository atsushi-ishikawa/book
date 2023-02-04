## 4. Theoretical Background of Electrocatalysis
#### Basic
* The reversible potential, also called equilibrium potential, is a voltage at which there is no net ion flux.

<!-- Norskov begin -->
* The potential difference between the anode and the cathode gives rise to a variation in the electrostatic potential through the cell. Since the electrolyte is conducting, there is no electric field there and the potential is constant. The potential variation happens in the so-called dipole layer. They are formed near the two electrodes and sets up strong electric fields there. The field is due to the electrons in the electrode (or holses for a positive electrode) and the counterions in the electrolyte.
* We will discuss three? important ways in which elementary surface electrochemical reactions may differ from their gas-phase counterparts:
	1. The chemical potential of the electrons entering the reaction is controlled by the potential of the electrode. This is by far the largest effect. Changing the potential by 1 V changes the reaction free energy by 1 eV, when one $e^{-}$ is involved, for example.
	2. The surface species will be solvated by the electrolyte. This effect is notably larger for water due to the hydrogen bonding network.
	3. The electronic field at the solid-electrolyte interface will change the adsorption energy.

* For elementary processes involving charge transfer at the surface, one can define a potential energy diagram for the process completely as for other surface processes.
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
which is the Tafel equation. A plot of the logarithm of the rate versus $U$ should give a linear plot, known as the **Tafel plot**, with the slope of $\gamma$.
* One can calculated the **polarization curve** (current density vs. potential) from the theoretical PED.[Hansen, JPC.C, 2014, 118, 6706]

#### Butler-Volmer, exchange current density, overpotential
* asdf

#### The overpotential in electrocatalytic processes
* Free energy diagrams for full electrocatalytic reactions are as useful in understanding surface electrocatalysis as those introduced for heterogeneous catalysis.
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

* We can write the free energy change for any elementary step $i$ in an electrochemical reaction with a transfer of one electron and a proton as 
$$
\Delta G_i(U) = \Delta G_{0,i} \pm eU
$$
where the sign of the last term depends on whether the electron transfer is from or to the surface.
* We can now define the limiting potential for elementary reaction step $i$ as the potential where the free energy difference for the reaction is zero:
$$
U_{L,i} = \frac{\mp\Delta G_{0,i}}{e}
$$
* For the ORR, the electrons are transferred from the surface to the reactant (thus plus sign). The minimum of the U for the elementary steps defines the potential where all steps are exergonic, and this potential is termed the *limiting potential* for the reaction
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
<!-- Norskov end -->

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

##### Example: Water electrolysis
* anode reaction: $\ce{H2 <=> 2H+ + 2e-}$
* cathode reaction: $\ce{1/2O2 + 2H+ + 2e- <=> H2O}$
$$
\begin{split}
E_1 &= E^0 - 2.303\frac{RT}{zF}\log\frac{1}{[H^+]^2} = - 0.059\times pH \\
E_2 &= E^0 - 2.303\frac{RT}{zF}\log\frac{1}{[H^+]^2} = 1.229 - 0.059\times pH
\end{split}
$$
* Below the $E_1$ line, the reduction is preferred in the anode reaction thus H2 is stable. Above it, H+ is formed.
* Similarly in the $E_2$ line, below it the H2O is stable but above it O2 is generated.

* **Surface Pourbaix diagrams** predict the thermodynamically most stable adsorbed species and thier coverage as a function of pH and the potential.
* The Pourbaix diagram can be also made for the molecular systems, which is interpreted as the phase diagram that shows the pH and potential at which certain molecular species are energetically stable.

##### Example: Oxygen reduction reaction (ORR)
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
* The potential for reversible water oxidation will show up as lines with a slope of $-k_B T/e\ln 10$ (or $-0.059 pH (\text{in V, at 298 K})$) in the Pourbaix diagram.
* From this Pourbaix diagram, following features are found;
  * at acidic condition, dissolution of the Ag electrode occurs at $U_{\rm SHE} > 0.8\ V$.
  * at moderately acidic condition, water oxidation occurs before the Ag dissolution.
  * at basic condition, the substitution of Ag atoms by O atoms are observed. This can be seen as the onset of the Ag electrode oxidation.
  * at pH = 14, water oxidation occurs at $U_{\rm SHE} = 0.10\ V$ (or $U_{\rm RHE} = 0.93\ V$). This value can be compared with cyclic voltammogram experiment, as confirmed in their work.

* The difference between the potential vs. SHE ($U_{\rm SHE}$)or the potential vs. reversible hydrogen electrode (RHE) ($U_{\rm RHE}$) is explained in Section X. The RHE easier to use in experiments thus more often used, and their values are related as follows;
$$
U_{\rm RHE} = U_{\rm SHE} + \frac{k_B T}{e} \ln10 
$$
