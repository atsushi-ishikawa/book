### 3-2. Computational Fluid Dynamics

* One of these characteristics is the mixing of fluids, which was addressed for both static mixing elements and more complex systems.
* Additionally, CFD allows to predict heat transfer and the outcome of reactions in chemical applications.
* CFD can provide fruitful information related to the formation of hot spots in an exothermic reaction environment.

#### Introduction
* Computational fluid dynamics (CFD) calculations have become increasingly important for obtaining insight into characteristics of flow reactors in particular, which are difficult to obtain experimentally. The basic equations in the CFD calculations are the conservation equations, which have been well known for over a century. However, these equations are inherently nonlinear, making exact solutions generally unavailable. There are certainly some fluid-mechanics problems that can be solved exactly, but these analytic solutions demand approximations such as constant properties. When chemical reactions are included, exaction solutions are essentially impossible. The CFD numerically solves these equations with computers, and makes quantitative analysis and/or visualization of the chemically complex three-dimensional problems.
* The present book considers only laminar flow, thus completely neglecting turbulent flow. The applications are generally concerned with internal flows with small  spatial dimensions (i.e. low Reynolds numbers). Many chemical reactors fall into this regime. Small-scale heat exchangers and fuel cells are based on the laminar flow in relatively small channels. Chemical reactors typically rely on catalysts, either in packed beds or as channel washcoats.

* Any material that deforms continuously under the influence of shearing forces is called a fluid.
* Generally speaking, a fluid can be a liquid or a gas, where an important difference is the fluid density, which plays a central role in fluid mechanics. For gases, however, one needs an equitation of state to establish the essential relationship among density, temperature, pressure, and species composition. A liquid is an incompressible fluid, meaning that the density can be treated as a constant with respect to the pressure and temperature changes.
* In addition to the equation of state, it is necessary to describe other thermodynamic properties of the fluid, such as specific heat, enthalpy, entropy, or free energy.
* Most descriptions of fluid behavior also depends on transport properties, including viscosity, thermal conductivity, and diffusion coefficients. These properties generally depend on temperature, pressure, and mixture composition.
* This book considers only fluids that are isotropic, meaning that the fluid properties are independent of direction. Such fluids are known as Newtonian fluids.

##### velocity
* Because a fluid is continuously deformable, defining its velocity takes a care. At the smallest scale the fluid is an ensemble of molecules. In principle, one can describe the velocity of a fluid in terms of the velocities of each molecule in the fluid. Obviously, this would be impractical owing to the extreme numbers of molecules that have to be considered. Instead, it is appropriate to use a velocity field that represents an average fluid velocity at very point within a macroscopic fluid domain.

#### Fluid properties

##### Ideal gas
* The ideal-gas equation of state for multicomponent mixtures depends on the species composition. Representing the composition as either mass fraction $Y_k$ or mole fraction $X_k$ leads to
$$
\rho = \frac{p}{RT}\frac{1}{\sum Y_k / W_k} = \frac{p}{RT}\sum_{k=1}^{K}X_k W_k = \frac{p}{RT}\bar{W}
$$
where $\rho$ is the mass density ($kg m^{-3}$), $p$ is the pressure ($N m^{-2}$), $T$ is the temperature ($K$). The universal gas constant is $R = 8.314 J mol^{-1} K^{-1}$. The species molecular mass are $W_k$ ($kg mol^{-1}$) and
$$
\bar{W} = \frac{1}{\sum_{k=1}^K Y_k / W_k} = \sum_{k=1}^K X_K W_k
$$
is the average molecular mass.

##### Thermodynamics
* The first law of thermodynamics for a system, which forms the fundamental basis for energy balances, may be stated as
$$
\frac{d E_t}{dt} = \frac{dW}{dt} + \frac{dQ}{dt}
$$
where $E_t$ is the total internal energy ($J$), $dW/dt$ represents the work *done on the system*, and $dQ/dt$ is the heat *added to the system*. The total internal energy may include kinetic and potential energy.
* The first law is often more conveniently represented on a per-unit-mass basis. The specific internal energy and enthalpy are $e = E/m$ and $h = H/m$ ($J kg^{-1}$), respectively. In this case, $h = e + pv$, where $v$ is the specific volume ($m^3 kg^{-1}$).
* It is often useful to represent the thermodynamics properties on a molar basis. Assuming a species with molecular mass $W$, the molar properties ($J mol^{-1}$) are
$$
E = eW, H = hW, S = sW, G = gW, A = aW
$$

* The thermodynamics properties of ideal gases have significant temperature dependencies.
* Broadly speaking, gas-phase heat capacity increases as a function of increasing temperature.
* Although a molecule's heat capacity depends on translational, rotational, and vibrational contributions, many CDF softwares typically represent the heat capacities as fits rather than directly evaluating them from fundamental theories.
* A widely used polynomial fit was first introduced in the NASA chemical equilibrium program [Ref]. In this case, the standard-state heat capacity $C_p^{\rm o}$ ($J mol^{-1} K^{-1}$) is represented as
$$
\frac{C_p^{\rm o}}{R} = \sum_{n=1}^N a_n T^{n-1}
$$
where $R$ is the gas constant and $T$ is temperature. The NASA program and databases used $N = 5$.
* The standard-state molar enthalpy ($J mol^{-1}$) is related to the standard-state heat capacity as
$$
H^{\rm o} = \int_{298}^{T}C_p^{\rm o}dT + H^{\rm o}(298)
$$
where $H^{\rm o}(298)$ is the standard-state enthalpy at $T = 298 K$. Assuming heat-capacity fits in the form of Eq.X, enthalpy fits follow as
$$
\frac{H^{\rm o}}{RT} = \sum_{n=1}^N \frac{a_n T^{n-1}}{n} + \frac{a_{N+1}}{T}
$$
where $a_{N+1}$ represents $H^{\rm o}(298)/R$.
* The standard-state molar entropy ($J mol^{-1} K^{-1}$) is related to the heat capacity as
$$
S^{\rm o} = \int_{298}^{T}\frac{C_p^{\rm o}}{T}dT + S^{\rm o}(298)
$$
where $S^{\rm o}(298)$ is the standard-state entropy at $T = 298K$. Again, assuming heat-capacity fits (Eq.X)
$$
\frac{S^{\rm o}}{R} = a_1 \ln T + \sum_{n=2}^N \frac{a_n T^{n-1}}{n-1} + a_{N+2}
$$
where $a_{N+2}$ represents $S^{\rm o}(298)/R$.
* Widely used databases catalog the thermodynamic property fits as
$$
\begin{align}
\frac{C_p^{\rm o}}{R} &= a_1 + a_2 T + a_3 T^2 + a_4 T^3 + a_5 T^4 \\
\frac{H^{\rm o}}{RT}  &= a_1 + \frac{a_2}{2}T + \frac{a_3}{3}T^2 
                             + \frac{a_4}{4}T^3 + \frac{a_5}{5}T^4 + \frac{a_6}{T} \\
\frac{S^{\rm o}}{R}   &= a_1 \ln T + a_2 T + \frac{a_3}{2}T^2 + \frac{a_4}{3}T^3 
                             + \frac{a_5}{4}T^4 + a_7
\end{align}
$$
* The fitting coefficients ($a_i (i=1\cdots7)$) are tabulated in the database.
* The fits are given over two temperature ranges (high- and low-temperature ranges? CHECK), and the properties are made smooth at the interfaces between the two temperature ranges.
* With heat capacities, enthalpies, and entropies in hand, all other ideal-gas thermodynamics properties follow in the usual ways as
$$
\begin{align}
C_v^{\rm o} &= C_p^{\rm o} - R \\
U^{\rm o} &= H^{\rm o} - RT \\
G^{\rm o} &= H^{\rm o} - TS^{\rm o} \\
A^{\rm o} &= U^{\rm o} - TS^{\rm o}
\end{align}
$$
where $C_v^{\rm o}$, $U^{\rm o}$, $G^{\rm o}$ and $A^{\rm o}$ are the molar constant-volume heat capacity, internal energy, Gibbs free energy, and Helmholtz energy, respectively.

#### Transport properties
* Chemically reacting flow must be generally concerned with diffusive molecular-transport processes. Specifically, properties of particular interest include dynamics viscosity $\mu$, thermal conductivity $\lambda$, and diffusion coefficients $D_{ij}$.
* There is an analytic expression fo the viscosity for dilute pure-species, based on the **Chapman-Enskog theory**. This is not covered in this book.
* Another approach is using somewhat empirical formula, where fitting coefficients are used. Gas viscosities are expressed with power law, as
$$
\frac{\mu}{\mu_{\rm o}} = \left( \frac{T}{T_{\rm o}} \right)^n
$$
Here viscosity is measured at a reference temperature $T_{\rm o}$, and the temperature dependence is given in terms of the exponent $n$. The value of $n$ is expected to be less than one; for air, $\mu_{\rm o} = 1.716\times 10^{-5} N s m^{-2}$ at $T_{\rm o} = 273 K$ and $n \approx 0.67$.
* The **Sutherland law** is also widely used to express the temperature dependence of viscosity
$$
\frac{\mu}{\mu_{\rm o}} = \left( \frac{T}{T_{\rm o}} \right)\frac{T_{\rm o}+S}{T+S}
$$
The reference viscosity $\mu_{\rm o}$ is measured for a specific gas at reference temperature, and the parameter $S$ is fit over an appropriate temperature range. For air, using $\mu_{\rm o} = 1.716\times 10^{-5} N s m^{-2}$ at $T_{\rm o} = 273 K$ and $S = 111 K$ provides accurate results over a large temperature range, $170 K < T < 1900 K$.
* Reacting flow inevitably involves mixture of fluids. Thus, mixing rules are needed to evaluate mixture properties from the individual-species properties. Mixture viscosity is often evaluated using the **Wilke formula** as
$$
\mu_{\rm mix} = \sum_{k=1}^K \frac{X_k\mu_k}{\sum_{j=1}^K X_j \Phi_{kj}}
$$
where
$$
\Phi_{kj} = \frac{1}{\sqrt{8}}\left(1+\frac{W_k}{W_j}\right)^{-1/2}
            \left[1+\left(\frac{\mu_k}{\mu_j}\right)^{1/2}\left(\frac{W_j}{W_k}\right)^{1/4}\right]^2
$$
In these expressions, $X_k$ are mole fractions, $\mu_k$ are the individual species viscosities, and $W_k$ are the molecular mass.

#### Diffusion coefficients

#### Thermal conductivity
* The diffusive transport of heat obeys the **Fourier's law** as
$$
{\bf q}'' = –\lambda \nabla T
$$
where ${\bf q}''$ represents the heat flux ($W m^{-2}$).
* Just as diffusive momentum transfer depends on a transport property of the fluid called viscosity $\mu$, diffusive heat transfer depends on a transport property called thermal conductivity $\lambda$.
* Even though there are a number of theories for estimating thermal conductivities, in practice, empirical curve fits are often used, such as power law
$$
\frac{\lambda}{\lambda_{\rm o}} = \left(\frac{T}{T_{\rm o}}\right)^n
$$
, the Sutherland form
$$
\frac{\lambda}{\lambda_{\rm o}} = \left(\frac{T}{T_{\rm o}}\right)\frac{T_{\rm o}+S}{T+S}
$$
or some other empirical form. In the Sutherland form, the value of $S$ is established empirically for each fluid. Assuming that the conductivity $\lambda_{\rm o}$ is known at some reference temperature $T_{\rm o}$, the conductivity at other temperature is easily evaluated.
* As with viscosity and other properties, the thermal conductivity of a fluid mixture must be derived in terms of the individual species conductivities. A commonly used averaging formula is expressed as
$$
\lambda_{\rm mix} = \frac{1}{2}\left(
                    \sum_{k=1}^K X_k\lambda_k + \frac{1}{\sum_{k=1}^K X_k/\lambda_k}
                    \right)
$$
where $X_k$ are mole fractions and $\lambda_k$ are individual species conductivities.

### Fluid kinematics
#### System and control volume
* The study of fluid mechanics is facilitated by understating and using the relationship between a *system* and a *control volume*.
* By definition, a system is a *certain mass of fluid*, that can move about in space. Moreover, the system is free to deform as it moves. As a result, it is practically impossible to follow and account for a particular mass of fluid in a flowing process.
* Nevertheless, because many of the basic physical laws are written in terms of a system (e.g. ${\bf F}=m{\bf a}$), it is convenient and traditional to take advantage of the notion of a system.
* An Eulerian control volume is a fixed region of space. Fluid may flow through the surfaces of the control volume (the control surfaces), carrying with it mass, momentum, energy, and chemical species. Equally important, momentum, energy, and chemical species can "diffuse" across the control surfaces, into and out of the control volume. There can also be creating or destruction of thermal energy and chemical species within a control volume.
* In deriving the conservation laws, it is useful to convert between the system and control-volume views, using both to advantage.
* (Lagrangian picture vs. Eulerian picture?)

#### Extensive and intensive variables (-> thermodynamics?)
* For a system, namely a uniquely identified mass of fluid, it is often appropriate to think of variables or properties that characterize the system as a whole.
* The total mass, momentum, or energy of the system are called *extensive* variables or properties.
* It is reasonable to expect that within a system there may be local spatial variation in variables or properties. The total system property is determined by integrating local distributions over the mass of the system. To accomplish the integration, it is useful to define an *intensive* variable, which is the extensive variable per unit mass. That is, if the extensive variable is called $N$, then the associated intensive variable $\eta$ is defined as
$$
\eta = \frac{N}{m}
$$
where $m$ is the mass.
* For our purposes, it is useful to integrate over a volume that encompasses the system at an instant of time. In this case, the mass density $\rho$ is used. The extensive property of a system is thus given as
$$
N_\text{system} = \int_\text{mass of system}\eta dm = \int_\text{volume of system}\rho\eta dV
$$
* To make this concept concrete, consider a few familiar examples. If $N$ is the mass of a system $m$, then $\eta = 1$; if $N$ is momentum vector ${\bf P}$, then $\eta = {\bf V}$, the velocity vector; and if $N$ is energy $E$ ($J$), then $\eta = e$ the specific internal energy ($J kg^{-1}$).

#### Reynolds transport theorem
* Consider the system and control volume as illustrated in Fig.X. The Eulerian control volume is fixed in an inertial reference frame, described by three independent, orthogonal coordinates, such as $x,y,z$ (Cartesian coordinate) or $z, r, \theta$ (cylindrical coordinate).
* At some initial time $t_0$, the system is defined to contain all the mass in the control volume. A flow field, described by the velocity vector ${\bf V}(t,z,r,\theta)$, carries the system's mass out of the control volume. In the limit of a vanishingly small time $\Delta t$, the relationship between the system and the control volume is know as the *Reynolds transport theorem*.
* The right-hand panel of Fig.X shows the control volume (dashed lines) in its original, fixed position, but the system has partially flowed out of the control volume.
* The figure identifies three regions as
  * Region I:   the volume of control volume that has been vacated by the system
  * Region II:  the volume of the control volume that is still occupied by some of the system mass
  * Region III: the portion of the system mass that has flowed out of the control volume
* The right-hand panel also indicates normal outward-pointing unit vector ${\bf n}$ that describe the local shape of the control surface. Since the control volume retains fixed in space, the ${\bf n}$ vectors also remain fixed in the reference frame.
* The subsequent derivation of partial differential equations that represent basic conservation laws (e.g. conservation of mass, momentum, and energy) are structured around a fixed differential control volume, meaning an *Eulerian* framework.
* The time rate of change of an extensive property $N$ of a system can be written quite generally as
$$
\frac{dN}{dt} = \lim_{\Delta t \rightarrow 0}\frac{N_{t_0+\Delta t} - N_0}{\Delta t}
$$
where $N_{t_0}$ represents the value of $N$ at tome time $t_0$ and $\Delta t$ is some small interval of time. By definition, the system fully occupies the control volume at $t_0$. In other words, the system's and the control volume's extensive property are the same at $t_0$.
$$
N_{t_0} = N_{CV, t_0}
$$
* At $t_0 + \Delta t$, the extensive property of the *system* can be written in terms of the three regions identified in Fig.X as
$$
\begin{align*}
N_{t_0+\Delta t} = N_{\rm II} + N_{\rm III} 
= (N_{CV} - N_{\rm I} + N_{\rm III})_{t_0+\Delta t} \\
= N_{CV,t_0+\Delta t} - N_{{\rm I},t_0+\Delta t} + N_{{\rm III},t_0+\Delta t}
\end{align*}
$$
* Substituting this into Eq.X yields
$$
\left(\frac{dN}{dt}\right)_{\rm system} = \lim_{\Delta t \rightarrow 0}\frac{N_{CV, t_0+\Delta t} - N_{{\rm I},t_0+\Delta t} + N_{{\rm III},t_0+\Delta t} - N_{CV, t_0}}{\Delta t}
$$
* Recognizing that the limit of a sum can be represented as the sum of the limits, the above equation can be rearranged as
$$
\left(\frac{dN}{dt}\right)_{\rm system} = \lim_{\Delta t \rightarrow 0}\frac{N_{CV, t_0+\Delta t} - N_{CV, t_0}}{\Delta t} + \lim_{\Delta t \rightarrow 0}\frac{N_{{\rm III}, t_0+\Delta t} - N_{{\rm I}, t_0}}{\Delta t}
$$
* Recall from the general relationships between intensive and extensive variables that an extensive variable is found by integrating the intensive variable over the mass of a system, or integrating over the volume of the system. The first term in the above equation can be rewritten as
$$
\lim_{\Delta t \rightarrow 0}\frac{\left[\int_{CV}\rho\eta dV\right]_{t_0+\Delta t} - \left[\int_{CV}\rho\eta dV\right]_{t_0}}{\Delta t} = \frac{\partial}{\partial t}\int_{CV}(\rho\eta) dV
$$
which describes the explicit time variation of the extensive property of the system.
* Consider now the second term on the right-hand side of Eq.X:
$$
\lim_{\Delta t \rightarrow 0}\frac{N_{{\rm III},t_0+\Delta t}-N_{{\rm I},t_0+\Delta t}}{\Delta t} = \lim_{\Delta t \rightarrow 0}\frac{\left[\int_{\rm III}\rho\eta dV\right]_{t_0+\Delta t} - \left[\int_{\rm I}\rho\eta dV\right]_{t_0}}{\Delta t}
$$
* In the limit of vanishingly small time interval, this term represents the rate at which the extensive property $N$ is transported convectively with the fluid motion across the control surface *out of* the control volume.
* Given that the fluid flow can be described by a velocity vector field ${\bf V}$, the convective transport flux across the area $A$ of the control surface can be written
$$
\lim_{\Delta t \rightarrow 0}\frac{N_{{\rm III},t_0+\Delta t} - N_{{\rm I},t_0+\Delta t}}{\Delta t} = \int_{CS}\rho\eta{\bf V}\cdot{\bf n}dA
$$
* The expression ${\bf V}\cdot{\bf n}dA$ is the scalar product (dot product) between the velocity vector and the outward-pointing normal unit vector that describes the control surface. Since ${\bf n}$ is defined as an *outward-normal* unit vector, a positive value of $\int_{CS}\rho\eta{\bf V}\cdot{\bf n}dA$ indicates that $N$ leaves the control volume, while $N$ remains in the *system*.
* Combining Eqs. 3.14, 3.15, and 3.17 yields the **Reynolds transport theorem**, which relates the time rate of change (net accumulation) of and extensive property in a flowing system to a fixed control volume that coincides with the system at an instant in time.
$$
\left(\frac{dN}{dt}\right)_{\rm system} = \int_{CV}\frac{\partial}{\partial t}(\rho\eta)dV + \int_{CS}\rho\eta{\bf V}\cdot{\bf n}dA
$$
* The left-hand side refers to the system, and the right-hand side refers to the control volume that is initially coincident with the system. The right-hand side has two terms. The volume integral term is concerned with the local time rate of change of intensive property within the control volume. This is the accumulation rate of $N$ within the control volume. The surface-integral term is concerned with the net rate at which $N$ is carried *out of the control volume* by convection with the fluid velocity ${\bf V}$ through the control surface.
* It should be recognized that Eq.X is not in itself a conservation equation, but is the relationship between the system and the control volume. It may be instructive, however, to anticipate how the relationship might be used to form a conservation equation. If $N$ represents mass, then $dN/dt = 0$, since a system, by definition, contains a fixed amount of mass. In the case where $N$ represents mass, the corresponding intensive variable is $\eta = 1$. Thus, Eq.X reduces to
$$
\int_{CV}\frac{\partial \rho}{\partial t}dV = -\int_{CS}\rho{\bf V}\cdot{\bf n}dA
$$
* Physically, this equation states that the rate of accumulation of mass (represented by density $\rho$) within the control volume is equal to the net amount of mass that flows across the control surfaces. The leading negative sign of the right-hand side accounts for the fact that a positive value of the control-surface integral indicates flow out of the control volume.
* It is possible, and very useful, to write the surface integral in terms of a volume integral via the Gauss divergence theorem, which states that
$$
\int_{CS}a{\bf G}\cdot{\bf n}dA = \int_{CV}(\nabla\cdot a{\bf G})dV = \int_{CV}(a\nabla\cdot{\bf G})dV
$$
where $a$ is a scaler, ${\bf G}$ is a vector, and ${\bf n}$ is the outward-pointing unit vector at the control surface. The *divergence* of the vector ${\bf G}$ is represented as $\nabla\cdot{\bf G}$, which produces a scalar.
* Using the Gauss divergence theorem, the Reynolds transport theorem (Eq.3.18) can be rewritten as
$$
\begin{align*}
\left(\frac{dN}{dt}\right)_{\rm system}
&= \int_{CV}\frac{\partial}{\partial t}(\rho\eta)dV + \int_{CS}\rho\eta{\bf V}\cdot{\bf n}dA \\
&= \int_{CV}\left(\frac{\partial (\rho\eta)}{\partial t} + \nabla\cdot\rho\eta{\bf V}\right)
\end{align*}
$$
* If the control volume is a vanishingly small one, meaning a differential control volume, then the integrand in Eq.3.21 can be viewed as being constant within the volume. Hence, carrying out the integral is rather simple, yielding
$$
\left(\frac{dN}{dt}\right)_{\rm system} = \left(\frac{\partial(\rho\eta)}{\partial t} + \nabla\cdot\rho\eta{\bf V}\right)\delta V
$$
where $\delta V$ is the volume of the differential control volume. For example, a Cartesian differential control volume as a volume $\delta V = dxdydz$. In a series of manipulations that follow, the terms withing the parentheses are defined as a differential operator called the **substantial derivative**.

#### Substantial derivative
* Although it will be derived in the later chapter, let us assume the mass-conservation equation
$$
\frac{\partial \rho}{\partial t} + \nabla\cdot(\rho{\bf V}) = 0
$$
holds. When the right-hand side of Eq.22 is expanded as
$$
\left(\frac{dN}{dt}\right)_{\rm system} = \left[\rho\frac{\partial \eta}{\partial t} + \eta\frac{\partial \rho}{\partial t} + \rho{\bf V}\cdot\nabla \eta + \eta\nabla\cdot\rho{\bf V}\right]\delta V
$$
it is apparent that the mass-continuity equation eliminates two terms exactly.
* The resulting expression represents the substantial derivative, which is defined as
$$
\frac{D\eta}{Dt} \equiv \frac{\partial \eta}{\partial t} + {\bf V}\cdot\nabla\eta
$$
* The fundamental relationship between a flowing system and an Eulerian control volume, which are coincident at an instant in time, is stated as
$$
\left(\frac{dN}{dt}\right)_{\rm system} = \left[\rho\frac{D\eta}{Dt}\right]_{CV}\delta V
$$
* This equation provides the relationship between the rate of the change of an extensive property $N$ for a system and the substantial derivative of the associated intensive variable $\eta$ in an Eulerian control volume $\delta V$ that is fixed in space.

### Conservation equations
* Overall the objective of this chapter is to cast the conservation equations in the form of partial differential equations (PDEs) in an Eulerian framework with the independent variable being time and the spatial coordinates. The approach combines the notions of conservation laws for systems with the behavior of control volumes fixed in space, through which fluid flows.

#### Mass continuity
* Regardless of what other conservation equations may be appropriate, a bulk-fluid mass-conservation equation is invariably required in representing any fluid-flow situation.
* Under most circumstances, mass cannot be created or destroyed within a control volume. Chemical reaction, for example, may produce or consume individual species, but overall no mass is created or destroyed. Furthermore, the only way that net mass can be transported across the control surface is by convection. While individual species may diffuse across the control surfaces by molecular actions, there can be no net mass transport by molecular diffusion.
* It may be noted that there may be circumstances in two-phase flow where for a given phase there are source or sink terms in describing mass conservation. The conservation of liquid to vapor (and vice versa) causes source terms, namely the creation (or destruction) of mass within a phase. Of course, there must be overall mass conservation between the phases, with not net creation of destruction of mass.
* The mass-continuity equation, actually which was already used in the derivation of the substantial-derivative, is
$$
\left(\frac{dm}{dt}\right)_{\rm system} 
  = \int_{CV}\left(\frac{\partial\rho}{\partial t} + \nabla\cdot\rho{\bf V}\right)dV = 0
$$
* Assuming a vanishingly small differential control volume, the integrand can be assumed to be uniform over the volume. Then the integrand in the above equation is
$$
\begin{align*}
\int_{CV}\left(\frac{\partial\rho}{\partial t} + \nabla\cdot\rho{\bf V}\right)dV
&=\left(\frac{\partial\rho}{\partial t} + \nabla\cdot{\bf V}\right) \int_{CV}dV \\
&=\left(\frac{\partial\rho}{\partial t} + \nabla\cdot{\bf V}\right)\delta V = 0
\end{align*}
$$
* Consequently, the continuity equation can be written in differential-equation form as
$$
\frac{\partial\rho}{\partial t} + \nabla\cdot\rho{\bf V} = 0
$$

#### Navier-Stokes equations
* The Navier-Stokes equations express the conservation of momentum. Together with the continuity equation, which expresses the conservation of mass, these equations are th fundamental underpinning of fluid mechanism.
* They are nonlinear partial differential equations that in general cannot be solved by analytical means. Although the Navier-Stokes equations defy general analytic solution, solution by numerical methods has become commonplace as a practical engineering analysis and design tool.
* The principle of momentum conservation for a system can be written generally as
$$
\frac{d{\bf P}}{dt} = \sum{\bf F}
$$
where the rate of change of momentum ${\bf P}$ is caused by the forces ${\bf F}$ acting on the system. For a solid body, momentum conservation is usually written as the familiar equation (Newton's second law)
$$
{\bf F} = m{\bf a}
$$
where ${\bf F}$ is the vector of forces acting on the mass $m$ and ${\bf a}$ is the mass' acceleration vector. The acceleration of a solid body can be simply expressed as
$$
{\bf a} = \frac{d {\bf V}}{dt}
$$
where ${\bf V}$ is the solid-body's velocity vector.
* For a fluid flow, the RTT establishes the relation between an system (where the momentum balance applies directly) and an Eulerian control volume (through which fluid flows).
* Here, the extensive variable $N$ is the momentum vector ${\bf P} = m{\bf V}$ and the intensive variable $\eta$ is the velocity vector ${\bf V}$. Thus, application of the RTT yields the following vector equation:
$$
\left[\rho\frac{D{\bf V}}{Dt}\right]\delta V = \frac{d{\bf P}}{dt}
= \sum{\bf F} = \sum{\bf F}_{\rm body} + \sum{\bf F}_{\rm surface}
$$
where $D{\bf V}/Dt$ is the substantial derivative and $\delta V$ is the volume of a differential control volume.
* It is useful to think of two different kinds of forces, one that acts over the volume of a fluid element and the other that acts on the elements' surface. The most common body force is exerted by the effect of gravity.
$$
{\bf F}_{\rm body} = m{\bf g}
$$
* For a differential control volume, it is convenient to divide Eq.X by the volume of the control volume $\delta V$. This leads to a general differential-equation statement of the Navier-Stokes equations as
$$
\rho\frac{D{\bf V}}{Dt} 
= {\bf f}_{\rm body} + {\bf f}_{\rm surface}
= \rho{\bf g} + \nabla\cdot{\textsf T} = \rho{\bf g} - \nabla p + \nabla\cdot{\textsf T}'
$$
where ${\bf f}$ is the force per unit volume ($N m^{-3}$), which is also the momentum-flux vector. The surface forces are determined from the stress tensor. The pressure-gradient term can be separated from the general stress term by using the deviatoric stress tensor.
* The gravitational body-force term becomes
$$
{\bf f}_{\rm body} = \rho{\bf g}
$$
* In the case of a constant-viscosity fluid, the Navier-Stokes equations can be significantly simplified as following general vector form,
$$
\begin{align*}
\rho\frac{D{\bf V}}{Dt}
&= \rho\left[\frac{\partial {\bf V}}{\partial t} + ({\bf V}\cdot\nabla){\bf V}\right] \\
&= \rho\left[\frac{\partial {\bf V}}{\partial t} + \nabla\left(\frac{{\bf V}\cdot{\bf V}}{2}\right)-{\bf V}\times(\nabla\times{\bf V})\right] \\
&= {\bf f} - \nabla p - \mu\nabla\times(\nabla\times{\bf V}) + (\kappa+2\mu)\nabla(\nabla\cdot{\bf V}) \\
&= {\bf f} - \nabla p + \mu\nabla^2{\bf V} + (\kappa+\mu)\nabla(\nabla\cdot{\bf V})
\end{align*}
$$
A vector identity defined the Laplacian of a vector as
$$
\begin{align*}
\nabla^2{\bf V} &\equiv \nabla(\nabla\cdot{\bf V}) - \nabla\times(\nabla\times{\bf V}) \\
                &= \nabla(\nabla\cdot{\bf V}) - \nabla\times\omega
\end{align*}
$$
These general vector operations can be expanded into any particular coordinate system of interest. In this form, the vorticity $\omega \equiv \nabla\times{\bf V}$ appears explicitly. Also recall that for incompressible flows $\nabla\cdot{\bf V}$ vanishes, thus eliminating the final term.

#### Species diffusion
* To derive the species-continuity equations that follow, it is important to establish some relationships between mass fluxes and species concentration fields.

##### Mass and mole measure
* Working with multicomponent mixtures requires quantifying the amounts of various chemical constituents that comprise the mixture. In the conservation equations the mass fraction is the most appropriate measure, since mass is a conserved quantity. By definition, the mass fraction is
$$
Y_k \equiv \frac{\rho_k}{\rho}
$$
where $\rho_k$ is the mass density of the k-th species and $\rho$ is the total density. Clearly,
$$
\sum_{k=1}^K\rho_k = \rho, \ \ \sum_{k=1}^K Y_k = 1
$$
* Chemical behaviors, such as chemical reactions, are usually best quantified on a molar basis. That is, a certain number of moles of one species reacts with a certain number of moles of another to produce a certain amount of product species. Here a mole fraction, not a mass fraction, is the most appropriate measure of the mixture composition. The mole fraction $X_k$ is the number of moles of species $k$ in a volume divided by the total moles in the volume. For an ideal gas, the mole fraction is related to mass fraction as
$$
X_k = Y_k\frac{\bar{W}}{W_k}
$$
where $W_k$ is the molecular mass of species $k$, and $\bar{W}$ is the average molecular mass defined by $\bar{W} = \sum_{k=1}^K X_k W_k$.
* The molar concentration $[X_k]$ is the measure of chemical composition that is most natural for the description of chemical reaction. The molar concentration of an ideal gas can be written as
$$
[X_k] = \frac{p}{RT}X_k = [X_k]_{\rm tot}X_k
$$
where $p/RT = [X_k]_{\rm tot}$ is the total concentration of a gas mixture. The concentration is measured as moles per unit volume ($mol m^{-3}$).
* The partial pressure of a species $p_k$ in an ideal-gas mixture is closely related to the mole fraction;
$$
\sum_{k=1}^K\frac{p_k}{p} = \sum_{k=1}^K X_k = 1
$$
Thus the mass fraction of $k$-th species $Y_k$ can be calculated from $p_k$ as
$$
Y_k = \frac{p_k}{p}\frac{W_k}{\bar{W}}
$$

##### Diffusive mass flux
* Whenever there are chemical-composition variations in a fluid, there is a tendency for chemical species to be transported by molecular diffusion from regions of higher concentration to regions of lower concentration.
* In the simplest theory, Fick's law, the diffusive mass flux of a species depends linearly on the negative concentration gradient of the species concentration with a proportionality constant called a diffusion coefficient. The negative sign sets the direction of the flux toward the low-concentration region.
* The diffusive mass-flux vector ($kg m^{-2} s^{-1}$) can be represented as
$$
{\bf j}_k = \rho Y_k {\bf V}_k
$$
where ${\bf V}_k$ is the *diffusion-velocity* vector for the k-the species.
* As represented by **Fick's law**,
$$
{\bf V}_k = -\frac{1}{X_k}D_{km}' \nabla X_k
$$
where $D_{km}'$ represents a "mixture-averaged" diffusion coefficient for species $k$ relative to the rest of the multicomponent mixture. The mixture-averaged diffusion coefficient can be evaluated in terms of the *binary diffusion coefficient ${\mathcal D}_{jk}$ as
$$
D_{km}' = \frac{1-Y_k}{\sum_{j \ne k}X_j/D_{jk}}
$$
* The species mass-flux vectors can be written in terms of the mole-fraction gradient as
$$
{\bf j}_k = -\rho\frac{Y_k}{X_k}D_{km}'\nabla X_k = –\rho\frac{W_k}{\bar{W}}D_{km}'\nabla X_k
$$
* It is important to realize that there cannot be a net transport of mass by diffusive action within a homogeneous multicomponent field. The transport of some species in one direction must be balanced by transport of other species in the other direction.
$$
\sum_{k=1}^K {\bf j}_k\cdot{\bf n}dA = \sum_{k=1}^K\rho Y_k{\bf V}_k\cdot{\bf n}dA = 0
$$
Here ${\bf n}dA$ is some differential area, with ints spatial orientation specified by an outward-normal unit vector ${\bf n}$. Since the equation is true for any differential area, it is generally true that
$$
\sum_{k=1}^K {\bf j}_k = 0
$$
* The discussion to this point in the section has considered only diffusive mass transport. It should be noted that the net mass transport of a species $k$ crossing a certain area $dA$ is the sum of diffusive and convective contributions. This is stated as
$$
\dot{m}_k = \rho Y_k({\bf V}+{\bf V}_k)\cdot{\bf n}dA
$$
where ${\bf V}$ is the bulk fluid velocity.

##### Stefan-Maxwell equations
* The evaluation of the diffusion velocities requires the evaluation of the multi-component diffusion coefficients from the binary diffusion coefficients.
* In **Stefan-Maxwell equations**, the diffusion velocities are related *implicitly* to the temperature, pressure, and mole-fraction gradients as
$$
\begin{align*}
\nabla X_k = \sum_{j=1}^K\frac{X_K X_J}{{\mathcal D}_{kj}} ({\bf V}_j - {\bf V}_k) + (Y_k-X_k)\frac{\nabla p}{p} \\
                       + \sum_{j=1}^K\frac{X_k X_j}{\rho{\mathcal D}_{kj}} \left(\frac{D_j^T}{Y_j} - \frac{D_k^T}{Y_k} \right) \frac{\nabla T}{T}
\end{align*}
$$
Note that the Stefan-Maxwell equations involve the *binary diffusion coefficients* ${\mathcal D}_{kj}$, and not hte ordinary multi-component diffusion coefficients $D_{jk}$.
* In the mixture-averaged formalism, the ordinary multicomponent diffusion coefficient must be evaluated from the binary diffusion coefficients. This requires the matrix inversion, which can add computational complexity. The Stefan-Maxwell formulation uses the binary diffusion coefficient directly, but requires the solution of the linear system to evaluate the diffusion velocities. At each point in a flow field, one could solve the system of equations (Eq.X, above) to determine the diffusion-velocity vector.

#### Species conservation
* The continuity equation is a statement of overall mas conservation, and no distinction is made as to the chemical identity of individual species in the flow.
* When considering the mass continuity of an individual species in a multicomponent mixture, there can be, and typically is, diffusive transport across the control surfaces and the production or destruction of an individual species by volumetric chemical reaction.

##### Conservation law for individual species
* In addition to overall mass conservation, the conservation law for individual chemical species is necessary in the simulation of multicomponent and chemically reacting flow.
* Begin with the system law, stated as
$$
\left(\frac{d m_k}{dt}\right)_{\rm system}
  = -\int_{CS}{\bf j}_k\cdot{\bf n}dA + \int_{CV}\dot{\omega}_k W_k dV
$$
In other words, the rate of change of the mass of species $k$ in the system is balanced by the rate at which the $k$-th species diffuses across the control surfaces and the rate at which the $k$-th species is created or consumed by chemical reaction.
* The species mass $m_k$ is the extensive variable and the associated intensive variable is the mass fraction $Y_k$. Based on the RTT, the convective flus of the $k$-th species is accommodated by the substantial derivative. Thus, exchanging the system view for the control-volume view,
$$
\left[\rho\frac{DY_k}{Dt}\right]\delta V = -\int_{CS}{\bf j}_k\cdot{\bf n}dA + \int_{CV}\dot{\omega}_k W_k dV
$$
* The first term on the right-hand side describes the net species mass flux that diffuses into the system. The lead minus sign is required to accommodate the fact that ${\bf n}$ is defined as an *outward-normal* unit vector. The second term represents the production/consumption of the $k$-th species by chemical reaction.
* The surface integral can be converted to a volume integral using the Gauss divergence theorem, yielding
$$
\left[\rho\frac{DY_k}{Dt}\right]\delta V = \int_{CV}(-\nabla\cdot{\bf j}_k + \dot{\omega}_k W_k)dV
$$
* For a vanishingly small differential control volume, the integrand can be considered constant. Thus, the integration yields
$$
\int_{CV}(-\nabla\cdot{\bf j}_k+\dot{\omega}_k W_k)dV = (-\nabla\cdot{\bf j}_k+\dot{\omega}_k W_k)\delta V
$$
where $\delta V$ is the volume of the differential control volume. After dividing each term by $\delta V$, a partial differential equation (the **species-continuity equation**) emerges as
$$
\rho\frac{DY_k}{Dt} = -\nabla\cdot{\bf j}_k + \dot{\omega}_k W_k
$$
* The variable $\dot{\omega}_k$ is used to denote the volumetric molar production/destruction rate of species $k$ by chemical reaction ($mol m^{-3} s^{-1}$). A great many reactions may participate in the production of $k$-th species that comprise the multicomponent mixture. While the molar production/destruction rate is natural from the chemical reaction formalism, the mass rate of change is more natural choice in the species mass balance. In this case, $\dot{\omega_k}W_k$ represents the mass rate of change ($kg m^{-3} s^{-1}$), where $W_k$ is the molar mass of species $k$. Chemical reaction converts some species to other species; hence the mass represented by individual species changes via chemical reaction. However, a homogeneous chemical reaction cannot create or destroy net mass. Therefore, in a homogeneous mixture
$$
\sum_{k=1}^K \dot{\omega}_k W_k = 0
$$

##### Summation of species continuity
* By expanding the substantial derivative of Eq.X, the species-continuity equation can be written as 
$$
\int_{CV}\left(\frac{\partial(\rho Y_k)}{\partial t} + \nabla\cdot\rho Y_k {\bf V}\right)dV = \int_{CV}(-\nabla\cdot{\bf j}_k + \dot{\omega}_k W_k)dV
$$
Integrating over the differential control volume and dividing by the differential volume yields the species-continuity equation in differential-equation form as
$$
\frac{\partial(\rho Y_k)}{\partial t} + \nabla\cdot\rho Y_k {\bf V} = -\nabla\cdot{\bf j}_k + \dot{\omega}_k W_k
$$
Summing over all species yields
$$
\begin{align*}
\sum_{k=1}^K \left(\frac{\partial(\rho Y_k)}{\partial t} + \nabla\cdot\rho Y_k {\bf V} \right)
&= \sum_{k=1}^K(-\nabla\cdot{\bf j}_k + \dot{\omega}_k W_k) \\
\frac{\partial\rho \sum_{k=1}^K Y_k}{\partial t} + \nabla\cdot\rho\sum_{k=1}^K Y_k {\bf V} &= -\nabla\cdot\sum_{k=1}^K{\bf j}_k + \sum_{k=1}^K\dot{\omega}_k W_k
\end{align*}
$$
As discussed earlier (Eqs.X and X), both therms on the right-hand side are zero. Also, since $\sum_{k=1}^K Y_k = 1$, the overall mass-continuity equation is recovered.
$$
\frac{\partial \rho}{\partial t} + \nabla\cdot\rho{\bf V} = 0
$$
It may also be noted that summing the system representation must also produce the starting point for derivation of the overall mass-continuity equation,
$$
\sum_{k=1}^K\frac{dm_k}{dt} = \frac{d}{dt}\sum_{k=1}^K m_k = \frac{dm}{dt} = 0
$$
In other words, by definition, the net mass in a system cannot change.

#### Conservation of energy
* Derivation of the energy equation begins with the first law of thermodynamics, which includes both thermal and mechanical energy. In fluid mechanics, however, the primary interest is in the thermal energy, which is represented by enthalpy or temperature fields. By subtracting the mechanical-energy components from the total energy equation, it is possible to derive a *thermal-energy equation* that serves as the basis for the most subsequent analysis.
* One important purpose of the energy equation is to describe and predict the fluid temperature fields. The energy equation must be closely coupled to the continuity and Navier-Stokes equations, which describe the velocity fields.
* The coupling comes through the convective terms in the substantial derivative, which, of course, involve the velocities.
* The Navier-Stokes equations are also coupled to the energy equation, since the density and other properties usually depend on temperature.
* Chemical reaction and molecular transport of chemical species can also have a major influence on the thermal energy of a flow.
* The energy equation is a statement of the first law of thermodynamics, just as the Navier-Stokes equations are statements of the Netwon's second law, ${\bf F} = m{\bf a}$.
* For a $system$, the first law states that the rate of the total energy change equals the rate of heat transferred to the system plus the rate of work done on the system. That is,
$$
\frac{dE_t}{dt} = \frac{dQ}{dt} + \frac{dW}{dt}
$$
where $E_t$ is the extensive variable that represents the total energy stored in a system ($J$). $Q$ represents heat added *to* the system, and $W$ represents work done *on* the system. The total energy includes internal, kinetic, and potential energy
$$
\frac{E_t}{m} = e_t = \left(e + \frac{1}{2}({\bf V}\cdot{\bf V}) - {\bf g}\cdot{\bf r}\right)
$$
* The intensive variable is the total specific energy $e_t$, where $m$ is the mass of the system ($kg$). There are three contributions to the total energy; the specific internal energy $e$ has contributions that represent the random motion of molecules associated with non-zero temperature. It also has contributions that represent the potential energy associated with chemical bonds. The second term represents directed kinetic energy of the fluid. The third term represents the potential energy associated with the downward-directed acceleration of gravity ${\bf g}$, and ${\bf r}$ is the displacement of a fluid packet relative to some reference.
* Recalling the general relationship between a system and a control volume, the left-hand side of the energy balance can be written as
$$
\left(\frac{dE_t}{dt}\right)_{\rm system} = \left[\rho\frac{\partial e_t}{\partial t} + \rho{\bf V}\cdot\nabla e_t \right]_{CV}\delta V = \left[\rho\frac{De_t}{Dt}\right]_{CV}\delta V
$$
* This equation represents the rate of change of the system's total energy in terms of hte substantial derivative for a flowing system applied to an Eulerian control volume fixed in space.
* Differentiating the definition of total energy yields an expression for the substantial derivative of the total energy as
$$
\frac{De_t}{Dt} = \left(\frac{De}{Dt} + {\bf V}\cdot\frac{D{\bf V}}{Dt} - {\bf g}\cdot{\bf V}\right)
$$
* The time derivative of the displacement vector ${\bf r}$ is the velocity ${\bf V}$, which assumes that the fluid system is moving with the fluid velocity.
* The left-hand side of the energy equation now represents the convective transport, and it remains to develop the heat-transfer and work terms on the right-hand side
$$
\rho\left(\frac{De}{Dt}+{\bf V}\cdot\frac{D{\bf V}}{Dt}-{\bf V}\cdot{\bf V}\right)\delta V = \frac{dQ}{dt}+\frac{dW}{dt}
$$

##### Heat-transfer rate
* The next task is to develop expressions for the heat-transfer and work terms in Eq.X. Consider two contributions to the heat transfer crosses the surfaces of a control volume.
* The first is thermal conduction via the Fourier's law, which behaves in the same way for a fluid as it does for a solid. The second contribution is associated with energy that crosses the control surfaces as chemical species diffuse into and out of the control volume.
* Chemical reaction is often thought to provide a "source of heat" or "heat release" in the control volume.
* Importantly, this contribution is *not* and internal heat source but rather represents a change in the internal energy by the breaking and forming of chemical bonds.
* Although the temperature may change due to chemical reactions, the total energy of the system does not change as a result.
* The Fourier's law states that *heat flux* ($W m^{-2} or J s^{-1}m^{-2}$) is proportional to the negative temperature gradient, with the constant of proportionality being the thermal conductivity $\lambda$,
$$
{\bf q} = -\lambda \nabla T
$$
* Clearly, the heat flux is a vector whose direction is defined by the temperature gradient. The net heat $dQ/dt$ ($W or J s^{-1}$) that crosses a control surface into the volume by thermal conduction is given by
$$
\left(\frac{dQ}{dt}\right)_{\rm conduction} = -\int_{CS}{\bf q}\cdot{\bf n}dA = \int_{CS}\lambda\nabla T\cdot{\bf n}dA
$$
where ${\bf n}$ is an outward-normal-pointing unit vector. When heat flows into the control volume, $dQ/dt$ is positive. Thus, the negative sign is needed because, when ${\bf q}$ flows into the control volume, its direction is opposite to ${\bf n}$, which points outward. Consequently, with the minus sign in front of the integral, $dQ/dt$ is positive when ${\bf q}$ is the opposite direction of ${\bf n}$.

* Turn now to the heat transfer associated with the species mass fluxes that diffuse across the control surfaces, which is stated as
$$
\left(\frac{dQ}{dt}\right)_{\rm species} = -\sum_{k=1}^K\int_{CS}h_k{\bf j}_k\cdot{\bf n}dA
$$
where $h_k$ is the enthalpy of species $k$. Each species carries with its energy as it diffuses across the control surface. The minus sign is needed, as it was for the thermal conduction, because when ${\bf n}$ and ${\bf j}_k$ have opposite directions, energy enters the control volume resulting in positive $dQ/dt$. The fact that enthalpy $h_k$, rather than internal energy $e_k$, represents the energy content is because "flow work" also contributes to the energy exchange. Using the definition of enthalpy $h_k = e_k + p_k / \rho_k$, and substituting for the mass flux ${\bf j}_k = \rho_k{\bf V}_k = \rho Y_k{\bf V}_k$, the above equation can be rewritten as
$$
\left(\frac{dQ}{dt}\right)_{\rm species} = -\sum_{k=1}^K\int_{CS}\left(e_k+\frac{p_k}{\rho_k}\right)\rho_k{\bf V}_k\cdot{\bf n}dA
$$
and then
$$
\left(\frac{dQ}{dt}\right)_{\rm species} = -\sum_{k=1}^K\left(\int_{CS}e_k\rho_k{\bf V}_k\cdot{\bf n}dA + \int_{CS}p_k{\bf V}_k\cdot{\bf n}dA \right)
$$
* The first term on the right-and side represents the internal energy that is carried across the control surface with the diffusion velocity. The second term represents the "$pV$" work caused by the force exerted at the control surface by the pressure as it acts on the fluid that is moving with the diffusion velocity. While this type of flow work could be grouped with the $dW/dt$ term, it is a long-standing convention to group it with $dQ/dt$ using the enthalpy as the energy measure.
* The net heat-transfer rate is the sum of two contributions,
$$
\begin{align*}
\frac{dQ}{dt} &= \left(\frac{dQ}{dt}\right)_{\rm conduction} + \left(\frac{dQ}{dt}\right)_{\rm species} \\
              &= \int_{CS}\lambda\nabla T\cdot{\bf n}dA - \sum_{k=1}^K\int_{CS}h_k{\bf j}_k\cdot{\bf n}dA
\end{align*}
$$
* Using the Gauss divergence theorem, the surface integrals can be rewritten as a volume integrals, yielding
$$
\frac{dQ}{dt} = \int_{CV}(\nabla\cdot\lambda\nabla T)dV - \sum_{k=1}^K\int_{CV}\nabla\cdot h_k{\bf j}_k dV
$$
* Assuming a vanishingly small control volume, so that there is no variation of the integrands within the control volume, the integrals are easily accomplished as 
$$
\frac{dQ}{dt} = \left( \nabla\cdot\lambda\nabla T - \sum_{k=1}^K\nabla\cdot h_k{\bf j}_k \right)\delta V
$$
where $\delta V$ is the volume of a differential control volume.

#### Rate of work
* Turn now to the work term $dW/dt$. The stress tensor causes forces on the surfaces of a control volume, through which fluid is moving, with the result being work.
* On any arbitrary surface $dA$, the resultant stress can be represented as a vector $\tau$.
* The velocity at the surface is represented as a vector ${\bf V}$.
* At any point in the flow field, the stress state is represented by a second-order tensor ${\textsf T}$.
* On a surface, which may represent some portion of the control surface that bounds a control volume, the stress is represented as a vector.
* The relationship between the stress tensor ${\textsf T}$ at a point and the stress vector $\tau$ on a particular surface that passes through the point is given as (Eq.3.113)
$$
\tau = {\bf n}\cdot{\textsf T}
$$
where ${\bf n}$ is the outward-directed unit vector of the surface $dA$.

* In general, the rate of work done at some surface moving with velocity ${\bf V}$ is
$$
$$
* Here, $\tau$ is the stress *vector*. Note that both normal stress and shear stress contribute ot work. That is, work is associated with both dilation and deformation.
* Now let us convert the surface integral to a volume integral, using the Gauss theorem. However, the work-rate integral does not appear to be in a form directly suitable for the divergence theorem, since it does not involve the scalar product of a vector with the normal component of the area. Therefore some operations are needed.
* First, recognize the associative property of the scalar product of the two vectors:
$$
\tau\cdot{\bf V} = {\bf V}\cdot\tau
$$
* Then, an identity permits the following less-than-obvious step:
$$
{\bf V}\cdot({\bf n}\cdot{\textsf T}) = {\bf n}\cdot({\bf V}\cdot{\textsf T}^T)
$$
* Here, ${\textsf T}^T$ is the transpose of the stress tensor. Since the stress tensor is symmetric, the transpose and the original tensor are identical. Thus
$$
\frac{dW}{dt} = \int_{CS}{\bf n}\cdot({\bf V}\cdot{\textsf T})dA = \int_{CS}({\bf V}\cdot{\textsf T})\cdot{\bf n}dA
$$
* Now we have ${\bf n}dA$ term thus the Gauss theorem can be used, as
$$
\int_{CS}({\bf V}\cdot{\textsf T})\cdot{\bf n} dA = \int_{CV}\nabla\cdot({\bf V}\cdot{\textsf T})dV
$$
* The rate of work done on a differential control volume $\delta V$ by the stress and velocity fields is expressed as
$$
\frac{dW}{dt} = \nabla\cdot({\bf V}\cdot{\textsf T})\delta V
$$

#### Total energy equation in vector form
* With all the individual terms in hand, the full energy equation can be assembled and represented in compact vector form as
$$
\begin{align*}
\rho\frac{De_t}{Dt} =& \frac{dQ}{dt} + \frac{dW}{dt} \\
\rho\left(\frac{De}{Dt}+{\bf V}\cdot\frac{D{\bf V}}{Dt}-{\bf g}\cdot{\bf V}\right) =& \nabla\cdot(\lambda\nabla T) \\
&- \sum_{k=1}^K\nabla\cdot h_k{\bf j}_k + \nabla\cdot({\bf V}\cdot{\textsf T})
\end{align*}
$$
* In fluid mechanics, it is unusual to formulate problems in the context of the total energy equation; rather, the *thermal-energy equation* is more practically useful. The thermal-energy equation is formed by subtracting the *mechanical-energy equation* from the total energy equation.

#### Mechanical energy
* The mechanical-energy equation is formed by the scalar product of the velocity vector and the momentum equations (Navier-Stokes equations).
* By a vector-tensor identity for symmetric tensors, the work-rate term $dW/dt$ can be expanded as
$$
\nabla\cdot({\bf V}\cdot{\textsf T}) = {\bf V}\cdot(\nabla\cdot{\textsf T}) + {\textsf T}:\nabla{\bf V}
$$
The dyadic product ($:$) of the stress tensor and the velocity-gradient tensor produces a scalar. Note that the work is a scalar quantity.
* The first term on the right-hand side of the above equation includes the divergence of the stress tensor, which also appears in the vector form of the momentum equations (Navier-Stokes equations). The momentum equation can be easily rearranged as
$$
\nabla\cdot{\textsf T} = \rho\left(\frac{D{\bf V}}{Dt} - {\bf g}\right)
$$
It then follows that the first term in Eq.X can be written as
$$
{\bf V}\cdot(\nabla\cdot{\textsf T}) = \rho\left({\bf V}\cdot\frac{D{\bf V}}{Dt} - {\bf g}\cdot{\bf V}\right)
$$
This equation describes the conservation of *mechanical energy*.
* It is apparent that all the terms in the above also appear directly on the left-hand side of the total energy equation (Eq.4.148). Therefore, subtraction is done by removing the mechanical-energy contribution from the total energy equation.

#### Thermal energy
* In the general vector form, the thermal-energy equation may be stated as
$$
\rho\frac{De}{Dt} = - \nabla\cdot{\bf q} - \sum_{k=1}^K\nabla\cdot h_k{\bf j}_k + {\textsf T}:\nabla{\bf V}
$$
where ${\bf q}$ is the heat flux. The objective of the following series of manipulations is to replace the internal energy $e$ on the left-hand side with the enthalpy $h$, which provides a form of the thermal-energy equation that is usually more convenient.
* By using the deviatroic stress tensor (Eq.3.141), the thermodynamic pressure can be separated from the ${\textsf T}:\nabla{\bf V}$ term as
$$
{\textsf T}:\nabla{\bf V} = {\textsf T'}:\nabla{\bf V} - p\nabla\cdot{\bf V}
$$
* To show that the $p\nabla\cdot{\bf V}$ term emerges, the "pressure tensor" may be written as
$$
{\textsf p} = 
\begin{pmatrix}
p &   &   \\
  & p &   \\
  &   & p
\end{pmatrix}
 = p{\textsf I}
$$
where ${\textsf I}$ is the identity matrix. Then, using the identity stated in Eq.4.149,
$$
{\textsf p}:\nabla{\bf V} = p\nabla\cdot({\bf V}\cdot{\textsf I})-p{\bf V}\cdot({\nabla\cdot{\textsf I}}) = p\nabla\cdot{\bf V}
$$
In this expression, note that ${\bf V}\cdot{\textsf I} = {\bf V}$ and $\nabla\cdot{\textsf I} = 0$.
* The overall mass-continuity equation (Eq. 4.10),
$$
\frac{1}{\rho}\frac{D\rho}{Dt} = -\nabla\cdot{\bf V}
$$
leads to an alternative way to express the $p\nabla\cdot{\bf V}$ term tat appears in Eq. 4.153. That is,
$$
p\nabla\cdot{\bf V} = -\frac{p}{\rho}\frac{D\rho}{Dt} = \rho\frac{D}{Dt}\left(\frac{p}{\rho}\right) - \frac{Dp}{Dt}
$$
With these substitutions the thermal-energy equation becomes
$$
\rho\frac{D}{Dt}\left(e+\frac{p}{\rho}\right) = \rho\frac{Dh}{Dt} = \frac{Dp}{Dt} - \nabla\cdot{\bf q} - \sum_{k=1}^K\nabla\cdot h_k{\bf j}_k + {\textsf T'}:\nabla{\bf V}
$$
where the enthalpy $h = e + p/\rho$ was used.
* $p\nabla\cdot{\bf V}$ term is called the "flow work", and it is the work associated with pressure around the control surfaces. As shown in Eq.X, this term can be expanded in two terms. One is combined with the internal energy $e$ to introduce the enthalpy $h$. The other term, $Dp/Dt$, is fortunately be neglected for most low-speed flows, as
$$
\frac{Dp}{Dt} = \frac{\partial p}{\partial t} + {\bf V}\cdot(\nabla p)
$$
and the pressure derivatives and velocities are often small there. However, since the magnitude of the pressure $p$ can often be large (recall that atmospheric pressure is $10^5 Pa$), $p\nabla\cdot{\bf V}$ needs to be retained.

#### Viscous dissipation
* The thermal-energy equation now has a single term that involves the viscosity; it is called the dissipation function or viscous dissipation
$$
\Phi = {\textsf T'}:\nabla{\bf V}
$$
* This term represents the irreversible conversion of kinetic energy into thermal energy, and must be always be positive. In other words, irreversible work must increase thermal energy in the flow.
* For low-speed flow of gases, viscous dissipation is rarely important. However, in high-speed flows such as supersonic flows, viscous dissipation is important. Also for the flow of high-viscosity fluids like oils in a journal bearing, viscous dissipation must be considered.

#### Thermal energy equation
* The thermal-energy equation is commonly written in the form
$$
\rho\frac{Dh}{Dt} = \frac{Dp}{Dt} + \nabla\cdot(\lambda\nabla T) - \sum_{k=1}^K\nabla\cdot h_k{\bf j}_k + \Phi
$$
* In this form, the Fourier's law is substituted for the heat flux. The thermal conductivity $\lambda$ is the average conductivity of the fluid mixture.
* The thermal-energy equation above has no explicit source term to describe the heat related associated with chemical reaction. Nevertheless, as stated, the thermal-energy equation does fully accommodate chemical reaction. As is described subsequently, the thermal effects of chemical heat release are captured in the enthalpy term on the left-hand side.

##### Ideal gas
* The majority of applications considered in this book is the mixtures of ideal gases. The thermodynamic properties of the mixture that appear in the energy equation are evaluated as mass-weighted averages of the individual species properties. Using the mass fraction $Y_k$ and the species' enthalpy $h_k$, the enthalpy can be written as $h=\sum_{k=1}^K Y_k h_k$ and this plays a central role in the thermal-energy equation. From the definition of $h$, the substantial derivative can be expanded as
$$
\frac{Dh}{Dt} = \sum_{k=1}^K\left(Y_k\frac{Dh_k}{Dt} + h_k\frac{DY_k}{Dt}\right)
$$
* For an ideal gas, where the specific heat is defined in terms of enthalpy as
$$
c_{pk} \equiv \left(\frac{\partial h_k}{\partial T}\right)_p
$$
or $dh_k = c_{pk}dT$. Therefore, the enthalpy derivative can be written in terms of a temperature derivative as
$$
\frac{Dh_k}{Dt} = c_{pk}\frac{DT}{Dt}
$$
* This equation *does not* imply that specific heats are constant. Indeed, they are generally the functions of temperature. Rather, there is not a specific-heat derivative, because of the definition of specific heat (Eq.X).
* Since the mixture specific heat for an ideal gas can be written as a mass-weighted sum of the specific heat of species, $c_p = \sum_{k=1}^K Y_k c_{pk}$
$$
\frac{Dh}{Dt} = c_p\frac{DT}{Dt} + \sum_{k=1}^K h_k \frac{DY_k}{Dt}
$$
* This equation shows that the rate of $Y_k$ change in the mixture contributes directly to the enthalpy change. Recall from the species-continuity equation that there are two contributions to the rate of chemical species change: molecular diffusion across the control surface and homogeneous chemical reaction within the control volume. Substituting the species-continuity equation (Eq.X) yields
$$
\rho c_p \frac{DT}{Dt} + \sum_{k=1}^K h_k(-\nabla\cdot{\bf j}_k + \dot{\omega}_k W_k) \\
= \frac{Dp}{Dt} + \nabla\cdot(\lambda \nabla T) - \sum_{k=1}^K\nabla\cdot h_k{\bf j}_k + \Phi
$$
* On expanding the enthalpy-flux term on the right-hand side as
$$
\sum_{k=1}^K \nabla\cdot h_k{\bf j}_k = \sum_{k=1}^K h_k\nabla\cdot{\bf j}_k + \sum_{k=1}^K{\bf j}_k\cdot\nabla h_k
$$
* The $h_k\nabla\cdot{\bf j}_k$ terms that appear on both sides cancel. Also, since $dh_k = c_{pk}dT$,
$$
\sum_{k=1}^K{\bf j}_k\cdot\nabla h_k = \sum_{k=1}^K c_{pk}{\bf j}_k\cdot\nabla T
$$
* The ideal-gas thermal-energy equation is finally simplified to
$$
\rho c_p \frac{DT}{Dt} = \frac{Dp}{Dt} + \nabla\cdot(\lambda\nabla T) - \sum_{k=1}^K c_{pk}{\bf j}_k\cdot\nabla T - \sum_{k=1}^K h_k\dot{\omega_k}W_k + \Phi
$$
* The physical interpretation of the thermal-energy equation in this form is that the rate of the temperature change is influenced through terms involving mechanical compression, heat conduction, diffusive flux of thermal enthalpy, heat of chemical reaction, and viscous dissipation.
* Note that the "heat source" due to chemical reaction is not really a source term per se. Chemical reaction breaks and forms chemical bonds, causing the temperature to increase or fall. However, the total energy in the system is not altered by the reaction. Rather, potential energy in the form of chemical bonds is converted to thermal energy in the form of temperature change. Thus, the temperature form of the energy equation has a chemical reaction term, whereas the enthalpy (or internal energy) form does not (e.g. Eq.X).
* For low-speed flows, both the mechanical compression $Dp/Dt$ and the viscous dissipation $\Phi$ are very small and can be safely neglected.

#### OLD
#### energy conservation
* The heat transport is described with the following equation.
$$
\frac{\partial(\rho h)}{\partial t} + \nabla\cdot(\rho h u) = \nabla\cdot(\lambda \nabla T) + Q_{GR}
$$
<!---
\rho C_p\frac{\partial T}{\partial t}+\rho C_p {\bf u}\cdot\nabla T = \nabla(\lambda \nabla T)
Here, $C_p$ is the heat capacity at constant pressure, and $\lambda$ is the thermal conductivity. $\rho$ is the density.
-->
Here, h is the relative? enthalpy, and $\lambda$ is the thermal conductivity. $\rho$ is the density. $Q_{GR}$ is the heat of reaction.

#### momentum transport
* The laminar flow in the fluid domain is described by
$$
\begin{split}
&\rho\frac{\partial {\bf u}}{\partial t} +\rho({\bf u}\cdot\nabla){\bf u}
	= \nabla\cdot( -{\bf pI+K+F} ) +\rho g \\
&\frac{\partial \rho}{\partial t} + \nabla\cdot(\rho u) = 0 \\
&{\bf K} = \mu \left(\nabla{\bf u} + \nabla{\bf u}^T \right)
\end{split}
$$
${\bf K}$ is called the viscous stress tensor. ${\bf u}$ is flow velocity vector. ${\bf I}$ is the unit momentum vector, and ${\bf F}$ is the volume force vector. $g$ is the gravity acceleration constant. $\mu$ is the dynamics viscosity of the fluid.
* The first equation is called the **Navier-Stokes equation**, and the second one is called **continuity equation**.
* The Navier-Stokes equation describes the conservation of momentum for an incompressible Newtonian fluid.

* For the momentum transport in the porous materials, the **Brinkman equation** is often used instead of the Navier-Stokes equation.
$$
\frac{\rho}{\varepsilon_P}\left({\bf u}\cdot\nabla\frac{\bf u}{\varepsilon_P}\right)
= \nabla
	\left[-{\bf pI}+\frac{\mu}{\varepsilon_P}(\nabla{\bf u}+\nabla{\bf u}^T)
		-\frac{2\mu}{3\varepsilon_P}(\nabla{\bf u})I-\frac{\mu}{K_{br}}{\bf u}
	\right]
$$
Here, $\varepsilon_P$ is the porosity of the material, and $K_{br}$ is the permeability of the bed.


#### species equation (= mass transport?)
$$
\frac{\partial (\rho Y_i)}{\partial t} + \nabla\cdot(\rho Y_i u_g) = \nabla\cdot(D_i\rho\nabla Y_i) + \nabla\cdot \left( \rho Y_i D_i^T + \frac{\nabla T}{T} \right) + Rxn_i
$$
<!---
$$
\frac{\partial c_i}{\partial t}+{\bf u}\cdot\nabla{\bf c}_i = \nabla \cdot(D_i\nabla c_i) + R_i
$$
where $c_i$, $D_i$, and $R_i$ is the concentration, the diffusion coefficient and the reaction rate of i-th species, respectively.
--->
where D is the diffusion constant, and $Rxn_i$ is the rate of reaction (defined in mass). $Y_i$ is the mass fraction of the species i, which is calculated with
$$
Y_k = \frac{s_k W_k}{\sum_k s_k W_k}
$$
where s_k is the site fraction ? and W_k is the molar mass.
* Gas density is often calculated from the ideal gas law.

* For the surface, we have to solve the site fraction conservation equation.
$$
\sigma \frac{\partial s_k}{\partial t} = RSk
$$
where $\sigma$ is the site density.

#### reaction
* In gas phase,
$$
\begin{split}
r_j = k_j \prod_m c_m^{\nu_{m,j}} \\
RRG_i = \sum \nu_{i,j} r_j W_j
\end{split}
$$
* On surface,
$$
\begin{split}
r_j = k_j \prod_m c_m^{\nu_{m,j}} \\
RRS_i = \sum \nu_{i,j} r_j W_j
\end{split}
$$
