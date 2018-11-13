# Project Plan

## 1. Program Design

### 1. Solve SE for specific quark mass and energy values

What is beta??

### 2. Calculate number of nodes and turning points

make specific functions for it

### 3. Guess new energy or beta and iterate for solution

### 4. Find energy or beta, whichever not found before

### 5. Plot resulting wavefunctions

## 2. Milestone Program

### Calculate energy levels and wavefns of hydrogen

#### Check for (n,l) = (1,0), (2,0), (2,1) states

#### Compare to analytic solns

### Calculate energy of (2,1) state and plot |u(nl)|^2

## 3. Program Extension

### Apply program to Charmonium

#### Calculate beta to 3dp, using spin average mass

#### Predict energy for (1,1) and (2,0)

#### Produce plots of normalised wavefns for (1,0), (1,1), (2,0) states

## 4. Further Research

## Context etc from Ben

- Quarkonium helps with QCD 
    - quantum field theory (qft) of the strong force
    - strong force binds protons and neutrons and nuclei etc
    - QCD in terms of quarks and gluons 
- helpful to compare QCD with QED
    - qft of emag
    - associated with things like electrons, positrons, and photons

### QED
- Feynman diagrams for exchange of photons between particles builds up classical Coulomb potential

$$
    V = -\frac{\alpha}{r}, \alpha \approx \frac{1}{137}
$$
- can calculate in QED using perturbation theory based on alpha being small

### QCD

- up, down, strange, charm, top, bottom quarks, and antiquarks like electrons and positrons
- gluons
    - like the photon in the sense of being a force carrier
    - different as it also couples to itself, i.e. there are forces between gluons
    - has chromomagnetic charge, $g_s$

$$
    \alpha_s = \frac{g_s^2}{4\pi} = \alpha(\mu)
$$

- $\mu$ is an energy scale

$$
    \begin{aligned}
    \alpha_s(\mu) &<< 1, \mu \geq 2\,\text{GeV} \\
    \alpha_s(\mu) &\approx O(1), \mu \approx \Lambda_{QCD} \approx 500\,\text{MeV}
    \end{aligned}
$$

- approx 1/r graph of alpha vs mu
- two regimes of QCD
    - infrared slavery, small mu, alpha of order 1 or higher
    - asymptotic freedom, large mu, small alpha, similar to QED can use perturbation theory

#### Heavy and light quarks

##### Light Quarks

- $m_s \approx 100\,MeV < \Lambda_{QCD}$
- $m_{u,d} \approx 2\,MeV < \Lambda_{QCD_{}}$
- bit more fiddly to solve with as they can't be solved using perturbation theory etc as alpha is very large

##### Heavy quarks

- $m_c \approx 1.5\,GeV \geq \Lambda_{QCD}$
- $m_b \approx 5\,GeV >> \Lambda_{QCD}$
- $m_t \approx 175\,GeV >>> \Lambda_{QCD_{}}$, similar mass to gold atom, decays too fast to form bound states

- Interested in charm and bottom bound states as they're heavy but not too heavy
- Charmonium and bottomonium are forms of quarkonium

#### Models for interquark potentials in quarkonium

Mass of a bound state:
$$
    M_{nl} = 2m_q + E_{nl}
$$
This is a non-relativistic approximation, holds if $E_{nl} << m_{q}$.

Non-relativistic means can use standard QM, i.e. Schrodinger eqn 

$$
    \left[-\frac{\hbar^2}{2m}\vec{\nabla}^2 + V(r)\right]\Psi_{nl} = E_{nl}\Psi_{nl}
$$
