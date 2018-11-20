# Notes for Pres

## Slide 1

- Study bound states then expand project and explore into further research

- Systems of quarks can be studied and solved using QCD, similar to QED for protons and electrons

- QCD has the strong force like QED has the electromagnetic force

- Charge carriers are the gluons - but unlike photons for EM, forces between gluons

- Quarks have colour charge, so from a mechanic of QCD known as colour confinement, they can only exist in bound systems and not on their own, must form systems of overall 0 colour charge

- 6 different "flavours" of quark - up, down, strange, charm, top, bottom

- Three generations split roughly by masses 
    1. up and down
    2. charm and strange
    3. top and bottom
- up, charm, top have same electric charge 2/3
- down, strange, bottom have same electric charge -1/3
- up and down have plus/minus half isospin, the rest zero
- four other properties defined for quarks - charm, strangeness, topness, bottomness
- values all zero apart from the eponymous flavours

## Slide 2

- Has two regimes, asymptotic freedom and infrared slavery

- Interested in asymptotic freedom, where the interactions between particles as the energy scale increases, allowing us to use perturbative calculations on them

- This means we can only consider higher mass quarks for this study, as lower mass ones are bound strongly together, e.g into a proton

- Toponium decays too quickly to be studied, so interested in the bound states of charm-anticharm and bottom-antibottom, known as charmonium and bottomonium respectively

- Charmonium and bottomonium can then be modeled similarly to a hydrogen atom, but using the strong force of QCD instead of EM from QED

- Leads us to use the 3D Schrodinger equation, using the Cornell potential as a model for quark potential compared to the EM potential in a hydrogen atom

- Cornell potential is just one of many models for interquark potential, but one of the most popular

## Slide 3

- The whole wavefunction is $r^2u_{nl}Y_{lm}$, but only interested in the radial function for this

- From the Schrodinger equation, can end up with a set of ODEs to be solved for the radial wavefunction

- We require that unl(0) = 0, and vnl(0) = 1, not caring for normalisation at this stage

- Then can solve this set of ODEs using `scipy.integrate`'s `odeint` function for unl

- Given either a set energy or beta to begin with, we want to iterate over values of other to find the correct solution for the values of n and l

- Guess three initial values, with second value being the mean of the outer two
- Find solution for three values, then evaluate expected number of nodes and turning points
- If nodes and tpts change between values, then set second value to outer value and new half way point between this and one of the other ones
- Continue until nodes and tpts converge
- If nothing in this range, try again over a different range
- Normalise the solution as normal

## Slide 4

- Milestone project to adapt program to values for hydrogen and solve its wavefunctions for (1,0),(2,0), and (2,1)

- Chose energies to iterate over to find their values once beta was known to be zero

- Energies found are quite close to actual values, E1 especially so

- The nature of the program is such that at larger values it will diverge to plus/minus infinity, but can solve up to a point

## Slide 5

- Looking forward with a functioning program for hydrogen, will be applying program to charmonium and bottomonium and modelling their wavefunctions

- Then looking to extend research into other others

1. Cornell potential just one way to model interquark potential, could look at other potential models and see what solutions these give
2. Could look at more exotic types of quark systems and see how these can be solved, such as a system with a charm and a bottom
3. Wavefunctions also vary in time, could look at the evolution of the system, and try to calculate an accurate figure of the lifetime of these systems as they decay
    - Possibly even try this with toponium, to confirm its short life
4. Similar to in a Helium atom, the quarks having spin half means they can be in different energy states inside n,l energy levels with slightly different values
5. Charmonium is expected to melt in high temperatures as a signal of the formation of a quark-gluon plasma. Unsure on feasability, but perhaps studying this "melting" and exploring high temperature interactions of these systems into a plasma.

