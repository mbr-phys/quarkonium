# Quarkonium Project

3rd Year Physics Problem Solving computing project to solve S.E. for quarkonium systems

## Project Plan

- [x] Program Design
    - [x] Solve SE for specific quark mass and energy values
    - [x] Calculate number of nodes and turning points
    - [x] Guess new energy or beta and iterate for solution
    - [x] Find energy or beta, whichever not found before
    - [x] Plot resulting wavefunctions

- [x] Milestone Program
    - [x] Calculate energy levels and wavefns of hydrogen
        - [x] Check for (n,l) = (1,0), (2,0), (2,1) states
        - [x] Compare to analytic solns
    - [x] ![My solution](/images/Hydro.png "My solution")
    - [x] Calculate energy of (2,1) state and plot |u(nl)|^2
        - ![Plot](/images/probs.png "Energies on plot as well")

- [x] Presentation
    - beamer slide show made with notes to read

- [x] Notate my code so I actually know what's going on
    - [x] possibly change up turning point counter, bit innaccurate right now

- [x] Program Extension
    - [x] Apply program to Charmonium - need to iterate over beta now instead of E
    - [x] Calculate beta to 3dp, using spin average mass
        - beta = 0.195
    - ![Plot for beta](/images/beta.png)
    - [x] Predict energy for (1,1) and (2,0)
        - ![Plot](/images/Charmonium.png)
    - [x] Produce plots of normalised wavefns for (1,0), (1,1), (2,0) states
        - ![Plot](/images/charmpdf.png)

- [x] Poster time
    - TeX poster using [this template](https://www.brian-amberg.de/uni/poster/) with some other people's modifications and my own

- [ ] Further Research
    - Extend spectrum to include (2,1) and (3,0)  
    - ![Plot](/images/fullspec.png)
    - Calculate spectrum using power law potential and `scipy.optimise` for better fitting
        - `scipy.optimize` doesn't quite work yet
    - ![Plot](/images/fullspec2.png)
    - Compare hyperfine splitting

