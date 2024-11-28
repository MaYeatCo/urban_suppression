This is a repository that contains codes for paper "A study of pattern formation for across-diffusion system appearing in social applications: hotspot formation and suppression" by M. Yerlanov, N. Rodriguez, Q. Wang.

Mathematica codes:

Case1D and Case2D outlines steps of nonlinear weakly analysis: separating the expansion by powers of epsilons, solving PDE and obtaining the necessary coefficients. ToMatlab extension is needed to transfer results to Matlab, where numerical simulations are done, see below.

Matlab codes: PDE toolbox is needed and familiarity with its usage is preferred. 

urban_1D_single and urban_1D_test produces results for 1D case and allows to plot graphs in figures 1(b) and 1(a) (or 2(b) and 2(a)) respectively. Dw (the police diffusion) needs to be changed to move between super- and sub- critical regimes. Note that the codes need to be run a couple of times and results saved under different names if you want to compare two different initial conditions. These codes contain both ODE and PDE parts. Avoid running ODE in subcritical regime if you know that there is no stable steady state solution to converge to.

If you want to do the same thing for 2d (i.e. produce figures 3 and 4), you need to use urban_2D_ODE_single, urban_2D_PDE_single for a single simulation and urban_2d_ODE_test, urban_2D_PDE_test for a range of epsilons.  Dw (the police diffusion) needs to be changed to move between super- and sub- critical regimes. Note that it is preferred to run PDE codes first, as they generate mesh. This mesh (which is not a standard equispaced mesh) is used to compute RMS for ODE solution. Again one needs to run multiple times and save the relevant under different names to plot results for different initial conditions. Avoid running ODE in subcritical regime if you know that there is no stable steady state solution to converge to.

urban_2D_PDE_urho is a code to simulate crackdown strategy and produce figures 5, 6, 7. This code is similar to urban_2D_PDE_single, but the PDE needs to be solved three times. Dw (the police diffusion) needs to be changed to move between super- and sub- critical regimes, however note that if epsilon is positive then it does not really matter which regime you are considering.

urban_2D_PDE_chi is a code to simulate reallocation strategy and produce figures 8, 9. This code is similar to urban_2D_PDE_single. Dw (the police diffusion) needs to be changed to move between super- and sub- critical regimes, however note that this strategy works irrespectively of the regime or sign epsilon.
