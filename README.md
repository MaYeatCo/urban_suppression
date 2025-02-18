This is a repository that contains codes for the paper "Formation and Suppression of Hotspots in Urban Crime
Models with Law Enforcement" by M. Yerlanov, N. Rodriguez, Q. Wang.

Mathematica codes:

Case1D and Case2D outline steps of nonlinear weakly analysis: separating the expansion by powers of epsilons, solving PDE and obtaining the necessary coefficients. ToMatlab extension is needed to transfer results to Matlab, where numerical simulations are done, see below.

Matlab codes: a PDE toolbox is needed and familiarity with its usage is preferred. 

urban_1D_single and urban_1D_test produce results for the 1D case and allow to plot graphs in figures 1(b) and 1(a) (or 3(b) and 3(a)) respectively, plus with small modification you can draw the error plot as well, see figure 2. Dw (the police diffusion) needs to be changed to move between super- and sub- critical regimes. Note that the codes need to be run a couple of times and results saved under different names if you want to compare two different initial conditions. These codes contain both ODE and PDE parts. Avoid running ODE in subcritical regime if you know that there is no stable steady-state solution to converge to.

Figure 1:

<img src="./pics/1deps_super_v3.png" alt="drawing" width="200"/><img src="./pics/1dsingle_super_v3.png" alt="drawing" width="200"/>

Figure 2:

<img src="./pics/1d_eps_super_error.png" alt="drawing" width="200"/>

Figure 3:

<img src="./pics/1deps_sub_v3.png" alt="drawing" width="200"/><img src="./pics/1dsingle_sub_v3.png" alt="drawing" width="200"/>

If you want to do the same thing for 2d (i.e. produce figures 3 and 4), you need to use urban_2D_ODE_single, urban_2D_PDE_single for a single simulation and urban_2d_ODE_test, urban_2D_PDE_test for a range of epsilons.  Dw (the police diffusion) needs to be changed to move between super- and sub- critical regimes. Note that it is preferred to run PDE codes first, as they generate mesh. This mesh (which is not a standard equispaced mesh) is used to compute RMS for ODE solution. Again one needs to run multiple times and save the relevant under different names to plot results for different initial conditions. Avoid running ODE in subcritical regime if you know that there is no stable steady-state solution to converge to.

Figure 4:

<img src="./pics/2deps_super_v5.png" alt="drawing" width="200"/><img src="./pics/2dsingle_super_v4.png" alt="drawing" width="200"/>

Figure 5:

<img src="./pics/2deps_sub_v3.png" alt="drawing" width="200"/><img src="./pics/2dsingle_sub_v3.png" alt="drawing" width="200"/>

urban_2D_PDE_urho is a code to simulate the crackdown strategy and produce figures 6, 7, 8, 9. This code is similar to urban_2D_PDE_single, but the PDE needs to be solved three times. Dw (the police diffusion) needs to be changed to move between super- and sub- critical regimes, however, note that if epsilon is positive then it does not really matter which regime you are considering.

Figure 6:

<img src="./pics/supp_urho_pos_v2.png" alt="drawing" width="200"/><img src="./pics/supp_urho_neg_v2.png" alt="drawing" width="200"/>

Figure 7:

<img src="./pics/supp_urho_pos_1_v2.png" alt="drawing" width="200"/><img src="./pics/supp_urho_pos_2_v2.png" alt="drawing" width="200"/><img src="./pics/supp_urho_pos_3_v2.png" alt="drawing" width="200"/>

Figure 8:

<img src="./pics/supp_urho_neg_1_v2.png" alt="drawing" width="200"/><img src="./pics/supp_urho_neg_2_v2.png" alt="drawing" width="200"/><img src="./pics/supp_urho_neg_3_v2.png" alt="drawing" width="200"/>

Figure 9:

<img src="./pics/supp_urtho_pos_sub_v1.png" alt="drawing" width="200"/><img src="./pics/supp_urtho_pos_super_v1.png" alt="drawing" width="200"/>

urban_2D_PDE_chi is a code to simulate a reallocation strategy and produce figures 10, 11. This code is similar to urban_2D_PDE_single. Dw (the police diffusion) needs to be changed to move between super- and sub- critical regimes, however, note that this strategy works irrespectively of the regime or sign epsilon.

Figure 10:

<img src="./pics/supp_chi_pos_v2.png" alt="drawing" width="200"/>

Figure 11:

<img src="./pics/supp_chi_pos_1_v2.png" alt="drawing" width="200"/><img src="./pics/supp_chi_pos_2_v2.png" alt="drawing" width="200"/><img src="./pics/supp_chi_pos_3_v2.png" alt="drawing" width="200"/>

urban_2D_PDE_arho is a code to simulate the restoration strategy and produce figure 12. This code is similar to urban_2D_PDE_single, but the PDE needs to be solved three times. Dw (the police diffusion) needs to be changed to move between super- and sub- critical regimes, however, note that if epsilon is positive then it does not really matter which regime you are considering.

Figure 12:

<img src="./pics/supp_arho_pos_v2.png" alt="drawing" width="200"/><img src="./pics/supp_arho_neg_v2.png" alt="drawing" width="200"/>
