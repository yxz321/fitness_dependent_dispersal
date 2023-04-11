# fitness_dependent_dispersal

Code for journal article "Maintenance of biodiversity in multitrophic metacommunities: dispersal mode matters" (Xiaozhou Ye &amp; Shaopeng Wang, 2023), with MATLAB.

Simulated data can be found in Zenodo: https://doi.org/10.5281/zenodo.7806048

* To repeat figure plotting for (most) figures: run the `%% draw figures` section in `main.m`
  * to be able to plot Fig. S1-3 and figures in supplementary text, first download /rawdata/, randPars_s50p3r1000_0218062416.mat, and randPars_s50p6r100_0815225637.mat from Zenodo and put in working directory
* To repeat data analysis (diversity calculation): download /rawdata/, randPars_s50p3r1000_0218062416.mat, and randPars_s50p6r100_0815225637.mat from Zenodo and put in working directory, then run the `%% data synthesis` section in `main.m`
* To repeat simulation (takes a long time) and subsequent data analysis:
  * from the same random foodwebs as in publication: download randPars_s50p3r1000_0218062416.mat, and randPars_s50p6r100_0815225637.mat from Zenodo and put in working directory, then run the entire `main.m`
  * from newly generated random foodwebs: change the `%% generate or load randPars (random basal foodwebs)` section in `sim_scenarios.m` by uncommenting `% if generates new random basal foodwebs` snippet and comment `% if use existing basal foodwebs` snippet, then run the entire `main.m`

## functions

The most important functions:
* `main.m`: main function for simulation, analysis, and plotting
* `odefunc.m`: define the model as ODE
* `sim_scenarios.m`
* `parasetgenerator.m`: generate random foodweb structures, together with parameters

### `main.m`

Main function, run the entire procedure of: simulation -> diversity calculation -> raw figure creation. Change `snr` to get results with modified plant dispersal scenario, results for more patches, and results without tendency to stay (epsilon=0).

### `paraSet = parasetgenerator(nrep, varargin)`

Generate #nrep of random foodwebs, each with different structure and parameter.

* `fw = foodwebgenerator(S,C)`

  Input parameters
  * S - number of species
  * C - network connectance of foodweb

  Output parameters
  ( fw - S*S adjacent matrix of predatory relationship

### `main_simulate(snr, vary, repList)`

Run simulation and save results. `vary` can take {'H', 'm0'}, indicating varying habitat heterogeneity (H) or connectivity (m0).

Define combinations of H, m0 and $\lambda$ (fitness-dependency), then loop through all replicates and call `sim_scenarios.m`. Each replicate represents a different basal foodweb structure. 

* `dBdt = odefunc(B,par)`

  * The ODE function describing the model.

* `[BfinalSum, BfinalCVSum] = sim_scenarios(parBase, lambList, mList, heterList, pltMigStr)`

  * For a given replicate, loop through defined combinations of H, m0 and $\lambda$, call `ode15s_withevent.m` to solve the equations and return equilibrium species abundance for metacommunity.

* `[res_t,res_b] = ode15s_withevent(B0,par,tspan,odeopt,varargin)`

  * Call ODE solver `ode15s_md.m` while dealing with species extinction event.

* `ode15s_md.m`

  * Originated from `ode15s.m` from MATLAB, slightly modified to stop integration and return 'TimeLimitReached' error when the adaptive solver has run for >10000 integration steps. The error occurs only rarely (~1%), when a wildly fluctuating system (e.g. chaotic) is encountered.
 
* functions in /funcs_from_matlab/
  * required for ode15s_md.m

* `extinEventFcn(t,B,par)`

  * The event function used in `ode15s_md.m`

### `processdata(snr)`

Calculate $\alpha$, $\beta$, and $\gamma$ diversity across replicates and for both varying H and varying m0 simulations. Save result as `sprintf('./%sres.mat',snr)`. 

### `drawfig(snr)`

Draw Fig. 1-5, and Fig. S4-14.

### `drawfig_sup()`

Draw Figs. S1-S3, and figures in supplementary text. 

* `plotfoodweb`
  * `trophiclevel`
  * `fwcylinder`
  * `fwmds`
  * `fwrotation`
  * `patchline`

## Datafiles

Available in Github: 'res.mat', 'p6_res.mat', 'pltiso_res.mat', 'pltrand_res.mat', './rawdata/%srep%03d_diff_H.mat'
* summary of \alpha, \beta, and \gamma diversity for different scenarios. Sufficient for generating Figs. 1-5 and Fig. S4-14.

Available only in Zenodo: 'randPars_s50p3r1000_0218062416.mat', 'randPars_s50p6r100_0815225637.mat', \rawdata\
