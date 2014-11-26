This collection of MATLAB scripts can be used to infer the parameters of a simple model of peptide presentation by MHC class I. Parameter inference is based on experimental time course measurements of intracellular peptide levels and cell-surface epitope levels. An example dataset is included in the files all_Kb_tryptic.csv (intracellular levels; scaled to have maximum value 1) and all_Kb_epitope.csv (cell-surface levels; absolute number of copies per cell).

*****************************************
***** Results ***************************
*****************************************
Results are saved to, and loaded from, the Results folder. An example results file (example_results.mat) is included. In general, results are given a filename containing the time and date of writing.

*****************************************
***** simplified_model_both_infer.m *****
*****************************************
To run the inference, run simplified_model_both_infer.m This script is also where you can choose the data file, the number of peptides being followed, the number of burn-in and sampling iterations and the initial state and step size for the sampler.

*****************************************
***** gibbs_sampler.m *******************
*****************************************
gibbs_sampler.m contains the MCMC sampler.

*****************************************
***** simplified_model.m ****************
*****************************************
The model itself, for intracellular peptide and cell-surface epitope levels, is contained in simplified_model.m. ODEs are solved using the SUNDIALS toolbox (see http://computation.llnl.gov/casc/sundials/main.html), which must be installed for the function to run.

*****************************************
***** LL.m ******************************
*****************************************
The log-likelihood of the data given some proposed parameter values is calculated using the function LL. Inputs are simulated data given current parameter values, data (intracellular and cell-surface), and a two-element vector [intracellular peptide error, tryptic error].

*****************************************
***** plot_results.m ********************
*****************************************
plot_results is a script to plot saved results. First load the results from the Results folder, then run the script to generate plots.

Plots of the results in the example file indicate problems with the inference including parameter collinearities and sample autocorrelation. The output should therefore be regarded with caution, and the code needs some work before it can be used for reliable inference.

Some particular improvements that could be made would include:
	* More informative command-line output.
	* Shortening the sampling code to include only one MCMC loop, with an "if" to differentiate burn-in and inference stages.
	* Imposing some convergence criteria which should be met for the burn-in stage to end and inference to begin. 
	* The option to impose prior probability distributions on parameters.
	* A mixed-effects (hierarchical) model of parameter values
	* More efficient sampling algorithms.

Joanna Lewis 
26 November 2014