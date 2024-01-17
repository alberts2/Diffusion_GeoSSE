#### Code for analyses used on **A Diffusion-Based Approach for Simulating Forward-in-Time State-Dependent Speciation and Extinction Dynamics**

## Overview 

# **``stationary time"**  folder 

This folder consists of code used to compute time to stationary state frequencies in a two-region GeoSSE model (Section 2.8 of the manuscript). 

- **find_stationary_time.R**: R code to compute the time to stationary give the model's rate parameters. 
- **DATA**: A compilation of sets of model parameters used on the manuscript.
   - **Dataset_1**: model parameters used for Figure 7 (left panel).
   - **Dataset_2**: model parameters used for Figure 7 (right panel). 
   - **Dataset_3**: model parameters used for Figure 8 (left panel). 
   - **Dataset_4**: model parameters used for Figure 8 (right panel)
   
# **``two-region GeoSSE"** folder

This folder consist of the following code 

- **compute_params_2reg.R**	: generate values for model parameters according to sets of rules from a given set of state frequencies (Section 3.3 of the manuscript). 
- **compute_1repli_2reg.R**	: simulate one replicate of state dynamics given the model parameters. 
- **compute_1repli_2reg.R**	: a shell command for **compute_1repli_2reg.R** to generate multiple simulation replicates. 
- **wrapper_2reg.R**			: a wrapper code to combine all the replicates. 
- **plot_diff_2reg.R**			: plot the dynamics (Figures 7 and 8 of the manuscript). 
