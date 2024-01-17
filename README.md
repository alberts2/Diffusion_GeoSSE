#### Code for analyses used on **A Diffusion-Based Approach for Simulating Forward-in-Time State-Dependent Speciation and Extinction Dynamics**

## Overview 

# *stationary time*  folder 

This folder consists of code used to compute time to stationary state frequencies in a two-region GeoSSE model (Section 2.8 of the manuscript). 

- **find_stationary_time.R**: R code to compute the time to stationary give the model's rate parameters. 
- **DATA**: A compilation of sets of model parameters used on the manuscript.
   - **Dataset_1**: model parameters used for Figure 7 (left panel).
   - **Dataset_2**: model parameters used for Figure 7 (right panel). 
   - **Dataset_3**: model parameters used for Figure 8 (left panel). 
   - **Dataset_4**: model parameters used for Figure 8 (right panel)
   
# *two-region GeoSSE* folder

This folder is used for analyses involving a two-region GeoSSE model and  it consists of the following code 

- **compute_params_2reg.R**	: generate values for model parameters according to sets of rules from a given set of state frequencies (Section 3.3 of the manuscript). 
- **compute_1repli_2reg.R**	: simulate one replicate of state dynamics given the model parameters. 
- **compute_1repli_2reg.R**	: a shell script for **compute_1repli_2reg.R** to generate multiple simulation replicates. 
- **wrapper_2reg.R**			: a wrapper script to combine all the replicates. 
- **plot_diff_2reg.R**			: plot the dynamics (Figures 7 and 8 of the manuscript). 

# *three-region GeoSSE* folder 

This folder is used for analyses involving a three-region GeoSSE model and it consists of the following sub-folders 

- **wrbet** folder					: consists of codes used for analyses using a GeoSSE model with only **within-region and between-region speciation events** included.
     - **wr_bet.R**					: generate one simulation replicate of state dynamics under a 3-region GeoSSE model with only **within-region and between-region speciation events** included.
     - **wr_bet.sh**					: a shell script to generate multiple replicates from **wr_bet.R**.
     - **wrapper_combined_wrbet.R**		: a wrapper script to combine all the replicates from **wr_bet.sh**.
     - **do_comparison_wrbet.R**		: plot the output (Figure 3) and do the statistical tests (Example 1, Table 1). 
- **wrdi** folder						: consists of codes used for analyses using a GeoSSE model with only **within-region speciation and dispersal events** included.
     - **wr_di_v2.R**					: generate one simulation replicate of state dynamics under a 3-region GeoSSE model with only **within-region speciation and dispersal events** included.
     - **wr_di_v2.sh**					: a shell script to generate multiple replicates from **wr_di_v2.R**.
     - **wrapper_combined_wrdi.R**		: a wrapper script to combine all the replicates from **wr_di_v2.sh**.
     - **do_comparison_wrdi.R**			: plot the output (Figure 4) and do the statistical tests (Example 2, Table 1). 
- **wrext** folder					: consists of codes used for analyses using a GeoSSE model with only **within-region speciation and extinction events** included.
     - **wr_ext.R**					: generate one simulation replicate of state dynamics under a 3-region GeoSSE model with only **within-region speciation and extinction events** included.
     - **wr_ext.sh**					: a shell script to generate multiple replicates from **wr_ext.R**.
      - **wrapper_combined_wrext.R**		: a wrapper script to combine all the replicates from **wr_ext.sh**.
      - **do_comparison_wrext.R**		: plot the output (Figure 4) and do the statistical tests (Example 3, Table 1). 
 - **geosse_all**	 folder			: consists of codes used for analyses using a GeoSSE model with **all the events** included.
     - **geosse_all.R**	 				: generate one simulation replicate of state dynamics under a full 3-region GeoSSE model. 
     - **geosse_all.sh**				: a shell script to generate multiple replicates from **geosse_all.R**.
     - **wrapper_combined_geosseall.R**	: a wrapper script to combine all the replicates from **geosse_all.sh**.
     - **do_comparison_geosseall.R**		: plot the output (Figure 6) and do the statistical tests (Example 4, Table 1). 

      
*Note: In each of these folders above, **MASTER Data** folder is used to store outputs from **MASTER Scripts** folder.* 

*Note: Inside **MASTER Scripts** folder, geosse.xml is used to generate one replicate of state dynamics using tree-based approach, and run_beast.sh is a shell script to simulate multiple replicates of geosse.xml output.*
 

