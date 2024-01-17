All reject the null hypothesis


UPDATE: Sep 5
this version is reference for other test as I have made the following changes

-. Avoid the case where Pi_(t) < 0 due to random draw from normal rv. Change made on wr_bet.R

-. Normalize Pi_(t) so its sum across range states = 1, can be not equal to 1 due to computation error. Change made on wr_bet.R

-. Change made on wrapper_combined_wrbet. Now we save pA_theor_list from wr_bet.R as pA_theor_list_norm so need to change the name accordingly in this wrapper code

-. Change made on do_comparison_wrbet.R so now it can plot stacked bar plot over time instead of separate plot

Already did run this code (Sep 6)
