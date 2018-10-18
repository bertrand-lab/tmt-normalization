These are some functions for doing the Internal Reference Standard / Sample Loading / TMM (EdgeR) Normalization.

The main function is `sl_irs_tmm_normalization()`, which also does the plotting showing how ion intensities shift by channel for each step.

The three inputs are 1) `protein_df`, which is the output of protein intensities. If you have multiple injections, perhaps the average of intensities should be taken beforehand. 2) `tmt_exp_columns`, which is a list. The list should be structured like this: `list(c(1, 2, 3), c(4, 5, 6))` where the numbers refer to the columns in `protein_df` separating different TMT experiments. The last thing is `tmt_common_channel_names`, which designates the names of internal reference standards within each TMT experiment.

Note that the original code was modified (into functions) from an excellent guide by Phil Wilmart (https://github.com/pwilmart/IRS_normalization). The main thing these scripts do is wrap into functions, allowing a variable number of TMT experiments.

The other key difference is that we do not take the _sum_ of IRS channels, rather the _mean_. This is in case there are different numbers of IRS channels (if for example you have one TMT experiment with just one IRS channel). 
