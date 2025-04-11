"""
    Write the results of the estimation as table code into a series of .tex-file.
    This script generates all the table files from the MATLAB simulation output.

"""
import csv
import numpy as np
import pandas as pd
from bld.project_paths import project_paths_join
pd.options.mode.chained_assignment = None  # default='warn'

# Set a different types of table row with placeholders.
table_row_header = ' & {val1} & {val2} & {val3} & {val4}  \\tabularnewline \n'
table_row_var = """{var} & {val0:.2f} & {val1:.4f}{str1} & {val2: .4f}{str2} & {val3: .4f}{str3} \\tabularnewline\n & & ({se1:.4f}) & ({se2:.4f}) & ({se3:.4f}) \\tabularnewline\n"""
table_row_likelihood = """{var} & & {val1:.1f}{str1} & {val2: .1f}{str2} & {val3: .1f}{str3} \\tabularnewline\n"""
table_row_SSres = """{var} & & {val1:.4f}{str1} & {val2: .4f}{str2} & {val3: .4f}{str3} \\tabularnewline\n"""

# Write Table D.1.
ols_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'ols_table_data.csv'))
# Write the results to a LaTeX table.
with open(project_paths_join('OUT_TABLES', 'estOLS_L3_new.tex'), 'w') as tex_file:
    # Top of table.
    tex_file.write('\\begin{tabular}{lc|ccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=50,000$', val3='$N=100,000$',val4='$N=1,000,000$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    for row_idx, row in ols_data.iterrows():
        if row_idx<=1:
            row['true']
            tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=ols_data['beta_hat_50000'][row_idx], val2=ols_data['beta_hat_100000'][row_idx],val3=ols_data['beta_hat_1000000'][row_idx], str1='',str2='',str3='',se1=ols_data['sd_50000'][row_idx],se2=ols_data['sd_100000'][row_idx],se3=ols_data['sd_1000000'][row_idx]))
        if row_idx == len(ols_data)-1:
            tex_file.write(table_row_SSres.format(var=row['Row'], val1=ols_data['beta_hat_50000'][row_idx], val2=ols_data['beta_hat_100000'][row_idx],val3=ols_data['beta_hat_1000000'][row_idx], str1='',str2='',str3=''))
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{14cm}}{\\footnotesize{\\textit{Notes: The parameters are estimated using our OLS estimator.  Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')

# Write Table D.2.
md_sophi_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'md_table_data.csv'))
with open(project_paths_join('OUT_TABLES', 'estMD_sophi_new.tex'), 'w') as tex_file:
    # Top of table.
    tex_file.write('\\begin{tabular}{lc|ccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=50,000$', val3='$N=100,000$',val4='$N=1,000,000$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    for row_idx, row in md_sophi_data.iterrows():
        row['true']
        if row_idx == len(md_sophi_data)-1:
            tex_file.write(table_row_SSres.format(var=row['Row'], val1=md_sophi_data['beta_hat_50000'][row_idx], val2=md_sophi_data['beta_hat_100000'][row_idx],val3=md_sophi_data['beta_hat_1000000'][row_idx], str1='',str2='',str3=''))
        else:
            tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=md_sophi_data['beta_hat_50000'][row_idx], val2=md_sophi_data['beta_hat_100000'][row_idx],val3=md_sophi_data['beta_hat_1000000'][row_idx], str1='',str2='',str3='',se1=md_sophi_data['sd_50000'][row_idx],se2=md_sophi_data['sd_100000'][row_idx],se3=md_sophi_data['sd_1000000'][row_idx]))
        
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{14cm}}{\\footnotesize{\\textit{Notes: The parameters are estimated using our ROLS estimator. Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')

# Write Table D.6.
md_rn_sophi_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'md_rn_table_data.csv'))
with open(project_paths_join('OUT_TABLES', 'estMD_rn_sophi_new.tex'), 'w') as tex_file:
    # Top of table.
    tex_file.write('\\begin{tabular}{lc|ccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=50,000$', val3='$N=100,000$',val4='$N=1,000,000$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    for row_idx, row in md_rn_sophi_data.iterrows():
        row['true']
        if row_idx == len(md_rn_sophi_data)-1:
            tex_file.write(table_row_SSres.format(var=row['Row'], val1=md_rn_sophi_data['beta_hat_50000'][row_idx], val2=md_rn_sophi_data['beta_hat_100000'][row_idx],val3=md_rn_sophi_data['beta_hat_1000000'][row_idx], str1='',str2='',str3=''))
        else:
            tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=md_rn_sophi_data['beta_hat_50000'][row_idx], val2=md_rn_sophi_data['beta_hat_100000'][row_idx],val3=md_rn_sophi_data['beta_hat_1000000'][row_idx], str1='',str2='',str3='',se1=md_rn_sophi_data['sd_50000'][row_idx],se2=md_rn_sophi_data['sd_100000'][row_idx],se3=md_rn_sophi_data['sd_1000000'][row_idx]))
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{14cm}}{\\footnotesize{\\textit{Notes: The parameters are estimated using our ROLS estimator with exclusion restriction. Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    tex_file.write('\\bottomrule\n\\end{tabular}\n')

# Write Table D.8.
md_naive_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'md_table_data_naive.csv'))
with open(project_paths_join('OUT_TABLES', 'estMD_naive_new.tex'), 'w') as tex_file:
    # Top of table.
    tex_file.write('\\begin{tabular}{lc|ccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=50,000$', val3='$N=100,000$',val4='$N=1,000,000$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    for row_idx, row in md_naive_data.iterrows():
        row['true']
        if row_idx == len(md_naive_data)-1:
            tex_file.write(table_row_SSres.format(var=row['Row'], val1=md_naive_data['beta_hat_50000'][row_idx], val2=md_naive_data['beta_hat_100000'][row_idx],val3=md_naive_data['beta_hat_1000000'][row_idx], str1='',str2='',str3=''))
        else:
            tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=md_naive_data['beta_hat_50000'][row_idx], val2=md_naive_data['beta_hat_100000'][row_idx],val3=md_naive_data['beta_hat_1000000'][row_idx], str1='',str2='',str3='',se1=md_naive_data['sd_50000'][row_idx],se2=md_naive_data['sd_100000'][row_idx],se3=md_naive_data['sd_1000000'][row_idx]))
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{14cm}}{\\footnotesize{\\textit{Notes: The parameters are estimated using our MD estimator. Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    tex_file.write('\\bottomrule\n\\end{tabular}\n')

# Write Table D.3.
mle_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'mle_table_data.csv'))
with open(project_paths_join('OUT_TABLES', 'estMLE_new.tex'), 'w') as tex_file:
    # Top of table.
    tex_file.write('\\begin{tabular}{lc|ccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=10,000$', val3='$N=20,000$',val4='$N=30,000$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    for row_idx, row in mle_data.iterrows():
        row['true']
        if row_idx == len(mle_data)-1:
            tex_file.write(table_row_likelihood.format(var=row['Row'], val1=mle_data['beta_hat_10000'][row_idx], val2=mle_data['beta_hat_20000'][row_idx],val3=mle_data['beta_hat_30000'][row_idx], str1='',str2='',str3=''))
        else:
            tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=mle_data['beta_hat_10000'][row_idx], val2=mle_data['beta_hat_20000'][row_idx],val3=mle_data['beta_hat_30000'][row_idx], str1='',str2='',str3='',se1=mle_data['sd_10000'][row_idx],se2=mle_data['sd_20000'][row_idx],se3=mle_data['sd_30000'][row_idx]))
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{12cm}}{\\footnotesize{\\textit{Notes: The parameters are estimated using MLE.  Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')

# Write Table D.9.
mle_naive_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'mle_naive_table_data.csv'))
# Write the results to a LaTeX table.
with open(project_paths_join('OUT_TABLES', 'estMLE_naive_new.tex'), 'w') as tex_file:
    tex_file.write('\\begin{tabular}{lc|ccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=10,000$', val3='$N=20,000$',val4='$N=30,000$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    for row_idx, row in mle_naive_data.iterrows():
        row['true']
        if row_idx == len(mle_naive_data)-1:
            tex_file.write(table_row_likelihood.format(var=row['Row'], val1=mle_naive_data['beta_hat_10000'][row_idx], val2=mle_naive_data['beta_hat_20000'][row_idx],val3=mle_naive_data['beta_hat_30000'][row_idx], str1='',str2='',str3=''))
        else:
            tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=mle_naive_data['beta_hat_10000'][row_idx], val2=mle_naive_data['beta_hat_20000'][row_idx],val3=mle_naive_data['beta_hat_30000'][row_idx], str1='',str2='',str3='',se1=mle_naive_data['sd_10000'][row_idx],se2=mle_naive_data['sd_20000'][row_idx],se3=mle_naive_data['sd_30000'][row_idx]))
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{12cm}}{\\footnotesize{\\textit{Notes:The parameters are estimated using MLE. Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')

# Write Table D.5.
mle_norm_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'mle_norm_table_data.csv'))
table_row_var = """{var} & {val0:.2f} & {val1:.4f}{str1} & {val2: .4f}{str2} & {val3: .4f}{str3} \\tabularnewline\n & & ({se1:.1e}) & ({se2:.1e}) & ({se3:.1e}) \\tabularnewline\n"""
with open(project_paths_join('OUT_TABLES', 'estMLE_norm_new.tex'), 'w') as tex_file:
    # Top of table.
    tex_file.write('\\begin{tabular}{lc|ccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=10,000$', val3='$N=20,000$',val4='$N=30,000$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    for row_idx, row in mle_norm_data.iterrows():
        row['true']
        if row_idx == len(mle_norm_data)-1:
            tex_file.write(table_row_likelihood.format(var=row['Row'], val1=mle_norm_data['beta_hat_10000'][row_idx], val2=mle_norm_data['beta_hat_20000'][row_idx],val3=mle_norm_data['beta_hat_30000'][row_idx], str1='',str2='',str3=''))
        else:
            tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=mle_norm_data['beta_hat_10000'][row_idx], val2=mle_norm_data['beta_hat_20000'][row_idx],val3=mle_norm_data['beta_hat_30000'][row_idx], str1='',str2='',str3='',se1=mle_norm_data['sd_10000'][row_idx],se2=mle_norm_data['sd_20000'][row_idx],se3=mle_norm_data['sd_30000'][row_idx]))
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{12cm}}{\\footnotesize{\\textit{Notes:  The parameters are estimated using MLE with the flow utility of action 0 being zero and exclusion restriction. Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')

# Write Table D.4.
mle_normTerm_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'mle_normTerm_table_data.csv'))
with open(project_paths_join('OUT_TABLES', 'estMLE_normTerm_new.tex'), 'w') as tex_file:
    # Top of table.
    tex_file.write('\\begin{tabular}{lc|ccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=10,000$', val3='$N=20,000$',val4='$N=30,000$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    for row_idx, row in mle_normTerm_data.iterrows():
        row['true']
        if row_idx == len(mle_normTerm_data)-1:
            tex_file.write(table_row_likelihood.format(var=row['Row'], val1=mle_normTerm_data['beta_hat_10000'][row_idx], val2=mle_normTerm_data['beta_hat_20000'][row_idx],val3=mle_normTerm_data['beta_hat_30000'][row_idx], str1='',str2='',str3=''))
        else:
            tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=mle_normTerm_data['beta_hat_10000'][row_idx], val2=mle_normTerm_data['beta_hat_20000'][row_idx],val3=mle_normTerm_data['beta_hat_30000'][row_idx], str1='',str2='',str3='',se1=mle_normTerm_data['sd_10000'][row_idx],se2=mle_normTerm_data['sd_20000'][row_idx],se3=mle_normTerm_data['sd_30000'][row_idx]))
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{12cm}}{\\footnotesize{\\textit{Notes:  The parameters are estimated using MLE with the flow utility of action 1 being zero and exclusion restriction. Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')

table_row_var = """{var} & {val0:.2f} & {val1:.4f}{str1} & {val2: .4f}{str2} & {val3: .4f}{str3} \\tabularnewline\n & & ({se1:.4f}) & ({se2:.4f}) & ({se3:.4f}) \\tabularnewline\n"""

# Write Table D.7.
mle_rn_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'mle_rn_table_data.csv'))
with open(project_paths_join('OUT_TABLES', 'estMLE_rn_new.tex'), 'w') as tex_file:
    # Top of table.
    tex_file.write('\\begin{tabular}{lc|ccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=10,000$', val3='$N=20,000$',val4='$N=30,000$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    for row_idx, row in mle_rn_data.iterrows():
        row['true']
        if row_idx == len(mle_rn_data)-1:
            tex_file.write(table_row_likelihood.format(var=row['Row'], val1=mle_rn_data['beta_hat_10000'][row_idx], val2=mle_rn_data['beta_hat_20000'][row_idx],val3=mle_rn_data['beta_hat_30000'][row_idx], str1='',str2='',str3=''))
        else:
            tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=mle_rn_data['beta_hat_10000'][row_idx], val2=mle_rn_data['beta_hat_20000'][row_idx],val3=mle_rn_data['beta_hat_30000'][row_idx], str1='',str2='',str3='',se1=mle_rn_data['sd_10000'][row_idx],se2=mle_rn_data['sd_20000'][row_idx],se3=mle_rn_data['sd_30000'][row_idx]))
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{12cm}}{\\footnotesize{\\textit{Notes: The parameters are estimated using MLE imposing the exclusion restriction. Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')

# Set a different types of table row with placeholders.
table_row_header = ' & {val1} & {val2} & {val3} & {val4} \\tabularnewline \n'
table_row_var = """{var} & {val0:.2f} & {val1:.4f}{str1} & {val2: .4f}{str2} & {val3: .4f}{str3} \\tabularnewline\n & & ({se1:.4f}) & ({se2:.4f}) & ({se3:.4f}) \\tabularnewline\n"""

# Write Table D.10.
mle_rn_eD_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'mle_rn_eD_table_data.csv'))
# Write the results to a LaTeX table.
with open(project_paths_join('OUT_TABLES', 'estMLE_rn_eD_new.tex'), 'w') as tex_file:
    # Top of table.
    tex_file.write('\\begin{tabular}{lc|cccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=20,000$', val3='$N=20,000 (Sophi)$',val4='$N=20,000(Naive)$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    for row_idx, row in mle_rn_eD_data.iterrows():
        row['true']
        if row_idx == len(mle_rn_eD_data)-1:
            tex_file.write(table_row_likelihood.format(var=row['Row'], val1=mle_rn_eD_data['beta_hat_20000'][row_idx], val2=mle_rn_eD_data['beta_hat_sophi_20000'][row_idx],val3=mle_rn_eD_data['beta_hat_naive_20000'][row_idx], str1='',str2='',str3=''))
        else:
            tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=mle_rn_eD_data['beta_hat_20000'][row_idx], val2=mle_rn_eD_data['beta_hat_sophi_20000'][row_idx],val3=mle_rn_eD_data['beta_hat_naive_20000'][row_idx], str1='',str2='',str3='',str4='',se1=mle_rn_eD_data['sd_20000'][row_idx],se2=mle_rn_eD_data['sd_sophi_20000'][row_idx],se3=mle_rn_eD_data['sd_naive_20000'][row_idx]))
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{15cm}}{\\footnotesize{\\textit{Notes: The third column set $\\beta$ to be 1 (and thus standard error is always zero). The fourth column (with Sophi in parentheses) estimates the sophisticated agent model. The fifth column (with Naive in parentheses) estimates the naive agent model. All parameters are estimated using MLE imposing the exclusion restriction. Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')