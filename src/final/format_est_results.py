"""

    Write the results of the simulation/estimation as table code into .tex-files.

    Input estimation files:
    - mle_norm_table_data.csv      
    - mle_normTerm_table_data.csv      
    - mle_table_data.csv
    - mle_naive_table_data.csv   
    - md_table_data_naive.csv      
    - ols_table_data.csv

    Output tables:
    - estMD_naive_new.tex
    - estOLS_L3_new.tex
    - estMLE_new.tex
    - estMLE_norm_new.tex
    - estMLE_normTerm_new.tex
    - estMLE_naive_new.tex
"""

import csv
import numpy as np
import pandas as pd
from bld.project_paths import project_paths_join
pd.options.mode.chained_assignment = None  # default='warn'

# Define different types of table row with placeholders.
table_row_header = ' & {val1} & {val2} & {val3} & {val4}  \\tabularnewline \n'
table_row_var = """{var} & {val0:.2f} & {val1:.4f}{str1} & {val2: .4f}{str2} & {val3: .4f}{str3} \\tabularnewline\n & & ({se1:.4f}) & ({se2:.4f}) & ({se3:.4f}) \\tabularnewline\n"""

# Load estimation results.
ols_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'ols_table_data.csv'))
# Write the results to a LaTeX table.
with open(project_paths_join('OUT_TABLES', 'estOLS_L3_new.tex'), 'w') as tex_file:

    # Top of table.
    tex_file.write('\\begin{tabular}{lc|ccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=50,000$', val3='$N=100,000$',val4='$N=1,000,000$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    # for i, var in enumerate(row_order):
    for row_idx, row in ols_data.iterrows():
        #print(row)
        # print(row_idx)
        # print(row['Row'])
        if row_idx<=1:
            row['true']
            tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=ols_data['beta_hat_50000'][row_idx], val2=ols_data['beta_hat_100000'][row_idx],val3=ols_data['beta_hat_1000000'][row_idx], str1='',str2='',str3='',se1=ols_data['sd_50000'][row_idx],se2=ols_data['sd_100000'][row_idx],se3=ols_data['sd_1000000'][row_idx]))
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{11cm}}{\\footnotesize{\\textit{Notes: Estimation results for the exponential discount factor $\\delta$ and the present-bias parameter $\\beta$ of the sophisticated agent for different simulated sample sizes (in the columns). The parameters are estimated using our OLS estimator.  Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')


md_naive_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'md_table_data_naive.csv'))
# Write the results to a LaTeX table.
with open(project_paths_join('OUT_TABLES', 'estMD_naive_new.tex'), 'w') as tex_file:

    # Top of table.
    tex_file.write('\\begin{tabular}{lc|ccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=50,000$', val3='$N=100,000$',val4='$N=1,000,000$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    # for i, var in enumerate(row_order):
    for row_idx, row in md_naive_data.iterrows():
        #print(row)
        # print(row_idx)
        # print(row['Row'])
        row['true']
        tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=md_naive_data['beta_hat_50000'][row_idx], val2=md_naive_data['beta_hat_100000'][row_idx],val3=md_naive_data['beta_hat_1000000'][row_idx], str1='',str2='',str3='',se1=md_naive_data['sd_50000'][row_idx],se2=md_naive_data['sd_100000'][row_idx],se3=md_naive_data['sd_1000000'][row_idx]))
        
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{11cm}}{\\footnotesize{\\textit{Notes: Estimation results for the exponential discount factor $\\delta$ and the present-bias parameter $\\beta$ of the naive agent for different simulated sample sizes (in the columns). The parameters are estimated using our MD estimator. Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')


mle_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'mle_table_data.csv'))
# Write the results to a LaTeX table.
with open(project_paths_join('OUT_TABLES', 'estMLE_new.tex'), 'w') as tex_file:

    # Top of table.
    tex_file.write('\\begin{tabular}{lc|ccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=5,000$', val3='$N=10,000$',val4='$N=20,000$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    # for i, var in enumerate(row_order):
    for row_idx, row in mle_data.iterrows():
        #print(row)
        # print(row_idx)
        # print(row['Row'])
        row['true']
        tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=mle_data['beta_hat_5000'][row_idx], val2=mle_data['beta_hat_10000'][row_idx],val3=mle_data['beta_hat_20000'][row_idx], str1='',str2='',str3='',se1=mle_data['sd_5000'][row_idx],se2=mle_data['sd_10000'][row_idx],se3=mle_data['sd_20000'][row_idx]))
        
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{11cm}}{\\footnotesize{\\textit{Notes: Estimation results for the exponential discount factor $\\delta$ and the present-bias parameter $\\beta$ and the flow utility parameters ($\\theta$) of the sophisticated agent for different simulated sample sizes (in the columns). The parameters are estimated using MLE.  Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')


mle_naive_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'mle_naive_table_data.csv'))
# Write the results to a LaTeX table.
with open(project_paths_join('OUT_TABLES', 'estMLE_naive_new.tex'), 'w') as tex_file:

    # Top of table.
    tex_file.write('\\begin{tabular}{lc|ccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=5,000$', val3='$N=10,000$',val4='$N=20,000$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    # for i, var in enumerate(row_order):
    for row_idx, row in mle_naive_data.iterrows():
        #print(row)
        # print(row_idx)
        # print(row['Row'])
        row['true']
        tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=mle_naive_data['beta_hat_5000'][row_idx], val2=mle_naive_data['beta_hat_10000'][row_idx],val3=mle_naive_data['beta_hat_20000'][row_idx], str1='',str2='',str3='',se1=mle_naive_data['sd_5000'][row_idx],se2=mle_naive_data['sd_10000'][row_idx],se3=mle_naive_data['sd_20000'][row_idx]))
        
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{11cm}}{\\footnotesize{\\textit{Notes: Estimation results for the exponential discount factor $\\delta$ and the present-bias parameter $\\beta$ and the flow utility parameters ($\\theta$) of the naive agent for different simulated sample sizes (in the columns). The parameters are estimated using MLE. Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')


mle_norm_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'mle_norm_table_data.csv'))

table_row_var = """{var} & {val0:.2f} & {val1:.4f}{str1} & {val2: .4f}{str2} & {val3: .4f}{str3} \\tabularnewline\n & & ({se1:.1e}) & ({se2:.1e}) & ({se3:.1e}) \\tabularnewline\n"""


# Write the results to a LaTeX table.
with open(project_paths_join('OUT_TABLES', 'estMLE_norm_new.tex'), 'w') as tex_file:

    # Top of table.
    tex_file.write('\\begin{tabular}{lc|ccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=5,000$', val3='$N=10,000$',val4='$N=20,000$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    # for i, var in enumerate(row_order):
    for row_idx, row in mle_norm_data.iterrows():
        #print(row)
        # print(row_idx)
        # print(row['Row'])
        row['true']
        tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=mle_norm_data['beta_hat_5000'][row_idx], val2=mle_norm_data['beta_hat_10000'][row_idx],val3=mle_norm_data['beta_hat_20000'][row_idx], str1='',str2='',str3='',se1=mle_norm_data['sd_5000'][row_idx],se2=mle_norm_data['sd_10000'][row_idx],se3=mle_norm_data['sd_20000'][row_idx]))
        
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{11cm}}{\\footnotesize{\\textit{Notes: Estimation results for the exponential discount factor $\\delta$ and the present-bias parameter $\\beta$ and the flow utility parameters ($\\theta$) of the sophisticated agent for different simulated sample sizes (in the columns). The parameters are estimated using MLE imposing that the flow utility of \\emph{not adopting} is normalized to zero. Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')

mle_normTerm_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'mle_normTerm_table_data.csv'))
# Write the results to a LaTeX table.
with open(project_paths_join('OUT_TABLES', 'estMLE_normTerm_new.tex'), 'w') as tex_file:

    # Top of table.
    tex_file.write('\\begin{tabular}{lc|ccc}\n\\toprule\n')
    # Header row with column names.
    tex_file.write(table_row_header.format(val1='True value', val2='$N=5,000$', val3='$N=10,000$',val4='$N=20,000$'))
    tex_file.write('\\midrule')
    # Write coefficients to table.
    # for i, var in enumerate(row_order):
    for row_idx, row in mle_normTerm_data.iterrows():
        #print(row)
        # print(row_idx)
        # print(row['Row'])
        row['true']
        tex_file.write(table_row_var.format(var=row['Row'], val0=row['true'], val1=mle_normTerm_data['beta_hat_5000'][row_idx], val2=mle_normTerm_data['beta_hat_10000'][row_idx],val3=mle_normTerm_data['beta_hat_20000'][row_idx], str1='',str2='',str3='',se1=mle_normTerm_data['sd_5000'][row_idx],se2=mle_normTerm_data['sd_10000'][row_idx],se3=mle_normTerm_data['sd_20000'][row_idx]))
        
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{11cm}}{\\footnotesize{\\textit{Notes: Estimation results for the exponential discount factor $\\delta$ and the present-bias parameter $\\beta$ and the flow utility parameters ($\\theta$) of the sophisticated agent for different simulated sample sizes (in the columns). The parameters are estimated using MLE imposing that the flow utility of \\emph{adopting} is normalized to zero. Standard deviations across simulations in parentheses.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')


