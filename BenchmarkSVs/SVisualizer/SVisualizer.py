#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import quickboard.base as qbb
from quickboard.app import start_app, deploy_app, app


# In[ ]:





# # User Config
# 
# Use this section to set variables that change the output and settings of the dashboard. Anything after this header's section should be untouched, and modifying the following code is advanced usage.

# In[ ]:


# (Optional) Covariate column names for benchmarking plots
# Should be numeric to plot along x-axis, e.g. average coverage
COVARIATE_X = 'Experiment'

# Option to skip second HWE plot for variant in Comp VCFs
# Useful to skip for single sample vs panel (Base) comparisons
SKIP_COMP_HWE = True

# Coerice missing FILTER to PASS
# Ensure '.' is interpreted as PASS FILTER in case base/comp have different conventions
# (Though this should be done before running through Wittyer...)
MAKE_MISSING_PASS_FILTER = True

# Use to fix explicit order for experiment groups
EXPERIMENT_ORDER = ['GATK-SV_HGSVC2'] + [f'Sniffles{i}x_vs_HGSVC2' for i in [2, 5, 8]]

# Use if running Truvari w/ dup-to-ins on
TRUVARI_DUP_TO_INS = True


# In[ ]:


# Toggle which categories of data to plot in Dashboard
INCLUDE_QC = False
INCLUDE_WITTYER = False    # Based on legacy code for generating Wittyer stats in WDL
INCLUDE_TRUVARI = True


# In[ ]:





# ## Write New User Config
# 
# DO NOT MODIFY CODE BELOW THIS SECTION
# 
# Code here writes above user config settings to appropriate file for tab generation

# In[ ]:


with open('user_config.py', 'w') as config:
    config.write(f'COVARIATE_X = "{COVARIATE_X}"\n')
    config.write(f'SKIP_COMP_HWE = {SKIP_COMP_HWE}\n')
    config.write(f'MAKE_MISSING_PASS_FILTER = {MAKE_MISSING_PASS_FILTER}\n')
    config.write(f'EXPERIMENT_ORDER = {EXPERIMENT_ORDER}\n')
    config.write(f'TRUVARI_DUP_TO_INS = {TRUVARI_DUP_TO_INS}\n\n')

    config.write(f'INCLUDE_QC = {INCLUDE_QC}\n')
    config.write(f'INCLUDE_WITTYER = {INCLUDE_WITTYER}\n')
    config.write(f'INCLUDE_TRUVARI = {INCLUDE_TRUVARI}\n')


# In[ ]:





# In[ ]:





# # Import Tabs

# In[ ]:


if INCLUDE_QC:
    from qc_tabs import counts_tab, histogram_tab, hwe_tab


# In[ ]:


if INCLUDE_WITTYER:
    from wittyer_tabs import basic_wittyer_tab, adv_wittyer_tab


# In[ ]:


if INCLUDE_TRUVARI:
    from truvari_tabs import truvari_bench_tab, truvari_errors_tab


# In[ ]:





# # Main Board

# In[ ]:


tab_list = []
if INCLUDE_QC:
    tab_list += [counts_tab, histogram_tab, hwe_tab]
if INCLUDE_WITTYER:
    tab_list += [basic_wittyer_tab, adv_wittyer_tab]
if INCLUDE_TRUVARI:
    tab_list += [truvari_bench_tab, truvari_errors_tab]

board = qbb.Quickboard(
    tab_list=tab_list
)


# In[ ]:


start_app(board, debug=False, app_title='SVisualizer', mode='external', port=8050)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




