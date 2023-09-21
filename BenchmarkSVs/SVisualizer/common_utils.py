import os
import pandas as pd

### Reading files
def read_and_postprocess(path, postprocessor):
    # Check file path exists (was computed in WDL) and if not toggled off read and postprocess
    if os.path.isfile(path):
        df = pd.read_csv(path, sep='\t')
        df = postprocessor(df)
        return df
    else:
        return pd.DataFrame()

### General Postprocessing Functions
# Add BBEND stats
def add_bbend_stats(df):
    interval_names = sorted(list(set([c.split('-')[0] for c in df.columns if '-count' in c])))
    for i in interval_names:
        df[f'{i}-BBEND-overlap'] = (df[f'{i}-LBEND-overlap'] + df[f'{i}-RBEND-overlap'])/2
        df[f'{i}-BBEND-count'] = df[f'{i}-LBEND-count'] + df[f'{i}-RBEND-count']

# Convert '.' FILTER to 'PASS'
def convert_missing_to_pass_filter(df, MAKE_MISSING_PASS_FILTER):
    if 'FILTER' in df.columns:
        df['FILTER'] = df['FILTER'].replace('.', 'PASS')

# Sort SVTYPE with internal hierarchy
def sort_svtypes(svtypes):
    sorted_types = []
    if 'ALL' in svtypes:
        sorted_types += ['ALL']
        svtypes.remove('ALL')

    if 'DEL' in svtypes:
        sorted_types += ['DEL']
        svtypes.remove('DEL')

    if 'INS' in svtypes:
        sorted_types += ['INS']
        svtypes.remove('INS')

    if 'DUP' in svtypes:
        sorted_types += ['DUP']
        svtypes.remove('DUP')

    sorted_types += svtypes
    return sorted_types

# Sort AF and SVLEN bins
def sort_bins(bins, sorter):
    sorted_bins = ['ALL']
    bins.remove('ALL')

    other_bins = sorted(bins, key=sorter)
    sorted_bins += other_bins

    return sorted_bins

def get_size_pct(p):
    lt = True if p[0] == '<' else False
    p = p.split('-')[0].strip('>')
    p = p.strip('<').strip('%').strip()   # Clear whitespace for '< ' case
    p = int(p)

    if lt:
        # Diambiguate between e.g. '< 1%' and '1-10%'
        p -= 1
    return p

def sort_af_bins(bins):
    return sort_bins(bins, sorter=get_size_pct)

def get_size_strnum(s):
    s = s.split('-')[0].strip('>')
    lt = True if s[0] == '<' else False
    s = s.strip('<')
    
    final_num = 0
    if s[-1] == 'k':
        final_num = float(s.strip('k')) * 1000
    elif s[-1] == 'M':
        final_num = float(s.strip('M')) * 1_000_000
    else:
        final_num = float(s)

    if lt:
        # Disambiguate between e.g. '<100' and '100'
        final_num -= 1
    return final_num

def sort_svlen_bins(bins):
    return sort_bins(bins, sorter=get_size_strnum)

