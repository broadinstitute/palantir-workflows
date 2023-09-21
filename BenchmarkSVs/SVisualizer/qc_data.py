import pandas as pd
from user_config import MAKE_MISSING_PASS_FILTER
from common_utils import add_bbend_stats, convert_missing_to_pass_filter, sort_af_bins, sort_svtypes, sort_svlen_bins


# Postprocessing
def postprocess_qc_df(qc_df):
    # Separate base vs comp experiment groups
    qc_df['Experiment_Suffix'] = qc_df['Experiment'].apply(lambda x: x.split('-')[-1])

    # Make categorical columns in proper order
    af_bins = sort_af_bins(list(qc_df['AF_Bin'].unique()))
    qc_df['AF_Bin'] = pd.Categorical(qc_df['AF_Bin'], ordered=True, categories=af_bins)

    svlen_bins = sort_svlen_bins(list(qc_df['SVLEN_Bin'].unique()))
    qc_df['SVLEN_Bin'] = pd.Categorical(qc_df['SVLEN_Bin'], ordered=True, categories=svlen_bins)

    svtypes = sort_svlen_bins(list(qc_df['SVTYPE'].unique()))
    qc_df['SVTYPE'] = pd.Categorical(qc_df['SVTYPE'], ordered=True, categories=svtypes)

    add_bbend_stats(qc_df)
    convert_missing_to_pass_filter(qc_df, MAKE_MISSING_PASS_FILTER)

    return qc_df
