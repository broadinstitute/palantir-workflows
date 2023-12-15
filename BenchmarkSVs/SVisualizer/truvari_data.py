import pandas as pd

from common_utils import read_and_postprocess, convert_missing_to_pass_filter, add_bbend_stats
from user_config import MAKE_MISSING_PASS_FILTER, TRUVARI_DUP_TO_INS


# Internal use only; restrict SVTYPEs globally in Truvari plots
RESTRICT_TYPES = ['ALL', 'DEL', 'INS']

def update_dup_to_ins(df):
    count_cols = ['tp_base_count', 'tp_base-BEND_count', 'tp_base-gt_count', 'tp_comp_count', 'tp_comp-BEND_count', 'tp_comp-gt_count',
                  'fp_count', 'fp-BEND_count', 'fn_count', 'fn-BEND_count']

    full_merge_columns = ['Base_Sample_Name', 'Comp_Sample_Name', 'Experiment', 'SVLEN_Bin', 'FILTER', 'Interval', 'Pct_Overlap']
    if 'Terra_workflow_id' in df.columns:
        full_merge_columns += ['Terra_workflow_id']
    exclude_cols = ['SVTYPE', 'Precision', 'Precision-BEND', 'Recall', 'Recall-BEND', 'F1_Score', 'F1_Score-BEND', 'GT_concordance']

    ins_df = df[df['SVTYPE'] == 'INS'].reset_index(drop=True).drop(columns=exclude_cols)
    dup_df = df[df['SVTYPE'] == 'DUP'].reset_index(drop=True).drop(columns=exclude_cols)
    new_ins = ins_df.merge(dup_df, on=full_merge_columns, how='outer')
    
    for col in count_cols:
        new_ins[col] = new_ins[f'{col}_x'].fillna(0) + new_ins[f'{col}_y'].fillna(0)
        new_ins = new_ins.drop(columns=[f'{col}_x', f'{col}_y'])
    
    new_ins['SVTYPE'] = 'INS'
    df = pd.concat([df[~df['SVTYPE'].isin(['INS', 'DUP'])], new_ins])

    # Regenerate Precision, Recall, F1 stats
    df['Precision'] = df['tp_comp_count'] / (df['tp_comp_count'] + df['fp_count'])
    df['Precision-BEND'] = df['tp_comp-BEND_count'] / (df['tp_comp-BEND_count'] + df['fp-BEND_count'])
    df['Recall'] = df['tp_base_count'] / (df['tp_base_count'] + df['fn_count'])
    df['Recall-BEND'] = df['tp_base-BEND_count'] / (df['tp_base-BEND_count'] + df['fn-BEND_count'])
    df['F1_Score'] = 2 * df['Precision'] * df['Recall'] / (df['Precision'] + df['Recall'])
    df['F1_Score-BEND'] = 2 * df['Precision-BEND'] * df['Recall-BEND'] / (df['Precision-BEND'] + df['Recall-BEND'])

    df['GT_concordance'] = df['tp_comp-gt_count'] / (df['tp_comp_count'])
    return df

def postprocess_truvari_bench(df):
    df['Experiment'] = df['Experiment'].astype(str)
    convert_missing_to_pass_filter(df, MAKE_MISSING_PASS_FILTER)
    if TRUVARI_DUP_TO_INS:
        df = update_dup_to_ins(df)

    if RESTRICT_TYPES is not None:
        df = df[df['SVTYPE'].isin(RESTRICT_TYPES)]

    return df

def postprocess_truvari_closest(df):
    df['LLEN'] = df['LLEN'].replace('.', 0).astype(int)
    df['RLEN'] = df['RLEN'].replace('.', 0).astype(int)
    df['SIZE_RATIO'] = df[['LLEN', 'RLEN']].min(axis=1) / df[['LLEN', 'RLEN']].max(axis=1)
    
    convert_missing_to_pass_filter(df, MAKE_MISSING_PASS_FILTER)

    k = 3  # Toggle this to compare other k-closest; should match WDL input
    df['kth_closest'] = df.index % k + 1

    return df
