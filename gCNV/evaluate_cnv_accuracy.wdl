version 1.0

workflow evaluate_cnv_accuracy {
    input {
        String version
        Array[File] cnv_vcfs
        Array[File] merged_vcfs
        Array[File] truth_beds
        File sample_id_mappings
        File sample_material_mappings
        File targets_interval_list
        Int truth_ref_panel_counts
    }

    call compute_accuracy {
        input:
            cnv_vcfs = cnv_vcfs,
            merged_vcfs = merged_vcfs,
            truth_beds = truth_beds,
            sample_id_mapping = sample_id_mappings,
            sample_material_mapping = sample_material_mappings,
            targets_interval_list = targets_interval_list,
            truth_ref_panel_count = truth_ref_panel_counts,
            version = version
    }

    output {
        File ppv = compute_accuracy.ppv
        File recall = compute_accuracy.recall
    }
}


task compute_accuracy {
    input {
        Array[File] cnv_vcfs
        Array[File] merged_vcfs
        Array[File] truth_beds
        File sample_id_mapping
        File sample_material_mapping
        File targets_interval_list
        Int truth_ref_panel_count
        String version
        Int mem_gb = 8
    }

    command <<<

        python << "EOF"
        import pandas as pd
        import numpy as np

        targets = pd.read_csv("~{targets_interval_list}", sep="\t", comment="@", names = ["contig","start","end","dummy1","dummy2"])
        targets['contig_idx'] = targets.groupby('contig').cumcount()
        targets = targets.set_index(targets.contig + "_" + targets.contig_idx.astype(str))

        df_truth = pd.concat([pd.read_csv(truth_bed, sep="\t") for truth_bed in ["~{sep='","' truth_beds}"]])
        df_truth_subset = df_truth.loc[:,["#chrom","start","end","name","samples","svtype","FILTER"]]
        df_truth_subset.loc[:,"name"]=df_truth_subset.apply(lambda x: x["name"].split(".")[0], axis=1)
        df_truth_subset["ac"]=df_truth_subset.apply(lambda x: len(x["samples"].split(",")), axis=1)
        sv_ref_panel_count = ~{truth_ref_panel_count}
        df_truth_subset["af"] = df_truth_subset["ac"]/(sv_ref_panel_count+1)
        df_truth_subset["case_sample_event"]=df_truth_subset.apply(lambda x: x["name"] in x["samples"], axis=1)
        df_truth_subset=df_truth_subset[df_truth_subset["case_sample_event"]]

        sample_id_mapping = pd.read_csv("~{sample_id_mapping}", sep="\t", dtype=str)
        sample_material_mapping = pd.read_csv("~{sample_material_mapping}", sep="\t", dtype=str)

        df_truth_subset = df_truth_subset.merge(sample_id_mapping, left_on="name", right_on="truth_sample_id", how="left")
        df_truth_subset = df_truth_subset.merge(sample_material_mapping, left_on="bge_sample_id", right_on="bge_sample_id", how="left")
        df_truth_subset=df_truth_subset.rename({"#chrom":"contig"}, axis=1)

        def add_exon_idxs(df, exons):
            contigs = set(df.contig)
            for contig in contigs:
                df.loc[df.contig==contig,"start_exon_idx"]=np.searchsorted(exons.loc[exons.contig==contig].end,
                                                                        df.loc[df.contig==contig].start,"left")
                df.loc[df.contig==contig,"end_exon_idx"]=np.searchsorted(exons.loc[exons.contig==contig].start,
                                                                        df.loc[df.contig==contig].end,"right")

        def get_exon_expanded_events(df, exons):
            add_exon_idxs(df, exons)
            df = df.loc[df.start_exon_idx != df.end_exon_idx].reset_index()
            df_expanded = df.loc[df.index.repeat(df.end_exon_idx-df.start_exon_idx)]
            df_expanded['exon_idx'] = df_expanded.groupby(df_expanded.index).cumcount() + df_expanded.start_exon_idx
            df_expanded = df_expanded.set_index(df_expanded.contig + "_" + df_expanded.exon_idx.astype(int).astype(str))
            df_expanded = df_expanded.join(exons[['start','end']], rsuffix='_exon')
            df_expanded['event_exon_start']=np.maximum(df_expanded.start, df_expanded.start_exon)
            df_expanded['event_exon_end']=np.minimum(df_expanded.end, df_expanded.end_exon)
            return df_expanded

        df_truth_subset = df_truth_subset.loc[(df_truth_subset.svtype=='DEL') | (df_truth_subset.svtype=='DUP')]
        df_truth_exon_expanded = get_exon_expanded_events(df_truth_subset, targets)
        df_truth_exon_expanded = df_truth_exon_expanded.set_index(df_truth_exon_expanded.index + 
                                "_" + 
                                df_truth_exon_expanded.bge_sample_id + 
                                "_" +
                                df_truth_exon_expanded.svtype)

        def standardize_vcf_df(df):
            df["svtype"] = df.ALT.str.replace("<","").str.replace(">","")
            df = df.query("svtype in ['DEL', 'DUP']").copy()
            df.index=list(range(len(df)))
            df['end'] = df['INFO'].str.replace('END=','').apply(lambda x: x.split(';')[0]).astype(int) #split('=').apply(lambda x: x[1])
            
            
            for num,f in enumerate(df.loc[0,'FORMAT'].split(':')):
                df[f] = df['SAMPLE'].str.split(':').apply(lambda x: x[num])
            df = df.astype({'NP':int,'CN':int,'QA':int,'QS':int, 'QUAL':float})
            return df

        def read_bge_sample_to_df(vcf):
            df = pd.read_csv(vcf, sep='\t',comment='#',compression='gzip',
                                names=['contig','start','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE'])
            df = standardize_vcf_df(df)
            df['sample_name']=vcf.split('/')[-1].split('.')[0]
            return df
        


        df_bge = pd.concat([read_bge_sample_to_df(vcf) for vcf in ["~{sep='","' cnv_vcfs}"]])
        df_bge_merged = pd.concat([read_bge_sample_to_df(vcf) for vcf in ["~{sep='","' merged_vcfs}"]])

        passing_events = df_bge.query("FILTER == 'PASS'")[["contig","start","end","ID","REF","ALT"]].copy()
        passing_events = passing_events.sort_values(list(passing_events.columns)).reset_index(drop=True)
        passing_events_merged = df_bge_merged[["contig","start","end","ID","REF","ALT"]].copy()
        passing_events_merged = passing_events_merged.sort_values(list(passing_events_merged.columns)).reset_index(drop=True)
        if not passing_events.equals(passing_events_merged):
            raise ValueError("The passing events in the individual VCFs do not match the passing events in the merged VCF.")
    
        df_bge = df_bge.merge(sample_material_mapping, left_on="sample_name", right_on="bge_sample_id", how="left")
        df_bge_expanded = get_exon_expanded_events(df_bge, targets)
        df_bge_expanded = df_bge_expanded.set_index(df_bge_expanded.index + "_" + df_bge_expanded.sample_name + "_" + df_bge_expanded.svtype)

        bge_truth_joined = df_bge_expanded.join(df_truth_exon_expanded, lsuffix='_bge', rsuffix='_truth', how='outer')

        bge_truth_joined['overlapping_bge_truth_exon_start'] = np.maximum(bge_truth_joined.event_exon_start_truth, bge_truth_joined.event_exon_start_bge)
        bge_truth_joined['overlapping_bge_truth_exon_end'] = np.minimum(bge_truth_joined.event_exon_end_truth, bge_truth_joined.event_exon_end_bge)

        bge_truth_joined['overlapping_bge_truth_exon_length'] = np.maximum(0, 
                                                                  bge_truth_joined.overlapping_bge_truth_exon_end -
                                                                  bge_truth_joined.overlapping_bge_truth_exon_start +1).fillna(0)
        
        df_bge_expanded['event_exon_length'] = df_bge_expanded['event_exon_end']-df_bge_expanded['event_exon_start'] + 1
        bge_length_df = (
            df_bge_expanded.rename({"contig":"contig_bge","start":"start_bge","svtype":"svtype_bge",
                                    "input_material_type":"input_material_type_bge","FILTER":"FILTER_bge"}, axis=1).
            groupby(['contig_bge','start_bge','svtype_bge','sample_name', 'input_material_type_bge', 'QUAL',
                'FILTER_bge','CN','NP']).
            agg(bge_length = pd.NamedAgg('event_exon_length',sum),
                n_exons = pd.NamedAgg('event_exon_length','count'))

        )
        bge_label_df = (
            bge_truth_joined.assign(overlapping_bge_truth_exon_length_pass = np.where(bge_truth_joined.FILTER_truth=="PASS",
                                                                 bge_truth_joined.overlapping_bge_truth_exon_length,
                                                                0)).
            groupby(['contig_bge','start_bge','svtype_bge','sample_name', 'input_material_type_bge', 'QUAL',
                'FILTER_bge','CN','NP']).
            agg(overlapping_truth_length_any=pd.NamedAgg('overlapping_bge_truth_exon_length',sum),
                overlapping_truth_length_pass=pd.NamedAgg('overlapping_bge_truth_exon_length_pass',sum)).
            join(bge_length_df)
        )
        bge_label_df['overlap_truth_frac_any'] = bge_label_df.overlapping_truth_length_any/bge_label_df.bge_length
        bge_label_df['overlap_truth_frac_pass'] = bge_label_df.overlapping_truth_length_pass/bge_label_df.bge_length

        ppv_dfs = [bge_label_df.query(f'FILTER_bge=="PASS" and n_exons >= {min_exon_count}').
               groupby(["svtype_bge","input_material_type_bge"]).
               agg(TP_frac = pd.NamedAgg('overlap_truth_frac_pass','sum'),
                  FP_frac = pd.NamedAgg('overlap_truth_frac_any',lambda x: (1-x).sum()))
                   for min_exon_count in range(1,10)]
        
        for ppv_df, min_exon_count in zip(ppv_dfs, list(range(1,10))):
            ppv_df['min_exon_count']=min_exon_count
        ppv_df = pd.concat(ppv_dfs)
        ppv_df['ppv']=ppv_df['TP_frac']/(ppv_df['TP_frac']+ppv_df['FP_frac'])
        ppv_df.reset_index().to_csv("ppv_~{version}.csv", sep="\t", index=False)

        df_truth_exon_expanded['event_exon_length'] = df_truth_exon_expanded['event_exon_end']-df_truth_exon_expanded['event_exon_start'] + 1
        truth_length_df = (
            df_truth_exon_expanded.rename({"contig":"contig_truth","start":"start_truth","svtype":"svtype_truth",
                                    "input_material_type":"input_material_type_truth","FILTER":"FILTER_truth",
                                    "bge_sample_id":"bge_sample_id_truth"}, axis=1).
            groupby(['contig_truth','start_truth','svtype_truth','bge_sample_id_truth', 'input_material_type_truth',
                'FILTER_truth','af']).
            agg(truth_length = pd.NamedAgg('event_exon_length',sum),
                n_exons = pd.NamedAgg('event_exon_length','count'))
        )

        truth_label_df = (
            bge_truth_joined.assign(overlapping_bge_truth_exon_length_pass = np.where(bge_truth_joined.FILTER_bge=="PASS",
                                                                 bge_truth_joined.overlapping_bge_truth_exon_length,
                                                                0)).
            groupby(['contig_truth','start_truth','svtype_truth','bge_sample_id_truth', 'input_material_type_truth',
                'FILTER_truth','af']).
            agg(overlapping_bge_length_pass=pd.NamedAgg('overlapping_bge_truth_exon_length_pass',sum)).
            join(truth_length_df)
        )
        truth_label_df['overlap_bge_frac_pass'] = truth_label_df.overlapping_bge_length_pass/truth_label_df.truth_length

        recall_dfs = [truth_label_df.query(f'FILTER_truth=="PASS" and n_exons >= {min_exon_count} and '
                                                     'af < 0.01 and contig_truth != "chrY"').groupby(["svtype_truth","input_material_type_truth"]).
                agg(TP_frac=pd.NamedAgg('overlap_bge_frac_pass', 'sum'),
                   FN_frac=pd.NamedAgg('overlap_bge_frac_pass',lambda x: (1-x).sum())) for min_exon_count in range(1,10)]
        for recall_df, min_exon_count in zip(recall_dfs, list(range(1,10))):
            recall_df['min_exon_count']=min_exon_count
        recall_df=pd.concat(recall_dfs)
        recall_df['recall']=recall_df.TP_frac/(recall_df.TP_frac + recall_df.FN_frac)
        recall_df.reset_index().to_csv("recall_~{version}.csv", sep="\t", index=False)
        EOF

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim@sha256:4c880fc4ca079f294dae4aef27c61fe58398684b2445caec5a484650f07a2616"
        preemptible: 3
        cpu: 2
        disks: "local-disk 50 HDD"
        memory: mem_gb + " GB"
    }

    output {
        File ppv = "ppv_~{version}.csv"
        File recall = "recall_~{version}.csv"

    }
}