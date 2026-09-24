version 1.0

import "single_sample_cnv_germline_case_filter_workflow.wdl" as cnv_case_and_filter
workflow CNVCallingAndMergeForFabric {
    input {
        File normal_bam
        File normal_bai
        File short_variant_vcf

        File contig_ploidy_model_tar
        File preprocessed_intervals
        File gcnv_model_tar
        Array[File]+ gcnv_panel_genotyped_segments
        Array[File]+ gcnv_panel_copy_ratios
        Array[File]+ gcnv_panel_read_counts

        Float overlap_thresh = 0.5

        String gatk_docker

        Int maximum_number_events_per_sample
        Int maximum_number_pass_events_per_sample
        Int ref_copy_number_autosomal_contigs

        File ref_fasta
        File ref_fasta_fai
        File ref_fasta_dict

        Array[String] allosomal_contigs
    }

    call cnv_case_and_filter.SingleSampleGCNVAndFilterVCFs {
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            contig_ploidy_model_tar = contig_ploidy_model_tar,
            preprocessed_intervals = preprocessed_intervals,
            gcnv_model_tar = gcnv_model_tar,
            pon_genotyped_segments_vcfs = gcnv_panel_genotyped_segments,
            gatk_docker = gatk_docker,
            maximum_number_events_per_sample = maximum_number_events_per_sample,
            maximum_number_pass_events_per_sample = maximum_number_pass_events_per_sample,
            ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            allosomal_contigs = allosomal_contigs,
            overlap_thresh = overlap_thresh
    }

    call ReformatAndMergeForFabric {
        input:
            cnv_vcf = SingleSampleGCNVAndFilterVCFs.filtered_vcf,
            short_variant_vcf = short_variant_vcf,
            gatk_docker = gatk_docker
    }

    call GCNVVisualzation {
        input:
            filtered_vcf = SingleSampleGCNVAndFilterVCFs.filtered_vcf,
            case_copy_ratios = SingleSampleGCNVAndFilterVCFs.denoised_copy_ratios,
            case_read_counts = SingleSampleGCNVAndFilterVCFs.read_counts,
            panel_copy_ratios = gcnv_panel_copy_ratios,
            panel_read_counts = gcnv_panel_read_counts,
            interval_lists = [SingleSampleGCNVAndFilterVCFs.interval_list]

    }

    output {
        File filtered_cnv_genotyped_segments_vcf = SingleSampleGCNVAndFilterVCFs.filtered_vcf
        File filtered_cnv_genotyped_segments_vcf_index = SingleSampleGCNVAndFilterVCFs.filtered_vcf_index
        File filtered_cnv_genotyped_segments_vcf_md5sum = SingleSampleGCNVAndFilterVCFs.filtered_vcf_md5sum

        File merged_vcf = ReformatAndMergeForFabric.merged_vcf
        File merged_vcf_index = ReformatAndMergeForFabric.merged_vcf_index
        File merged_vcf_md5sum = ReformatAndMergeForFabric.merged_vcf_md5sum

        Boolean qc_passed = SingleSampleGCNVAndFilterVCFs.qc_passed
        File cnv_metrics = SingleSampleGCNVAndFilterVCFs.cnv_metrics
        File cnv_event_report = GCNVVisualzation.cnv_event_report

    }
}

#Fabric doesn't seem to like ./. genotypes on
task ReformatGCNVForFabric {
    input {
        File cnv_vcf
        Int disk_size_gb = 20
        Int mem_gb = 4
    }

    String output_basename = basename(cnv_vcf, ".filtered.genotyped-segments.vcf.gz")

    command <<<
        set -euo pipefail

        python << CODE
        from pysam import VariantFile

        with VariantFile("~{cnv_vcf}") as cnv_vcf_in:
            header_out = cnv_vcf_in.header
            header_out.info.add("CN", "A", "Integer", "Copy number associated with <CNV> alleles")
            with VariantFile("~{output_basename}.reformatted_for_fabric.vcf.gz",'w', header = header_out) as cnv_vcf_out:
                for rec in cnv_vcf_in.fetch():
                    if 'PASS' in rec.filter:
                        if rec.alts and rec.alts[0] == "<DUP>":
                            for rec_sample in rec.samples.values():
                                ploidy = len(rec_sample.alleles)
                                rec_sample.allele_indices = (None,)*(ploidy - 1) + (1,)
                        cnv_vcf_out.write(rec)
        CODE
    >>>

    runtime {
            docker: "us.gcr.io/broad-dsde-methods/pysam:v1.1"
            preemptible: 3
            cpu: 2
            disks: "local-disk " + disk_size_gb + " HDD"
            memory: mem_gb + " GB"
        }

    output {
        File reformatted_vcf = "~{output_basename}.reformatted_for_fabric.vcf.gz"
    }
}

task MergeVcfs {
    input {
        File cnv_vcf
        File short_variant_vcf

        String gatk_docker
        Int mem_gb=4
        Int disk_size_gb = 100
    }

    String output_basename = basename(short_variant_vcf, ".hard-filtered.vcf.gz")
    String output_cnv_basename = basename(cnv_vcf, ".reformatted_for_fabric.vcf.gz")
    command <<<
        set -euo pipefail

        if [ "~{output_basename}" != "~{output_cnv_basename}" ]; then
            echo "input vcf names do not agree"
            exit 1
        fi

        gatk --java-options "-Dsamjdk.create_md5=true" MergeVcfs -I ~{short_variant_vcf} -I ~{cnv_vcf} -O ~{output_basename}.merged.vcf.gz

        mv ~{output_basename}.merged.vcf.gz.md5 ~{output_basename}.merged.vcf.gz.md5sum
    >>>

    output {
        File merged_vcf = "~{output_basename}.merged.vcf.gz"
        File merged_vcf_index = "~{output_basename}.merged.vcf.gz.tbi"
        File merged_vcf_md5sum = "~{output_basename}.merged.vcf.gz.md5sum"
    }

     runtime {
        docker: gatk_docker
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 5
    }
}

task ReformatAndMergeForFabric {
    input {
            File cnv_vcf
            File short_variant_vcf


            String gatk_docker
            Int mem_gb=4
            Int disk_size_gb = 100
        }

    String output_basename = basename(cnv_vcf, ".filtered.genotyped-segments.vcf.gz")
    String output_short_variant_basename = basename(short_variant_vcf, ".hard-filtered.vcf.gz")

    command <<<
        set -euo pipefail

        if [ "~{output_basename}" != "~{output_short_variant_basename}" ]; then
            echo "input vcf names do not agree"
            exit 1
        fi

        python << CODE
        from pysam import VariantFile

        with VariantFile("~{cnv_vcf}") as cnv_vcf_in:
            header_out = cnv_vcf_in.header
            header_out.info.add("CN", "A", "Integer", "Copy number associated with <CNV> alleles")
            with VariantFile("~{output_basename}.reformatted_for_fabric.vcf.gz",'w', header = header_out) as cnv_vcf_out:
                for rec in cnv_vcf_in.fetch():
                    if 'PASS' in rec.filter:
                        if rec.alts and rec.alts[0] == "<DUP>":
                            for rec_sample in rec.samples.values():
                                ploidy = len(rec_sample.alleles)
                                rec_sample.allele_indices = (None,)*(ploidy - 1) + (1,)
                        cnv_vcf_out.write(rec)
        CODE

        gatk --java-options "-Dsamjdk.create_md5=true" MergeVcfs -I ~{short_variant_vcf} -I ~{output_basename}.reformatted_for_fabric.vcf.gz -O ~{output_basename}.merged.vcf.gz

        mv ~{output_basename}.merged.vcf.gz.md5 ~{output_basename}.merged.vcf.gz.md5sum
    >>>

    output {
        File merged_vcf = "~{output_basename}.merged.vcf.gz"
        File merged_vcf_index = "~{output_basename}.merged.vcf.gz.tbi"
        File merged_vcf_md5sum = "~{output_basename}.merged.vcf.gz.md5sum"
    }

    runtime {
        docker: gatk_docker
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 5
    }
}

task GCNVVisualzation {
    input {
        File filtered_vcf
        File case_copy_ratios
        File case_read_counts
        Array[File]+ panel_copy_ratios
        Array[File]+ panel_read_counts
        Array[File]+ interval_lists
        Int mem_gb=4
    }

    String output_prefix = basename(filtered_vcf, ".filtered.genotyped-segments.vcf.gz")

    # Produces a self-contained HTML report with plots for each passing CNV call and an Interval Browser tab in which
    # the user can enter any interval to plot it. The case and panel data for every interval are embedded in the HTML,
    # so no server or internet connection is needed to view it.
    command <<<
        set -euo pipefail

        cat << 'EOF' > gcnv_visualization.py
        """Build a self-contained interactive HTML report of gCNV calls.

        The report embeds the case and panel denoised copy ratios and adjusted read counts for every
        interval, so plots are drawn in the browser both for each passing CNV call and for any interval
        the user enters in the Interval Browser tab. No server or internet connection is needed to view it.
        """
        import argparse
        import base64
        import gzip
        import html
        import json
        import re
        import warnings

        import h5py
        import numpy as np
        import pandas as pd

        # Panel values are stored as uint8 round(value * PANEL_SCALE), capped at 254; 255 marks a missing value. Each
        # interval's panel values are sorted, which makes the embedded data several times smaller after compression; the
        # plots do not need to know which panel sample each value came from.
        PANEL_SCALE = 32
        PANEL_MISSING = 255


        def read_list(path):
            with open(path) as f:
                return [line.strip() for line in f if line.strip()]


        def read_header(path):
            header = []
            with open(path) as f:
                for line in f:
                    if not line.startswith("@"):
                        break
                    header.append(line.rstrip("\n"))
            return header


        def header_field(header, record_type, tag):
            values = []
            for line in header:
                if line.startswith(record_type):
                    match = re.search(tag + r":([^\t]+)", line)
                    if match:
                        values.append(match.group(1))
            return values


        def read_copy_ratios(path):
            header = read_header(path)
            sample_names = header_field(header, "@RG", "SM")
            table = pd.read_csv(path, sep="\t", skiprows=len(header), dtype={"CONTIG": str})
            return table, (sample_names[0] if sample_names else path), header_field(header, "@SQ", "SN")


        def to_str(value):
            return value.decode("utf-8") if isinstance(value, bytes) else str(value)


        def read_normalized_counts(path):
            """Return contig, start, end and counts normalized by the sample's mean count, as in the original R report."""
            with h5py.File(path, "r") as f:
                contig_names = np.array([to_str(c) for c in f["intervals/indexed_contig_names"][:]])
                index_start_end = f["intervals/transposed_index_start_end"][:]
                counts = f["counts/values"][:].ravel().astype(np.float64)
                sample_name = to_str(f["sample_metadata/sample_name"][0])
            if index_start_end.shape[0] != 3:
                index_start_end = index_start_end.T
            return (contig_names[index_start_end[0].astype(int)], index_start_end[1], index_start_end[2],
                    counts / counts.mean(), sample_name)


        def quantize(values):
            q = np.clip(np.round(values * PANEL_SCALE), 0, PANEL_MISSING - 1)
            return np.where(np.isfinite(values), q, PANEL_MISSING).astype(np.uint8)


        def parse_info(info):
            fields = {}
            for entry in info.split(";"):
                key, _, value = entry.partition("=")
                fields[key] = value
            return fields


        def to_number(value, kind=float):
            try:
                return kind(value)
            except (TypeError, ValueError):
                return None


        def read_calls(vcf_path):
            """Return the non-reference records of the genotyped-segments VCF."""
            opener = gzip.open if vcf_path.endswith(".gz") else open
            calls = []
            with opener(vcf_path, "rt") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    chrom, pos, record_id, _, alt, qual, filters, info = fields[:8]
                    if alt == ".":
                        continue
                    info = parse_info(info)
                    sample = dict(zip(fields[8].split(":"), fields[9].split(":"))) if len(fields) > 9 else {}
                    calls.append({
                        "contig": chrom,
                        "start": int(pos),
                        "end": to_number(info.get("END"), int) or int(pos),
                        "id": record_id,
                        "alt": alt.strip("<>"),
                        "qual": to_number(qual),
                        "filter": filters,
                        "pass": filters == "PASS",
                        "cn": sample.get("CN"),
                        "gt": sample.get("GT"),
                        "np": to_number(sample.get("NP"), int),
                        "panel_freq": to_number(info.get("PANEL_FREQ")),
                        "panel_count": to_number(info.get("PANEL_COUNT"), int),
                    })
            return calls


        def pack(array):
            data = gzip.compress(np.ascontiguousarray(array).tobytes(), compresslevel=6, mtime=0)
            return base64.b64encode(data).decode("ascii")


        def main():
            parser = argparse.ArgumentParser(description=__doc__)
            parser.add_argument("--vcf", required=True, help="filtered genotyped-segments VCF")
            parser.add_argument("--case-copy-ratios", required=True)
            parser.add_argument("--case-read-counts", required=True)
            parser.add_argument("--panel-copy-ratios-list", required=True, help="file listing panel denoised copy ratio TSVs")
            parser.add_argument("--panel-read-counts-list", required=True, help="file listing panel read count HDF5s")
            parser.add_argument("--interval-lists-list", required=True, help="file listing gCNV interval_list.tsv files")
            parser.add_argument("--title", required=True)
            parser.add_argument("--output", required=True)
            args = parser.parse_args()

            # Case copy ratios define the intervals that are plotted.
            case_cr, case_sample, sequence_order = read_copy_ratios(args.case_copy_ratios)
            contig_rank = {c: i for i, c in enumerate(sequence_order)}
            for contig in pd.unique(case_cr["CONTIG"]):
                contig_rank.setdefault(contig, len(contig_rank))
            case_cr = (case_cr.assign(rank=case_cr["CONTIG"].map(contig_rank))
                       .sort_values(["rank", "START"], kind="stable").reset_index(drop=True))
            n_intervals = len(case_cr)
            interval_index = pd.MultiIndex.from_arrays([case_cr["CONTIG"].to_numpy(), case_cr["START"].to_numpy(),
                                                        case_cr["END"].to_numpy()])

            def align(contig, start, end):
                return interval_index.get_indexer(pd.MultiIndex.from_arrays([np.asarray(contig), np.asarray(start),
                                                                             np.asarray(end)]))

            # Panel denoised copy ratios.
            panel_cr_paths = read_list(args.panel_copy_ratios_list)
            panel_cr = np.full((n_intervals, len(panel_cr_paths)), PANEL_MISSING, dtype=np.uint8)
            panel_cr_samples = []
            for j, path in enumerate(panel_cr_paths):
                table, sample_name, _ = read_copy_ratios(path)
                idx = align(table["CONTIG"], table["START"], table["END"])
                found = idx >= 0
                panel_cr[idx[found], j] = quantize(table["LINEAR_COPY_RATIO"].to_numpy()[found])
                panel_cr_samples.append(sample_name)

            # Read counts: normalize each sample by its mean count, then scale each interval by the panel mean so that
            # copy number 2 sits at 2.
            panel_rc_paths = read_list(args.panel_read_counts_list)
            panel_norm = np.full((len(panel_rc_paths), n_intervals), np.nan, dtype=np.float32)
            panel_rc_samples = []
            for j, path in enumerate(panel_rc_paths):
                contig, start, end, norm_counts, sample_name = read_normalized_counts(path)
                idx = align(contig, start, end)
                found = idx >= 0
                panel_norm[j, idx[found]] = norm_counts[found]
                panel_rc_samples.append(sample_name)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                panel_mean = np.nanmean(panel_norm, axis=0)
                panel_adjusted = 2 * panel_norm / panel_mean
            del panel_norm
            panel_adjusted[~np.isfinite(panel_adjusted)] = np.nan
            panel_adj = np.ascontiguousarray(quantize(panel_adjusted).T)
            del panel_adjusted

            contig, start, end, norm_counts, _ = read_normalized_counts(args.case_read_counts)
            idx = align(contig, start, end)
            found = idx >= 0
            case_adj = np.full(n_intervals, np.nan, dtype=np.float32)
            with np.errstate(divide="ignore", invalid="ignore"):
                case_adj[idx[found]] = 2 * norm_counts[found] / panel_mean[idx[found]]
            case_adj[~np.isfinite(case_adj)] = np.nan

            # GC content is only present when gCNV was run with annotated intervals.
            interval_tables = []
            for path in read_list(args.interval_lists_list):
                table = pd.read_csv(path, sep="\t", skiprows=len(read_header(path)), dtype={"CONTIG": str})
                if "GC_CONTENT" in table.columns:
                    interval_tables.append(table[["CONTIG", "START", "END", "GC_CONTENT"]])
            gc = None
            if interval_tables:
                intervals = pd.concat(interval_tables).drop_duplicates(["CONTIG", "START", "END"])
                idx = align(intervals["CONTIG"], intervals["START"], intervals["END"])
                found = idx >= 0
                gc = np.full(n_intervals, np.nan, dtype=np.float32)
                gc[idx[found]] = intervals["GC_CONTENT"].to_numpy()[found]

            starts = case_cr["START"].to_numpy().astype(np.int64)
            start_delta = np.diff(starts, prepend=0).astype(np.int32)
            lengths = (case_cr["END"].to_numpy() - starts).astype(np.int32)
            contigs = []
            for _, group in case_cr.groupby("rank", sort=True):
                contigs.append({"name": group["CONTIG"].iloc[0], "offset": int(group.index[0]), "count": len(group)})

            calls = read_calls(args.vcf)
            arrays = {
                "start_delta": start_delta,
                "length": lengths,
                "case_cr": case_cr["LINEAR_COPY_RATIO"].to_numpy().astype(np.float32),
                "case_adj": case_adj,
                "panel_cr": np.sort(panel_cr, axis=1),
                "panel_adj": np.sort(panel_adj, axis=1),
            }
            if gc is not None:
                arrays["gc"] = gc

            meta = {
                "title": args.title,
                "sample": case_sample,
                "n_intervals": n_intervals,
                "contigs": contigs,
                "panel_cr_samples": panel_cr_samples,
                "panel_adj_samples": panel_rc_samples,
                "panel_scale": PANEL_SCALE,
                "panel_missing": PANEL_MISSING,
                "has_gc": gc is not None,
                "calls": calls,
            }
            data_tags = "\n".join(f'<script type="application/octet-stream" id="gcnv-{name}">{pack(array)}</script>'
                                  for name, array in arrays.items())
            meta_json = json.dumps(meta, separators=(",", ":")).replace("</", "<\\/")
            report = (HTML_TEMPLATE.replace("__TITLE__", html.escape(args.title))
                      .replace("__META__", meta_json)
                      .replace("__DATA__", data_tags))
            with open(args.output, "w", encoding="utf-8") as f:
                f.write(report)
            n_pass = sum(c["pass"] for c in calls)
            print(f"Wrote {args.output}: {n_intervals} intervals, {len(panel_cr_samples)} panel copy ratio samples, "
                  f"{len(panel_rc_samples)} panel read count samples, {n_pass} passing calls, "
                  f"{len(calls) - n_pass} filtered calls")


        HTML_TEMPLATE = r"""<!DOCTYPE html>
        <html lang="en">
        <head>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1">
        <title>__TITLE__</title>
        <style>
          :root {
            --fg: #1f2937; --muted: #6b7280; --border: #e5e7eb; --soft: #f8f9fa;
            --accent: #1d4ed8; --err: #b91c1c;
          }
          * { box-sizing: border-box; }
          body { margin: 0; background: #fff; color: var(--fg);
                 font: 14px/1.45 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
          .wrap { max-width: 1500px; margin: 0 auto; padding: 0 24px; }
          header { padding-top: 18px; padding-bottom: 12px; }
          h1 { font-size: 22px; margin: 0 0 4px; word-break: break-all; }
          h2 { font-size: 17px; margin: 22px 0 8px; }
          h3 { font-size: 15px; margin: 0; }
          .muted { color: var(--muted); }
          nav { position: sticky; top: 0; z-index: 10; background: var(--soft); border-bottom: 2px solid var(--accent); }
          nav .wrap { display: flex; gap: 4px; padding-top: 6px; padding-bottom: 6px; }
          nav button { font: inherit; border: 0; background: none; padding: 6px 14px; border-radius: 16px; cursor: pointer; color: var(--fg); }
          nav button:hover { background: #e5e7eb; }
          nav button.active { background: var(--accent); color: #fff; }
          main { padding-top: 12px; padding-bottom: 60px; }
          .tab { display: none; }
          .tab.active { display: block; }
          .plot { position: relative; user-select: none; }
          .plot canvas { display: block; width: 100%; }
          .tip { position: absolute; display: none; pointer-events: none; z-index: 5; background: rgba(17,24,39,.93); color: #fff;
                 font-size: 12px; line-height: 1.4; padding: 6px 9px; border-radius: 4px; white-space: nowrap; }
          .guide { position: absolute; display: none; pointer-events: none; width: 1px; background: rgba(17,24,39,.35); }
          .sel { position: absolute; display: none; pointer-events: none; background: rgba(29,78,216,.12);
                 border-left: 1px solid var(--accent); border-right: 1px solid var(--accent); }
          .pair { display: grid; grid-template-columns: repeat(auto-fit, minmax(480px, 1fr)); gap: 8px 20px; }
          .pair h4 { margin: 6px 0 0; font-size: 13px; font-weight: 600; }
          .legend { display: flex; flex-wrap: wrap; gap: 6px 18px; font-size: 13px; color: var(--muted); margin: 6px 0; }
          .legend span { display: inline-flex; align-items: center; gap: 6px; }
          .sw { display: inline-block; width: 12px; height: 12px; border-radius: 50%; }
          .sw.rect { border-radius: 2px; }
          table { border-collapse: collapse; font-size: 13px; }
          th, td { padding: 4px 10px; border-bottom: 1px solid var(--border); text-align: left; white-space: nowrap; }
          th { background: var(--soft); font-weight: 600; }
          td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
          .table-scroll { overflow-x: auto; }
          a, .link { color: var(--accent); cursor: pointer; text-decoration: none; background: none; border: 0; font: inherit; padding: 0; }
          a:hover, .link:hover { text-decoration: underline; }
          .event { border: 1px solid var(--border); border-radius: 6px; margin: 16px 0; padding: 12px 16px; min-height: 420px; }
          .event-head { display: flex; flex-wrap: wrap; align-items: baseline; gap: 6px 16px; }
          .badge { font-size: 12px; padding: 1px 8px; border-radius: 10px; background: #e5e7eb; }
          .badge.DEL { background: #fee2e2; color: #991b1b; }
          .badge.DUP { background: #dcfce7; color: #166534; }
          .facets { display: grid; grid-template-columns: repeat(auto-fill, minmax(230px, 1fr)); gap: 8px; }
          .facet { border: 1px solid var(--border); cursor: pointer; }
          .facet:hover { border-color: var(--accent); }
          .facet .label { background: #e5e7eb; font-size: 12px; text-align: center; padding: 1px 0; }
          .controls { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; margin: 8px 0; }
          .controls input[type=text] { font: inherit; padding: 6px 9px; width: 340px; max-width: 100%;
                                       border: 1px solid #d1d5db; border-radius: 4px; }
          .controls select, .controls button { font: inherit; padding: 5px 10px; border: 1px solid #d1d5db; border-radius: 4px;
                                               background: #fff; cursor: pointer; }
          .controls button.primary { background: var(--accent); border-color: var(--accent); color: #fff; }
          .controls button:disabled { opacity: .5; cursor: default; }
          .err { color: var(--err); }
          #status { padding: 40px 0; }
          details summary { cursor: pointer; margin: 18px 0 8px; font-weight: 600; }
          tr.pass td:first-child { font-weight: 600; }
        </style>
        </head>
        <body>
        <header class="wrap">
          <h1>__TITLE__</h1>
          <div class="muted" id="summary"></div>
        </header>
        <nav><div class="wrap">
          <button data-tab="overview" class="active">Overview</button>
          <button data-tab="events">Passing CNV Calls</button>
          <button data-tab="browser">Interval Browser</button>
        </div></nav>
        <main class="wrap">
          <div id="status" class="muted">Loading data&hellip;</div>

          <section class="tab" id="tab-overview">
            <h2>Denoised copy ratio by chromosome</h2>
            <div class="muted">Case sample only. Shaded bars mark passing calls. Click a chromosome to open it in the Interval Browser.</div>
            <div class="facets" id="facets"></div>
            <h2>Adjusted read counts by GC content</h2>
            <div id="gc-host"></div>
          </section>

          <section class="tab" id="tab-events">
            <div id="events-summary"></div>
            <div id="events-list"></div>
            <div id="filtered-calls"></div>
          </section>

          <section class="tab" id="tab-browser">
            <h2>Interval Browser</h2>
            <div class="muted">Enter any interval to plot the case and panel data there, for example
              <span id="example-locus"></span>. Accepted formats: <code>chr:start-end</code>, <code>chr start end</code>,
              <code>chr:position</code> or a whole contig name. Drag across a plot to zoom in.</div>
            <form class="controls" id="locus-form">
              <input type="text" id="locus-input" spellcheck="false" autocomplete="off" aria-label="Interval">
              <label>Flanking context
                <select id="context-select">
                  <option value="auto">Auto</option>
                  <option value="0">None</option>
                  <option value="0.5">50% of width each side</option>
                  <option value="1">100% of width each side</option>
                  <option value="5">500% of width each side</option>
                </select>
              </label>
              <button type="submit" class="primary">Plot</button>
            </form>
            <div class="controls" id="nav-controls">
              <button type="button" data-nav="left" title="Pan left">&larr;</button>
              <button type="button" data-nav="right" title="Pan right">&rarr;</button>
              <button type="button" data-nav="in" title="Zoom in">Zoom in</button>
              <button type="button" data-nav="out" title="Zoom out">Zoom out</button>
              <button type="button" data-nav="reset" title="Back to the entered interval">Reset</button>
              <span id="view-readout" class="muted"></span>
            </div>
            <div id="locus-error" class="err"></div>
            <div id="browser-plots"></div>
            <div id="browser-calls"></div>
          </section>
        </main>

        <script type="application/json" id="gcnv-meta">__META__</script>
        __DATA__

        <script>
        (function () {
        "use strict";

        const META = JSON.parse(document.getElementById("gcnv-meta").textContent);
        const PANEL_RGB = [37, 99, 235];
        const CASE_COLOR = "#111827";
        const CALL_RGB = { DEL: "220,38,38", DUP: "22,163,74" };
        const FILTERED_RGB = "107,114,128";
        const QUERY_RGB = "202,138,4";
        const REGION_YMAX = 7;
        const OVERVIEW_YMAX = 5;
        // Flanking context rules for call plots, ported from the original R report.
        // Widths below minCallWidth / minQueryWidth are treated as that width, so short calls and single positions get context.
        const CONTEXT = { minPoints: 10, minWidthFactor: 2.5, maxWidthFactor: 10, minCallWidth: 1000, minQueryWidth: 10000 };

        let D = null;
        const contigByName = new Map();
        let TRACKS = [];

        // ---------------------------------------------------------------- helpers
        const fmt = v => (Math.round(v) || 0).toLocaleString("en-US");
        const esc = s => String(s).replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
        const el = (tag, attrs, html) => { const e = document.createElement(tag); Object.assign(e, attrs || {}); if (html != null) e.innerHTML = html; return e; };
        const locusText = (contig, start, end) => contig + ":" + fmt(start) + "-" + fmt(end);
        function fmtSize(bp) {
          if (bp >= 1e6) return (bp / 1e6).toFixed(bp >= 1e7 ? 1 : 2) + " Mb";
          if (bp >= 1e3) return (bp / 1e3).toFixed(bp >= 1e4 ? 1 : 2) + " kb";
          return fmt(bp) + " bp";
        }
        function firstIndex(lo, hi, pred) {
          while (lo < hi) { const m = (lo + hi) >> 1; if (pred(m)) hi = m; else lo = m + 1; }
          return lo;
        }
        // Indices [lo, hi) of intervals on the contig that overlap [ws, we].
        function indexRange(ci, ws, we) {
          const a = ci.offset, b = ci.offset + ci.count;
          const lo = firstIndex(a, b, i => D.ends[i] >= ws);
          return [lo, firstIndex(lo, b, i => D.starts[i] > we)];
        }
        function niceTicks(a, b, target) {
          const raw = (b - a) / Math.max(1, target);
          const p = Math.pow(10, Math.floor(Math.log10(raw)));
          const f = raw / p;
          const step = (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p;
          const out = [];
          for (let k = Math.ceil(a / step); k * step <= b; k++) out.push(k * step);
          return out;
        }

        // ---------------------------------------------------------------- data
        function b64ToBytes(b64) {
          const bin = atob(b64), out = new Uint8Array(bin.length);
          for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
          return out;
        }
        async function loadArray(name, Type) {
          const node = document.getElementById("gcnv-" + name);
          if (!node) return null;
          const stream = new Blob([b64ToBytes(node.textContent.trim())]).stream().pipeThrough(new DecompressionStream("gzip"));
          const buffer = await new Response(stream).arrayBuffer();
          node.textContent = "";
          return new Type(buffer);
        }
        async function loadData() {
          const [startDelta, lengths, caseCr, caseAdj, panelCr, panelAdj, gc] = await Promise.all([
            loadArray("start_delta", Int32Array), loadArray("length", Int32Array),
            loadArray("case_cr", Float32Array), loadArray("case_adj", Float32Array),
            loadArray("panel_cr", Uint8Array), loadArray("panel_adj", Uint8Array), loadArray("gc", Float32Array)]);
          const n = META.n_intervals;
          const starts = new Float64Array(n), ends = new Float64Array(n), mids = new Float64Array(n);
          let pos = 0;
          for (let i = 0; i < n; i++) {
            pos += startDelta[i];
            starts[i] = pos; ends[i] = pos + lengths[i]; mids[i] = pos + lengths[i] / 2;
          }
          return { starts, ends, mids, caseCr, caseAdj, panelCr, panelAdj, gc };
        }

        // ---------------------------------------------------------------- windows
        // Grow the window around [start, end] one interval at a time (taking whichever side is the smaller step) until it
        // is at least minWidthFactor x the width and holds minPoints flanking intervals, without exceeding maxWidthFactor x.
        function autoWindow(ci, start, end, minWidth) {
          const a = ci.offset, b = ci.offset + ci.count, S = D.starts;
          const firstInside = firstIndex(a, b, i => S[i] >= start);
          const firstAfter = firstIndex(a, b, i => S[i] > end);
          const totalBefore = firstInside - a, totalAfter = b - firstAfter;
          const before = n => n > 0 ? S[firstInside - n] : start;
          const after = n => n > 0 ? S[firstAfter + n - 1] : end;
          const width = Math.max(end - start, minWidth);
          let nb = 0, na = 0, prevNb = 0, prevNa = 0;
          while ((nb < totalBefore || na < totalAfter) &&
                 (after(na) - before(nb) < CONTEXT.minWidthFactor * width || na + nb < CONTEXT.minPoints) &&
                 after(na) - before(nb) < CONTEXT.maxWidthFactor * width) {
            prevNb = nb; prevNa = na;
            if (nb === totalBefore) na++;
            else if (na === totalAfter) nb++;
            else if (before(nb) - before(nb + 1) < after(na + 1) - after(na)) nb++;
            else na++;
          }
          if (after(na) - before(nb) > CONTEXT.maxWidthFactor * width) { nb = prevNb; na = prevNa; }
          return { ws: Math.min(before(nb), start), we: Math.max(na > 0 ? D.ends[firstAfter + na - 1] : end, end) };
        }
        function contextWindow(ci, start, end, mode) {
          if (mode === "auto") return autoWindow(ci, start, end, CONTEXT.minQueryWidth);
          const pad = (end - start + 1) * Number(mode);
          let ws = start - pad, we = end + pad;
          if (we - ws < 200) { const c = (ws + we) / 2; ws = c - 100; we = c + 100; }
          return { ws: Math.max(1, ws), we };
        }

        // ---------------------------------------------------------------- drawing
        function diskOffsets(r) {
          const out = [], R = Math.ceil(r);
          for (let dy = -R; dy <= R; dy++) for (let dx = -R; dx <= R; dx++) if (dx * dx + dy * dy <= r * r + 0.3) out.push(dx, dy);
          return out.length ? out : [0, 0];
        }
        // Accumulates overlapping semi-transparent points into a per-pixel count so that millions of panel points can be
        // drawn quickly; each pixel's opacity matches stacking that many points of the given alpha.
        class Density {
          constructor(w, h, radius) { this.w = Math.max(0, w); this.h = Math.max(0, h); this.cnt = new Uint16Array(this.w * this.h); this.offs = diskOffsets(radius); }
          add(x, y) {
            const xi = Math.round(x), yi = Math.round(y), w = this.w, h = this.h, c = this.cnt, o = this.offs;
            for (let k = 0; k < o.length; k += 2) {
              const xx = xi + o[k], yy = yi + o[k + 1];
              if (xx >= 0 && xx < w && yy >= 0 && yy < h) { const p = yy * w + xx; if (c[p] < 65535) c[p]++; }
            }
          }
          draw(ctx, dx, dy, rgb, alpha) {
            if (!this.w || !this.h) return;
            const lut = new Uint8ClampedArray(256);
            for (let k = 1; k < 256; k++) lut[k] = Math.round(255 * (1 - Math.pow(1 - alpha, k)));
            const img = new ImageData(this.w, this.h), d = img.data, c = this.cnt;
            for (let p = 0; p < c.length; p++) {
              const v = c[p];
              if (v) { const q = p * 4; d[q] = rgb[0]; d[q + 1] = rgb[1]; d[q + 2] = rgb[2]; d[q + 3] = lut[v > 255 ? 255 : v]; }
            }
            const off = document.createElement("canvas");
            off.width = this.w; off.height = this.h;
            off.getContext("2d").putImageData(img, 0, 0);
            ctx.save(); ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.drawImage(off, dx, dy); ctx.restore();
          }
        }

        function makePlotHost(height) {
          const host = el("div", { className: "plot" });
          host.dataset.height = height;
          host.append(el("canvas"), el("div", { className: "guide" }), el("div", { className: "sel" }), el("div", { className: "tip" }));
          return host;
        }
        function releasePlot(host) {
          const c = host.querySelector("canvas");
          c.width = 0; c.height = 0;
          host._s = null;
        }
        function setupCanvas(host) {
          const W = host.clientWidth, H = Number(host.dataset.height);
          if (!W) return null;
          const dpr = window.devicePixelRatio || 1, canvas = host.querySelector("canvas");
          canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr); canvas.style.height = H + "px";
          const ctx = canvas.getContext("2d");
          ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
          ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, W, H);
          return { ctx, dpr, W, H };
        }
        function drawAxes(ctx, M, pw, ph, xTicks, X, fmtX, yTicks, Y, xLabel, yLabel, fontSize) {
          ctx.save();
          ctx.font = fontSize + "px sans-serif";
          ctx.strokeStyle = "#374151"; ctx.fillStyle = "#374151"; ctx.lineWidth = 1;
          ctx.beginPath();
          ctx.moveTo(M.l + 0.5, M.t); ctx.lineTo(M.l + 0.5, M.t + ph + 0.5); ctx.lineTo(M.l + pw, M.t + ph + 0.5);
          ctx.stroke();
          ctx.textAlign = "center"; ctx.textBaseline = "top";
          for (const v of xTicks) {
            const x = Math.round(X(v)) + 0.5;
            ctx.beginPath(); ctx.moveTo(x, M.t + ph); ctx.lineTo(x, M.t + ph + 4); ctx.stroke();
            ctx.fillText(fmtX(v), x, M.t + ph + 6);
          }
          ctx.textAlign = "right"; ctx.textBaseline = "middle";
          for (const v of yTicks) {
            const y = Math.round(Y(v)) + 0.5;
            ctx.beginPath(); ctx.moveTo(M.l - 4, y); ctx.lineTo(M.l, y); ctx.stroke();
            ctx.fillText(String(v), M.l - 6, y);
          }
          if (xLabel) { ctx.textAlign = "center"; ctx.textBaseline = "bottom"; ctx.fillText(xLabel, M.l + pw / 2, M.t + ph + M.b - 2); }
          if (yLabel) {
            ctx.translate(12, M.t + ph / 2); ctx.rotate(-Math.PI / 2);
            ctx.textAlign = "center"; ctx.textBaseline = "middle"; ctx.fillText(yLabel, 0, 0);
          }
          ctx.restore();
        }
        function drawGrid(ctx, M, pw, ph, yTicks, Y, refLines) {
          ctx.save(); ctx.lineWidth = 1;
          ctx.strokeStyle = "#eef0f3";
          for (const v of yTicks) { const y = Math.round(Y(v)) + 0.5; ctx.beginPath(); ctx.moveTo(M.l, y); ctx.lineTo(M.l + pw, y); ctx.stroke(); }
          ctx.strokeStyle = "rgba(55,65,81,.7)"; ctx.setLineDash([5, 4]);
          for (const v of refLines) { const y = Math.round(Y(v)) + 0.5; ctx.beginPath(); ctx.moveTo(M.l, y); ctx.lineTo(M.l + pw, y); ctx.stroke(); }
          ctx.restore();
        }

        // Shaded rectangles for calls and the queried interval in view.
        function callMarks(ci, ws, we, focus) {
          const marks = [];
          for (const c of META.calls) {
            if (c.contig !== ci.name || c.end < ws || c.start > we) continue;
            const rgb = c.pass ? (CALL_RGB[c.alt] || FILTERED_RGB) : FILTERED_RGB;
            const isFocus = focus && focus.call === c;
            marks.push({ start: c.start, end: c.end, fill: "rgba(" + rgb + "," + (c.pass ? 0.16 : 0.1) + ")",
                         stroke: c.pass ? null : "rgba(" + rgb + ",.6)", dash: c.pass ? null : [4, 3], bold: isFocus,
                         label: c.alt + " CN=" + (c.cn == null ? "?" : c.cn) + (c.pass ? "" : " (" + c.filter + ")") });
          }
          if (focus && focus.query) {
            marks.push({ start: focus.query.start, end: focus.query.end, fill: "rgba(" + QUERY_RGB + ",.1)",
                         stroke: "rgba(" + QUERY_RGB + ",.9)", dash: [6, 3], label: "entered interval" });
          }
          return marks;
        }

        // Scatter of panel values (blue) with the case (black points joined by a line) over a genomic window.
        function renderRegion(host, track, ci, view, marks) {
          const cv = setupCanvas(host);
          if (!cv) return;
          const { ctx, dpr, W, H } = cv;
          const M = { l: 54, r: 12, t: 22, b: 40 }, pw = W - M.l - M.r, ph = H - M.t - M.b, yMax = REGION_YMAX;
          const pad = Math.max(view.we - view.ws, 1) * 0.02, x0 = view.ws - pad, x1 = view.we + pad;
          const X = v => M.l + (v - x0) / (x1 - x0) * pw, Y = v => M.t + ph - v / yMax * ph;
          const [lo, hi] = indexRange(ci, x0, x1), n = hi - lo;
          const yTicks = [0, 1, 2, 3, 4, 5, 6, 7];

          ctx.save(); ctx.beginPath(); ctx.rect(M.l, M.t, pw, ph); ctx.clip();
          for (const m of marks) {
            const a = X(m.start), w = Math.max(X(m.end) - a, 1.5);
            ctx.fillStyle = m.fill; ctx.fillRect(a, M.t, w, ph);
            if (m.stroke || m.bold) {
              ctx.strokeStyle = m.stroke || "rgba(17,24,39,.55)"; ctx.lineWidth = m.bold ? 1.5 : 1; ctx.setLineDash(m.dash || []);
              ctx.strokeRect(a + 0.5, M.t + 0.5, w - 1, ph - 1);
            }
          }
          ctx.restore();
          drawGrid(ctx, M, pw, ph, yTicks, Y, [1, 2, 3]);

          if (track.P > 0 && n > 0) {
            const radius = n > 3000 ? 0.8 : n > 800 ? 1.4 : 2.3;
            const dens = new Density(Math.round(pw * dpr), Math.round(ph * dpr), radius * dpr);
            const P = track.P, pv = track.panel, sc = META.panel_scale, miss = META.panel_missing;
            const kx = pw * dpr / (x1 - x0), ky = ph * dpr / yMax;
            for (let i = lo; i < hi; i++) {
              const px = (D.mids[i] - x0) * kx, base = i * P;
              for (let j = 0; j < P; j++) {
                const q = pv[base + j];
                if (q === miss) break;  // sorted, missing values last
                const v = q / sc;
                if (v > yMax) break;
                dens.add(px, (yMax - v) * ky);
              }
            }
            dens.draw(ctx, Math.round(M.l * dpr), Math.round(M.t * dpr), PANEL_RGB, 0.2);
          }

          ctx.save(); ctx.beginPath(); ctx.rect(M.l, M.t, pw, ph); ctx.clip();
          const cvals = track.caseVals;
          ctx.strokeStyle = CASE_COLOR; ctx.lineWidth = 1; ctx.beginPath();
          let pen = false;
          for (let i = lo; i < hi; i++) {
            const v = cvals[i];
            if (!Number.isFinite(v)) { pen = false; continue; }
            const x = X(D.mids[i]), y = Y(Math.min(v, yMax + 1));
            if (pen) ctx.lineTo(x, y); else ctx.moveTo(x, y);
            pen = true;
          }
          ctx.stroke();
          const pr = n > 4000 ? 0 : n > 1000 ? 1.2 : 2.4;
          ctx.fillStyle = CASE_COLOR; ctx.beginPath();
          for (let i = lo; i < hi; i++) {
            const v = cvals[i];
            if (!Number.isFinite(v)) continue;
            const x = X(D.mids[i]);
            if (v > yMax) { ctx.moveTo(x - 4, M.t + 7); ctx.lineTo(x + 4, M.t + 7); ctx.lineTo(x, M.t + 1); ctx.closePath(); }
            else if (pr) { const y = Y(v); ctx.moveTo(x + pr, y); ctx.arc(x, y, pr, 0, 2 * Math.PI); }
          }
          ctx.fill();
          ctx.restore();

          // Labels go in the first lane (above the plot, then inside its top edge) where they do not collide.
          ctx.save(); ctx.font = "11px sans-serif"; ctx.textBaseline = "bottom";
          const laneY = [M.t - 4, M.t + 14, M.t + 28], lanes = laneY.map(() => []);
          for (const m of marks) {
            if (!m.label) continue;
            const tw = ctx.measureText(m.label).width, x = Math.min(Math.max(X(m.start), M.l) + 2, W - M.r - tw);
            const lane = lanes.findIndex(used => used.every(([a, b]) => x > b + 8 || x + tw < a - 8));
            if (lane < 0) continue;
            lanes[lane].push([x, x + tw]);
            ctx.fillStyle = m.stroke || "#374151";
            ctx.fillText(m.label, x, laneY[lane]);
          }
          ctx.restore();

          if (n === 0) {
            ctx.save(); ctx.fillStyle = "#6b7280"; ctx.font = "13px sans-serif"; ctx.textAlign = "center";
            ctx.fillText("No gCNV intervals in this region", M.l + pw / 2, M.t + ph / 2); ctx.restore();
          }
          drawAxes(ctx, M, pw, ph, niceTicks(Math.max(x0, 0), x1, Math.max(2, Math.floor(pw / 115))), X, fmt, yTicks, Y,
                   "Position on " + ci.name, track.yLabel, 11);
          host._s = { kind: "region", ci, track, x0, x1, M, pw, ph, lo, hi };
        }

        // Small multiple of the case copy ratio across a whole contig.
        function renderFacet(host, ci) {
          const cv = setupCanvas(host);
          if (!cv) return;
          const { ctx, dpr, W, H } = cv;
          const M = { l: 22, r: 6, t: 4, b: 6 }, pw = W - M.l - M.r, ph = H - M.t - M.b, yMax = OVERVIEW_YMAX;
          const a = ci.offset, b = ci.offset + ci.count;
          const x0 = D.starts[a], x1 = Math.max(D.ends[b - 1], x0 + 1);
          const X = v => M.l + (v - x0) / (x1 - x0) * pw, Y = v => M.t + ph - v / yMax * ph;
          for (const c of META.calls) {
            if (!c.pass || c.contig !== ci.name) continue;
            ctx.fillStyle = "rgba(" + (CALL_RGB[c.alt] || FILTERED_RGB) + ",.3)";
            ctx.fillRect(X(c.start), M.t, Math.max(X(c.end) - X(c.start), 2), ph);
          }
          drawGrid(ctx, M, pw, ph, [0, 1, 2, 3, 4, 5], Y, []);
          const dens = new Density(Math.round(pw * dpr), Math.round(ph * dpr), 0.6 * dpr);
          const kx = pw * dpr / (x1 - x0), ky = ph * dpr / yMax;
          for (let i = a; i < b; i++) {
            const v = D.caseCr[i];
            if (Number.isFinite(v) && v <= yMax) dens.add((D.mids[i] - x0) * kx, (yMax - v) * ky);
          }
          dens.draw(ctx, Math.round(M.l * dpr), Math.round(M.t * dpr), [17, 24, 39], 0.2);
          drawAxes(ctx, M, pw, ph, [], X, fmt, [0, 1, 2, 3, 4, 5], Y, null, null, 9);
        }

        function renderGc(host) {
          const cv = setupCanvas(host);
          if (!cv) return;
          const { ctx, dpr, W, H } = cv;
          const M = { l: 54, r: 12, t: 10, b: 40 }, pw = W - M.l - M.r, ph = H - M.t - M.b, yMax = OVERVIEW_YMAX;
          const X = v => M.l + v * pw, Y = v => M.t + ph - v / yMax * ph;
          drawGrid(ctx, M, pw, ph, [0, 1, 2, 3, 4, 5], Y, []);
          const dens = new Density(Math.round(pw * dpr), Math.round(ph * dpr), 0.6 * dpr);
          for (let i = 0; i < META.n_intervals; i++) {
            const g = D.gc[i], v = D.caseAdj[i];
            if (Number.isFinite(g) && Number.isFinite(v) && v <= yMax) dens.add(g * pw * dpr, (yMax - v) * ph * dpr / yMax);
          }
          dens.draw(ctx, Math.round(M.l * dpr), Math.round(M.t * dpr), [17, 24, 39], 0.2);
          drawAxes(ctx, M, pw, ph, [0, 0.2, 0.4, 0.6, 0.8, 1], X, v => v.toFixed(1), [0, 1, 2, 3, 4, 5], Y,
                   "GC content", "Adjusted read counts", 11);
        }

        // ---------------------------------------------------------------- hover and drag-to-zoom
        function panelSummary(track, i) {
          if (!track.P) return null;
          const vals = [];
          for (let j = 0; j < track.P; j++) {
            const q = track.panel[i * track.P + j];
            if (q === META.panel_missing) break;  // each interval's values are sorted, missing values last
            vals.push(q / META.panel_scale);
          }
          if (!vals.length) return null;
          const m = vals.length >> 1;
          const median = vals.length % 2 ? vals[m] : (vals[m - 1] + vals[m]) / 2;
          return "median " + median.toFixed(2) + ", range " + vals[0].toFixed(2) + "&ndash;" + vals[vals.length - 1].toFixed(2) + " (n=" + vals.length + ")";
        }
        function tooltipHtml(s, i) {
          const f = v => Number.isFinite(v) ? v.toFixed(3) : "NA";
          const rows = ["<b>" + esc(s.ci.name) + ":" + fmt(D.starts[i]) + "-" + fmt(D.ends[i]) + "</b>",
                        "Case copy ratio: " + f(D.caseCr[i]), "Case adjusted counts: " + f(D.caseAdj[i])];
          if (D.gc) rows.push("GC content: " + f(D.gc[i]));
          const p = panelSummary(s.track, i);
          if (p) rows.push("Panel " + (s.track.key === "cr" ? "copy ratio" : "adjusted counts") + ": " + p);
          return rows.join("<br>");
        }
        function attachInteractions(host, onZoom) {
          const tip = host.querySelector(".tip"), guide = host.querySelector(".guide"), sel = host.querySelector(".sel");
          let dragX = null;
          const rel = e => { const r = host.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };
          const toPos = (s, x) => s.x0 + (x - s.M.l) / s.pw * (s.x1 - s.x0);
          const hide = () => { tip.style.display = "none"; guide.style.display = "none"; };
          host.addEventListener("mousemove", e => {
            const s = host._s;
            if (!s) return;
            const [mx, my] = rel(e);
            if (dragX !== null) {
              const x = Math.min(Math.max(mx, s.M.l), s.M.l + s.pw);
              Object.assign(sel.style, { display: "block", left: Math.min(dragX, x) + "px", width: Math.abs(x - dragX) + "px",
                                         top: s.M.t + "px", height: s.ph + "px" });
            }
            if (mx < s.M.l || mx > s.M.l + s.pw || my < s.M.t || my > s.M.t + s.ph || s.hi <= s.lo) { hide(); return; }
            const pos = toPos(s, mx);
            let i = firstIndex(s.lo, s.hi, k => D.mids[k] >= pos);
            if (i === s.hi || (i > s.lo && pos - D.mids[i - 1] < D.mids[i] - pos)) i--;
            const gx = s.M.l + (D.mids[i] - s.x0) / (s.x1 - s.x0) * s.pw;
            Object.assign(guide.style, { display: "block", left: Math.round(gx) + "px", top: s.M.t + "px", height: s.ph + "px" });
            tip.innerHTML = tooltipHtml(s, i);
            tip.style.display = "block";
            const tw = tip.offsetWidth, left = mx + 14 + tw > host.clientWidth ? mx - 14 - tw : mx + 14;
            tip.style.left = Math.max(0, left) + "px";
            tip.style.top = Math.max(0, my - tip.offsetHeight - 8) + "px";
          });
          host.addEventListener("mouseleave", hide);
          if (!onZoom) return;
          host.style.cursor = "crosshair";
          host.addEventListener("mousedown", e => {
            const s = host._s;
            if (!s || e.button !== 0) return;
            const [mx] = rel(e);
            if (mx < s.M.l || mx > s.M.l + s.pw) return;
            dragX = mx;
            e.preventDefault();
          });
          window.addEventListener("mouseup", e => {
            if (dragX === null) return;
            const s = host._s, [mx] = rel(e), a = dragX;
            dragX = null; sel.style.display = "none";
            if (!s || Math.abs(mx - a) < 6) return;
            const x = Math.min(Math.max(mx, s.M.l), s.M.l + s.pw);
            onZoom(toPos(s, Math.min(a, x)), toPos(s, Math.max(a, x)));
          });
        }

        function legendHtml(focusText) {
          const nCr = META.panel_cr_samples.length, nAdj = META.panel_adj_samples.length;
          const panel = nCr === nAdj ? nCr + " samples" : nCr + " copy ratio / " + nAdj + " read count samples";
          return '<div class="legend">' +
            '<span><i class="sw" style="background:#111827"></i>Case (' + esc(META.sample) + ')</span>' +
            '<span><i class="sw" style="background:rgba(37,99,235,.45)"></i>Panel (' + panel + ')</span>' +
            '<span><i class="sw rect" style="background:rgba(' + CALL_RGB.DEL + ',.3)"></i>Passing DEL</span>' +
            '<span><i class="sw rect" style="background:rgba(' + CALL_RGB.DUP + ',.3)"></i>Passing DUP</span>' +
            '<span><i class="sw rect" style="background:rgba(' + FILTERED_RGB + ',.25);border:1px dashed #6b7280"></i>Filtered call</span>' +
            (focusText || "") + "</div>";
        }

        // A pair of plots (copy ratio and read counts) for one window.
        function makeRegionPair(parent, onZoom) {
          const pair = el("div", { className: "pair" });
          const hosts = TRACKS.map(t => {
            const cell = el("div");
            cell.append(el("h4", {}, esc(t.title)));
            const host = makePlotHost(300);
            attachInteractions(host, onZoom);
            cell.append(host);
            pair.append(cell);
            return host;
          });
          parent.append(pair);
          return hosts;
        }

        // ---------------------------------------------------------------- tables
        function callRow(c, i, actions) {
          const size = c.end - c.start + 1;
          return "<tr class='" + (c.pass ? "pass" : "") + "'>" +
            (i != null ? "<td class='num'>" + (i + 1) + "</td>" : "") +
            "<td>" + esc(locusText(c.contig, c.start, c.end)) + "</td>" +
            "<td><span class='badge " + esc(c.alt) + "'>" + esc(c.alt) + "</span></td>" +
            "<td class='num'>" + esc(c.cn == null ? "" : c.cn) + "</td>" +
            "<td class='num'>" + (c.qual == null ? "" : fmt(c.qual)) + "</td>" +
            "<td class='num'>" + fmtSize(size) + "</td>" +
            "<td class='num'>" + (c.np == null ? "" : c.np) + "</td>" +
            "<td class='num'>" + (c.panel_freq == null ? "" : c.panel_freq.toFixed(3)) + "</td>" +
            "<td class='num'>" + (c.panel_count == null ? "" : c.panel_count) + "</td>" +
            "<td>" + esc(c.filter) + "</td>" +
            "<td>" + actions + "</td></tr>";
        }
        function callTable(calls, numbered, actionsFor) {
          return "<div class='table-scroll'><table><thead><tr>" + (numbered ? "<th class='num'>#</th>" : "") +
            "<th>Interval</th><th>Type</th><th class='num'>CN</th><th class='num'>QUAL</th><th class='num'>Size</th>" +
            "<th class='num'>Intervals</th><th class='num'>PANEL_FREQ</th><th class='num'>PANEL_COUNT</th><th>FILTER</th><th></th>" +
            "</tr></thead><tbody>" + calls.map((c, i) => callRow(c, numbered ? i : null, actionsFor(c, i))).join("") +
            "</tbody></table></div>";
        }
        function browseLink(c) {
          return "<button class='link' data-browse='" + esc(c.contig + ":" + c.start + "-" + c.end) + "'>Interval Browser</button>";
        }

        // ---------------------------------------------------------------- overview tab
        let overviewBuilt = false;
        function renderOverview() {
          const facets = document.getElementById("facets"), gcHost = document.getElementById("gc-host");
          if (!overviewBuilt) {
            overviewBuilt = true;
            for (const ci of META.contigs) {
              const box = el("div", { className: "facet", title: "Open " + ci.name + " in the Interval Browser" });
              box.append(el("div", { className: "label" }, esc(ci.name)));
              const host = makePlotHost(130);
              box.append(host);
              box.addEventListener("click", () => openBrowser(ci.name, "0"));
              facets.append(box);
              host._render = () => renderFacet(host, ci);
            }
            if (D.gc) {
              const host = makePlotHost(380);
              host.style.maxWidth = "760px";
              gcHost.append(host);
              host._render = () => renderGc(host);
            } else {
              gcHost.innerHTML = "<div class='muted'>GC content is not available: gCNV was run without annotated intervals.</div>";
            }
          }
          document.querySelectorAll("#tab-overview .plot").forEach(h => h._render && h._render());
        }

        // ---------------------------------------------------------------- passing calls tab
        function buildEvents() {
          const passing = META.calls.filter(c => c.pass), filtered = META.calls.filter(c => !c.pass);
          const summary = document.getElementById("events-summary"), list = document.getElementById("events-list");
          summary.innerHTML = "<h2>Passing CNV calls (" + passing.length + ")</h2>" + (passing.length
            ? callTable(passing, true, (c, i) => "<a href='#event-" + (i + 1) + "'>Plots</a> &middot; " + browseLink(c))
            : "<div class='muted'>No CNV calls passed filters. Use the Interval Browser to inspect any region.</div>");
          if (filtered.length) {
            document.getElementById("filtered-calls").innerHTML = "<details><summary>Filtered calls (" + filtered.length +
              ")</summary>" + callTable(filtered, false, c => browseLink(c)) + "</details>";
          }
          const observer = new IntersectionObserver(entries => {
            for (const entry of entries) {
              const hosts = entry.target._hosts;
              if (entry.isIntersecting && !entry.target._rendered) { entry.target._rendered = true; hosts.forEach(h => h._render()); }
              else if (!entry.isIntersecting && entry.target._rendered) { entry.target._rendered = false; hosts.forEach(releasePlot); }
            }
          }, { rootMargin: "900px 0px" });
          passing.forEach((c, i) => {
            const ci = contigByName.get(c.contig);
            const box = el("div", { className: "event", id: "event-" + (i + 1) });
            box.append(el("div", { className: "event-head" },
              "<h3>" + (i + 1) + ". " + esc(locusText(c.contig, c.start, c.end)) + "</h3>" +
              "<span class='badge " + esc(c.alt) + "'>" + esc(c.alt) + "</span>" +
              "<span class='muted'>CN " + esc(c.cn == null ? "?" : c.cn) + " &middot; QUAL " + (c.qual == null ? "NA" : fmt(c.qual)) +
              " &middot; " + fmtSize(c.end - c.start + 1) + (c.panel_freq == null ? "" : " &middot; PANEL_FREQ " + c.panel_freq.toFixed(3)) +
              "</span>" + browseLink(c)));
            list.append(box);
            if (!ci) { box.append(el("div", { className: "err" }, "Contig " + esc(c.contig) + " has no gCNV intervals.")); return; }
            box.append(el("div", {}, legendHtml()));
            const view = autoWindow(ci, c.start, c.end, CONTEXT.minCallWidth);
            const marks = callMarks(ci, view.ws, view.we, { call: c });
            const hosts = makeRegionPair(box, null);
            hosts.forEach((h, k) => { h._render = () => renderRegion(h, TRACKS[k], ci, view, marks); });
            box._hosts = hosts;
            observer.observe(box);
          });
        }
        function rerenderEvents() {
          document.querySelectorAll("#events-list .event").forEach(b => { if (b._rendered) b._hosts.forEach(h => h._render()); });
        }

        // ---------------------------------------------------------------- interval browser tab
        const browser = { query: null, view: null, hosts: null };
        function resolveContig(name) {
          if (contigByName.has(name)) return contigByName.get(name);
          const alt = name.toLowerCase().startsWith("chr") ? name.slice(3) : "chr" + name;
          if (contigByName.has(alt)) return contigByName.get(alt);
          for (const ci of META.contigs) if (ci.name.toLowerCase() === name.toLowerCase() || ci.name.toLowerCase() === alt.toLowerCase()) return ci;
          return null;
        }
        function parseLocus(text) {
          const s = text.trim().replace(/[\u2013\u2014]/g, "-");
          if (!s) throw new Error("Enter an interval, for example " + exampleLocus() + ".");
          let ci = resolveContig(s);
          if (ci) return { ci, start: D.starts[ci.offset], end: D.ends[ci.offset + ci.count - 1] };
          const m = s.replace(/,/g, "").match(/^(.+?)(?::|\s+)(\d+)(?:\s*(?:-|\s)\s*(\d+))?$/);
          if (!m) throw new Error("Could not read \"" + text.trim() + "\". Use chr:start-end, for example " + exampleLocus() + ".");
          ci = resolveContig(m[1].trim());
          if (!ci) {
            const names = META.contigs.map(c => c.name);
            throw new Error("Contig \"" + m[1].trim() + "\" has no gCNV intervals. Available: " + names.slice(0, 30).join(", ") + (names.length > 30 ? ", ..." : "") + ".");
          }
          let start = Number(m[2]), end = m[3] == null ? start : Number(m[3]);
          if (end < start) [start, end] = [end, start];
          return { ci, start, end };
        }
        function exampleLocus() {
          const c = META.calls.find(x => x.pass) || null;
          if (c) return c.contig + ":" + fmt(c.start) + "-" + fmt(c.end);
          const ci = META.contigs[0], i = ci.offset + (ci.count >> 1);
          return locusText(ci.name, D.starts[i], D.starts[i] + 100000);
        }
        function openBrowser(text, mode) {
          if (mode != null) document.getElementById("context-select").value = mode;
          document.getElementById("locus-input").value = text;
          activateTab("browser");
          browseTo(text);
        }
        function browseTo(text) {
          const err = document.getElementById("locus-error");
          let q;
          try { q = parseLocus(text); } catch (e) { err.textContent = e.message; return; }
          err.textContent = "";
          browser.query = q;
          const canonical = q.ci.name + ":" + q.start + "-" + q.end;
          history.replaceState(null, "", "#locus=" + encodeURIComponent(canonical));
          const input = document.getElementById("locus-input");
          if (document.activeElement !== input) input.value = locusText(q.ci.name, q.start, q.end);
          resetView();
        }
        function resetView() {
          const q = browser.query;
          if (!q) return;
          browser.view = contextWindow(q.ci, q.start, q.end, document.getElementById("context-select").value);
          renderBrowser();
        }
        function setView(ws, we) {
          if (we - ws < 200) { const c = (ws + we) / 2; ws = c - 100; we = c + 100; }
          browser.view = { ws: Math.max(1, ws), we };
          renderBrowser();
        }
        function renderBrowser() {
          const q = browser.query, v = browser.view;
          document.querySelectorAll("#nav-controls button").forEach(b => { b.disabled = !q; });
          if (!q) return;
          const container = document.getElementById("browser-plots");
          if (!browser.hosts) {
            container.innerHTML = legendHtml('<span><i class="sw rect" style="background:rgba(' + QUERY_RGB + ',.2);border:1px dashed rgb(' + QUERY_RGB + ')"></i>Entered interval</span>');
            browser.hosts = makeRegionPair(container, setView);
          }
          const whole = q.start === D.starts[q.ci.offset] && q.end === D.ends[q.ci.offset + q.ci.count - 1];
          const focus = whole ? {} : { query: q };
          const marks = callMarks(q.ci, v.ws, v.we, focus);
          browser.hosts.forEach((h, k) => renderRegion(h, TRACKS[k], q.ci, v, marks));
          const [lo, hi] = indexRange(q.ci, v.ws, v.we);
          document.getElementById("view-readout").textContent =
            "Showing " + locusText(q.ci.name, v.ws, v.we) + " (" + fmtSize(v.we - v.ws + 1) + ", " + fmt(hi - lo) + " gCNV intervals)";
          const calls = META.calls.filter(c => c.contig === q.ci.name && c.end >= v.ws && c.start <= v.we);
          document.getElementById("browser-calls").innerHTML = "<h2>Calls in view (" + calls.length + ")</h2>" +
            (calls.length ? callTable(calls, false, c => browseLink(c)) : "<div class='muted'>No CNV calls overlap this view.</div>");
        }
        function navigate(action) {
          const v = browser.view;
          if (!v) return;
          const span = v.we - v.ws, c = (v.ws + v.we) / 2;
          if (action === "reset") resetView();
          else if (action === "left") setView(v.ws - span / 2, v.we - span / 2);
          else if (action === "right") setView(v.ws + span / 2, v.we + span / 2);
          else if (action === "in") setView(c - span / 4, c + span / 4);
          else if (action === "out") setView(c - span * 1.5, c + span * 1.5);
        }

        // ---------------------------------------------------------------- tabs and startup
        let activeTab = "overview";
        function activateTab(name) {
          activeTab = name;
          document.querySelectorAll("nav button").forEach(b => b.classList.toggle("active", b.dataset.tab === name));
          document.querySelectorAll(".tab").forEach(s => s.classList.toggle("active", s.id === "tab-" + name));
          renderActive();
        }
        function renderActive() {
          if (!D) return;
          if (activeTab === "overview") renderOverview();
          else if (activeTab === "events") rerenderEvents();
          else renderBrowser();
        }

        async function start() {
          const status = document.getElementById("status");
          if (typeof DecompressionStream === "undefined") {
            status.innerHTML = "<span class='err'>This browser is too old to open the report (DecompressionStream is not supported). Please use a current version of Chrome, Edge, Firefox or Safari.</span>";
            return;
          }
          try { D = await loadData(); } catch (e) {
            status.innerHTML = "<span class='err'>Failed to load report data: " + esc(e.message) + "</span>";
            throw e;
          }
          status.remove();
          META.contigs.forEach(ci => contigByName.set(ci.name, ci));
          TRACKS = [
            { key: "cr", title: "Denoised Copy Ratio", yLabel: "Denoised linear copy ratio", caseVals: D.caseCr, panel: D.panelCr, P: META.panel_cr_samples.length },
            { key: "adj", title: "Adjusted Read Counts", yLabel: "Adjusted read counts", caseVals: D.caseAdj, panel: D.panelAdj, P: META.panel_adj_samples.length },
          ];
          const nPass = META.calls.filter(c => c.pass).length;
          document.getElementById("summary").textContent = "Sample " + META.sample + " \u00b7 " + nPass + " passing and " +
            (META.calls.length - nPass) + " filtered CNV calls \u00b7 " + fmt(META.n_intervals) + " gCNV intervals";
          document.querySelector("nav button[data-tab=events]").textContent = "Passing CNV Calls (" + nPass + ")";
          document.getElementById("example-locus").textContent = exampleLocus();
          document.getElementById("locus-input").placeholder = exampleLocus();

          document.querySelectorAll("nav button").forEach(b => b.addEventListener("click", () => activateTab(b.dataset.tab)));
          document.getElementById("locus-form").addEventListener("submit", e => { e.preventDefault(); browseTo(document.getElementById("locus-input").value); });
          document.getElementById("context-select").addEventListener("change", resetView);
          document.querySelectorAll("#nav-controls button").forEach(b => b.addEventListener("click", () => navigate(b.dataset.nav)));
          document.addEventListener("click", e => {
            const t = e.target.closest("[data-browse]");
            if (t) { e.preventDefault(); openBrowser(t.dataset.browse, "auto"); }
          });
          let resizeTimer = null;
          window.addEventListener("resize", () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(renderActive, 150); });

          buildEvents();
          const hash = decodeURIComponent(location.hash || "");
          if (hash.startsWith("#locus=")) { openBrowser(hash.slice(7)); }
          else if (hash.startsWith("#event-")) { activateTab("events"); document.getElementById(hash.slice(1))?.scrollIntoView(); }
          else activateTab("overview");
        }
        start();
        })();
        </script>
        </body>
        </html>
        """

        if __name__ == "__main__":
            main()
        EOF

        python3 gcnv_visualization.py \
            --vcf ~{filtered_vcf} \
            --case-copy-ratios ~{case_copy_ratios} \
            --case-read-counts ~{case_read_counts} \
            --panel-copy-ratios-list ~{write_lines(panel_copy_ratios)} \
            --panel-read-counts-list ~{write_lines(panel_read_counts)} \
            --interval-lists-list ~{write_lines(interval_lists)} \
            --title "~{output_prefix}" \
            --output ~{output_prefix}_cnv_event_report.html
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-h5py@sha256:7f8d59658a06c585f005ee0cc33356ac029ce742cb561117f43a9345e3c05536"
        disks: "local-disk 100 HDD"
        memory: mem_gb + " GB"
    }

    output {
        File cnv_event_report = "~{output_prefix}_cnv_event_report.html"
    }
}
