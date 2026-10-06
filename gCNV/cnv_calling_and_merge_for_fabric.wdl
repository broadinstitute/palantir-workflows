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

        # Panel samples to draw as lines, alongside the case, in the plots of the gCNV visualization report.
        Array[String] representative_sample_ids = []
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

    call GCNVVisualization {
        input:
            filtered_vcf = SingleSampleGCNVAndFilterVCFs.filtered_vcf,
            case_copy_ratios = SingleSampleGCNVAndFilterVCFs.denoised_copy_ratios,
            case_read_counts = SingleSampleGCNVAndFilterVCFs.read_counts,
            panel_copy_ratios = gcnv_panel_copy_ratios,
            panel_read_counts = gcnv_panel_read_counts,
            interval_lists = [SingleSampleGCNVAndFilterVCFs.interval_list],
            representative_sample_ids = representative_sample_ids
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
        File cnv_event_report = GCNVVisualization.cnv_event_report

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

task GCNVVisualization {
    input {
        File filtered_vcf
        File case_copy_ratios
        File case_read_counts
        Array[File]+ panel_copy_ratios
        Array[File]+ panel_read_counts
        Array[File]+ interval_lists
        Array[String] representative_sample_ids = []
        Int mem_gb=4
    }

    String output_prefix = basename(filtered_vcf, ".filtered.genotyped-segments.vcf.gz")

    # Produces a self-contained HTML report with plots for each passing CNV call and an Interval Browser tab in which
    # the user can enter any interval to plot it. The case and panel data for every interval, and the plotly.js library
    # used to draw the plots, are embedded in the HTML, so no server or internet connection is needed to view it.
    command <<<
        set -euo pipefail

        cat << 'EOF' > gcnv_visualization.py
        """Build a self-contained interactive HTML report of gCNV calls.

        The report embeds the case and panel denoised copy ratios and adjusted read counts for every
        interval, so plots are drawn in the browser with plotly.js (also embedded) both for each passing
        CNV call and for any interval the user enters in the Interval Browser tab. No server or internet
        connection is needed to view it.
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
        from plotly.offline import get_plotlyjs

        # Panel values are stored as uint8 round(value * PANEL_SCALE), capped at 254; 255 marks a missing value. Each
        # interval's panel values are sorted, which makes the embedded data several times smaller after compression; the
        # plots do not need to know which panel sample each value came from. Representative panel samples, which are drawn
        # as lines, are embedded separately at full precision.
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
            parser.add_argument("--representative-samples-list",
                                help="file listing panel sample names to draw as lines alongside the case")
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

            case_rc = read_normalized_counts(args.case_read_counts)
            # The case is left out of the panel if it is also a panel sample, so that the panel does not include its values.
            case_names = {case_sample, case_rc[4]}

            representatives = []
            for name in read_list(args.representative_samples_list) if args.representative_samples_list else []:
                if name in case_names:
                    print(f"Representative sample {name} is the case sample, which is always drawn; ignoring it")
                elif name not in representatives:
                    representatives.append(name)

            # Panel denoised copy ratios.
            panel_cr_paths = read_list(args.panel_copy_ratios_list)
            panel_cr = np.full((n_intervals, len(panel_cr_paths)), PANEL_MISSING, dtype=np.uint8)
            panel_cr_samples = []
            rep_cr = {}
            for path in panel_cr_paths:
                table, sample_name, _ = read_copy_ratios(path)
                if sample_name in case_names:
                    continue
                idx = align(table["CONTIG"], table["START"], table["END"])
                found = idx >= 0
                values = np.full(n_intervals, np.nan, dtype=np.float32)
                values[idx[found]] = table["LINEAR_COPY_RATIO"].to_numpy()[found]
                if sample_name in representatives:
                    rep_cr.setdefault(sample_name, values)
                panel_cr[:, len(panel_cr_samples)] = quantize(values)
                panel_cr_samples.append(sample_name)
            panel_cr = panel_cr[:, :len(panel_cr_samples)]

            # Read counts: normalize each sample by its mean count, then scale each interval by the panel mean so that
            # copy number 2 sits at 2.
            panel_rc_paths = read_list(args.panel_read_counts_list)
            panel_norm = np.full((len(panel_rc_paths), n_intervals), np.nan, dtype=np.float32)
            panel_rc_samples = []
            for path in panel_rc_paths:
                contig, start, end, norm_counts, sample_name = read_normalized_counts(path)
                if sample_name in case_names:
                    continue
                idx = align(contig, start, end)
                found = idx >= 0
                panel_norm[len(panel_rc_samples), idx[found]] = norm_counts[found]
                panel_rc_samples.append(sample_name)
            panel_norm = panel_norm[:len(panel_rc_samples)]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                panel_mean = np.nanmean(panel_norm, axis=0)
                panel_adjusted = 2 * panel_norm / panel_mean
            del panel_norm
            panel_adjusted[~np.isfinite(panel_adjusted)] = np.nan
            rep_adj = {name: panel_adjusted[panel_rc_samples.index(name)].copy()
                       for name in representatives if name in panel_rc_samples}
            panel_adj = np.ascontiguousarray(quantize(panel_adjusted).T)
            del panel_adjusted

            n_excluded = len(panel_cr_paths) - len(panel_cr_samples) + len(panel_rc_paths) - len(panel_rc_samples)
            if n_excluded:
                print(f"Left {n_excluded} file(s) of the case sample out of the panel")
            missing = [name for name in representatives if name not in rep_cr and name not in rep_adj]
            if missing:
                raise ValueError("Representative samples not found in the panel: " + ", ".join(missing))

            contig, start, end, norm_counts, _ = case_rc
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
            # A representative sample missing from one of the panels has no array for that panel.
            for k, name in enumerate(representatives):
                if name in rep_cr:
                    arrays[f"rep_cr_{k}"] = rep_cr[name]
                if name in rep_adj:
                    arrays[f"rep_adj_{k}"] = rep_adj[name]

            meta = {
                "title": args.title,
                "sample": case_sample,
                "n_intervals": n_intervals,
                "contigs": contigs,
                "panel_cr_samples": panel_cr_samples,
                "panel_adj_samples": panel_rc_samples,
                "panel_excludes_case": n_excluded > 0,
                "representatives": representatives,
                "panel_scale": PANEL_SCALE,
                "panel_missing": PANEL_MISSING,
                "has_gc": gc is not None,
                "calls": calls,
            }
            data_tags = "\n".join(f'<script type="application/octet-stream" id="gcnv-{name}">{pack(array)}</script>'
                                  for name, array in arrays.items())
            meta_json = json.dumps(meta, separators=(",", ":")).replace("</", "<\\/")
            # plotly.js is inlined last so that its source is not searched for the other placeholders.
            report = (HTML_TEMPLATE.replace("__TITLE__", html.escape(args.title))
                      .replace("__META__", meta_json)
                      .replace("__DATA__", data_tags)
                      .replace("__PLOTLY__", get_plotlyjs()))
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
            <div class="muted">Case sample only. Shaded bars mark passing calls. Click a chromosome name to open it in the Interval Browser.</div>
            <div id="facets"></div>
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
              <code>chr:position</code> or a whole contig name. Drag across a plot to zoom in; double-click to go back.</div>
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

        <script>__PLOTLY__</script>
        <script>
        (function () {
        "use strict";

        const META = JSON.parse(document.getElementById("gcnv-meta").textContent);
        const PANEL_COLOR = "rgb(37,99,235)";
        const CASE_COLOR = "#111827";
        const REP_COLOR = "#9ca3af";
        const CALL_RGB = { DEL: "220,38,38", DUP: "22,163,74" };
        const FILTERED_RGB = "107,114,128";
        const QUERY_RGB = "202,138,4";
        const REGION_YMAX = 7;
        const OVERVIEW_YMAX = 5;
        // Views holding more intervals than this merge neighbouring intervals into this many columns of panel points.
        const MAX_PANEL_COLUMNS = 1500;
        // Flanking context rules for call plots, ported from the original R report.
        // Widths below minCallWidth / minQueryWidth are treated as that width, so short calls and single positions get context.
        const CONTEXT = { minPoints: 10, minWidthFactor: 2.5, maxWidthFactor: 10, minCallWidth: 1000, minQueryWidth: 10000 };

        const PLOT_CONFIG = { displaylogo: false, responsive: true, doubleClick: false,
                              modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d", "resetScale2d"] };
        const BASE_LAYOUT = { paper_bgcolor: "#fff", plot_bgcolor: "#fff", showlegend: false,
                              font: { family: "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif", size: 11, color: "#374151" },
                              hoverlabel: { bgcolor: "rgba(17,24,39,.93)", bordercolor: "rgba(17,24,39,.93)", font: { color: "#fff", size: 12 }, align: "left" } };
        const AXIS = { showline: true, linecolor: "#374151", ticks: "outside", tickcolor: "#374151", zeroline: false, gridcolor: "#eef0f3" };
        const REF_LINES = [1, 2, 3].map(v => ({ type: "line", layer: "below", xref: "paper", x0: 0, x1: 1, yref: "y", y0: v, y1: v,
                                                line: { color: "rgba(55,65,81,.7)", width: 1, dash: "dash" } }));

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
          const reps = META.representatives.map((_, k) => ["rep_cr_" + k, "rep_adj_" + k]).flat();
          const [startDelta, lengths, caseCr, caseAdj, panelCr, panelAdj, gc, ...repArrays] = await Promise.all([
            loadArray("start_delta", Int32Array), loadArray("length", Int32Array),
            loadArray("case_cr", Float32Array), loadArray("case_adj", Float32Array),
            loadArray("panel_cr", Uint8Array), loadArray("panel_adj", Uint8Array), loadArray("gc", Float32Array),
            ...reps.map(name => loadArray(name, Float32Array))]);
          const n = META.n_intervals;
          const starts = new Float64Array(n), ends = new Float64Array(n), mids = new Float64Array(n);
          let pos = 0;
          for (let i = 0; i < n; i++) {
            pos += startDelta[i];
            starts[i] = pos; ends[i] = pos + lengths[i]; mids[i] = pos + lengths[i] / 2;
          }
          const repCr = repArrays.filter((_, k) => k % 2 === 0), repAdj = repArrays.filter((_, k) => k % 2 === 1);
          return { starts, ends, mids, caseCr, caseAdj, panelCr, panelAdj, gc, repCr, repAdj };
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
        // Views are at least 200 bp wide and start at position 1 or later.
        function clampView(ws, we) {
          if (we - ws < 200) { const c = (ws + we) / 2; ws = c - 100; we = c + 100; }
          return { ws: Math.max(1, ws), we };
        }
        function contextWindow(ci, start, end, mode) {
          if (mode === "auto") return autoWindow(ci, start, end, CONTEXT.minQueryWidth);
          const pad = (end - start + 1) * Number(mode);
          return clampView(start - pad, end + pad);
        }

        // ---------------------------------------------------------------- plotting
        // Plotly.react draws into an empty div or updates an existing plot. Event handlers are attached again after a purge.
        function plot(host, traces, layout) {
          Plotly.react(host, traces, Object.assign({}, BASE_LAYOUT, layout), PLOT_CONFIG);
          if (!host._listening) {
            host._listening = true;
            for (const [name, fn] of Object.entries(host._handlers || {})) host.on(name, fn);
          }
        }
        // Frees the plot's WebGL context; browsers only allow a few at a time.
        function releasePlot(host) {
          Plotly.purge(host);
          host._listening = false;
        }

        // Shaded rectangles for calls and the queried interval in view.
        function callShapes(ci, ws, we, focus) {
          const span = Math.max(we - ws, 1), shapes = [];
          // Narrow marks in wide views are left unlabelled so that labels do not pile up; the table below lists every call.
          const band = (start, end, fillcolor, line, text, textposition) => shapes.push({
            type: "rect", layer: "below", xref: "x", yref: "paper", y0: 0, y1: 1,
            x0: start, x1: Math.max(end, start + span * 0.002), fillcolor, line,
            label: (end - start) / span >= 0.03 ? { text, textposition, font: { size: 11, color: line.color || "#374151" } } : undefined });
          for (const c of META.calls) {
            if (c.contig !== ci.name || c.end < ws || c.start > we) continue;
            const rgb = c.pass ? (CALL_RGB[c.alt] || FILTERED_RGB) : FILTERED_RGB;
            const line = focus && focus.call === c ? { color: "rgba(17,24,39,.55)", width: 1.5 }
                       : c.pass ? { width: 0 } : { color: "rgba(" + rgb + ",.6)", width: 1, dash: "dash" };
            band(c.start, c.end, "rgba(" + rgb + "," + (c.pass ? 0.16 : 0.1) + ")", line,
                 c.alt + " CN=" + (c.cn == null ? "?" : c.cn) + (c.pass ? "" : " (" + c.filter + ")"), "top left");
          }
          if (focus && focus.query) {
            band(focus.query.start, focus.query.end, "rgba(" + QUERY_RGB + ",.1)",
                 { color: "rgba(" + QUERY_RGB + ",.9)", width: 1, dash: "dash" }, "entered interval", "bottom left");
          }
          return shapes;
        }

        // Each interval's panel values are sorted and quantized, so equal values are drawn once, with the opacity that
        // stacking that many points of opacity 0.2 would give.
        function panelTrace(track, lo, hi, x0, x1, yMax) {
          const P = track.P, pv = track.panel, sc = META.panel_scale;
          const n = hi - lo, cols = Math.min(n, MAX_PANEL_COLUMNS), binned = n > cols;
          const qMax = Math.min(META.panel_missing - 1, Math.floor(yMax * sc));
          const counts = new Uint32Array(cols * 256);
          for (let i = lo; i < hi; i++) {
            const col = binned ? Math.min(cols - 1, Math.max(0, Math.floor((D.mids[i] - x0) / (x1 - x0) * cols))) : i - lo;
            for (let j = 0, base = i * P; j < P; j++) {
              const q = pv[base + j];
              if (q > qMax) break;  // sorted, so values above the plot and missing values come last
              counts[col * 256 + q]++;
            }
          }
          const x = [], y = [], opacity = [];
          for (let col = 0; col < cols; col++) {
            for (let q = 0; q <= qMax; q++) {
              const k = counts[col * 256 + q];
              if (!k) continue;
              x.push(binned ? x0 + (col + 0.5) * (x1 - x0) / cols : D.mids[lo + col]);
              y.push(q / sc);
              opacity.push(1 - Math.pow(0.8, k));
            }
          }
          return { type: "scattergl", mode: "markers", x, y, hoverinfo: "skip",
                   marker: { color: PANEL_COLOR, opacity, size: binned || n > 3000 ? 2 : n > 800 ? 3 : 5 } };
        }

        function panelSummary(track, i) {
          const P = track.P, base = i * P, pv = track.panel, sc = META.panel_scale;
          const n = firstIndex(0, P, j => pv[base + j] === META.panel_missing);  // sorted, missing values last
          if (!n) return null;
          const v = j => pv[base + j] / sc, m = n >> 1;
          const median = n % 2 ? v(m) : (v(m - 1) + v(m)) / 2;
          return "median " + median.toFixed(2) + ", range " + v(0).toFixed(2) + "\u2013" + v(n - 1).toFixed(2) + " (n=" + n + ")";
        }
        function tooltipHtml(track, ci, i) {
          const f = v => Number.isFinite(v) ? v.toFixed(3) : "NA";
          const rows = ["<b>" + esc(ci.name) + ":" + fmt(D.starts[i]) + "-" + fmt(D.ends[i]) + "</b>",
                        "Case copy ratio: " + f(D.caseCr[i]), "Case adjusted counts: " + f(D.caseAdj[i])];
          if (D.gc) rows.push("GC content: " + f(D.gc[i]));
          const p = panelSummary(track, i);
          if (p) rows.push("Panel " + (track.key === "cr" ? "copy ratio" : "adjusted counts") + ": " + p);
          for (const r of track.reps) rows.push(esc(r.name) + ": " + f(r.vals[i]));
          return rows.join("<br>");
        }

        // A sample is drawn as points joined by a line; values above the plot are drawn as triangles just below its top.
        function sampleTrace(vals, lo, hi, yMax, color, hover) {
          const n = hi - lo, dot = n > 4000 ? 0 : n > 1000 ? 2.5 : 5;
          const x = [], y = [], size = [], symbol = [];
          for (let i = lo; i < hi; i++) {
            const v = vals[i], over = v > yMax;
            x.push(D.mids[i]);
            y.push(Number.isFinite(v) ? (over ? yMax - 0.2 : v) : null);
            size.push(over ? 9 : dot);
            symbol.push(over ? "triangle-up" : "circle");
          }
          return Object.assign({ type: "scattergl", mode: "lines+markers", x, y,
                                 line: { color, width: 1 }, marker: { color, size, symbol } }, hover);
        }
        // Only the case trace has a tooltip, which also lists the panel and representative sample values.
        function caseTrace(track, ci, lo, hi, yMax) {
          const text = [];
          for (let i = lo; i < hi; i++) text.push(tooltipHtml(track, ci, i));
          return sampleTrace(track.caseVals, lo, hi, yMax, CASE_COLOR, { text, hovertemplate: "%{text}<extra></extra>" });
        }

        // Scatter of panel values (blue) with representative panel samples (grey) and the case (black) over a genomic window.
        function renderRegion(host, track, ci, view, shapes) {
          const yMax = REGION_YMAX, pad = Math.max(view.we - view.ws, 1) * 0.02, x0 = view.ws - pad, x1 = view.we + pad;
          const [lo, hi] = indexRange(ci, x0, x1);
          const traces = track.reps.map(r => sampleTrace(r.vals, lo, hi, yMax, REP_COLOR, { hoverinfo: "skip" }));
          if (track.P > 0 && hi > lo) traces.unshift(panelTrace(track, lo, hi, x0, x1, yMax));
          traces.push(caseTrace(track, ci, lo, hi, yMax));
          plot(host, traces, {
            height: 300, margin: { l: 58, r: 12, t: 16, b: 44 }, dragmode: "zoom", hovermode: "x", hoverdistance: -1,
            xaxis: Object.assign({}, AXIS, { range: [x0, x1], tickformat: ",d", showgrid: false, title: { text: "Position on " + ci.name },
                                             showspikes: true, spikemode: "across", spikesnap: "data", spikedash: "solid",
                                             spikethickness: 1, spikecolor: "rgba(17,24,39,.35)" }),
            yaxis: Object.assign({}, AXIS, { range: [0, yMax], dtick: 1, fixedrange: true, title: { text: track.yLabel } }),
            shapes: shapes.concat(REF_LINES),
            annotations: hi > lo ? [] : [{ text: "No gCNV intervals in this region", xref: "paper", yref: "paper", x: 0.5, y: 0.5,
                                           showarrow: false, font: { size: 13, color: "#6b7280" } }],
          });
        }

        // Small multiples of the case copy ratio across each contig, drawn as one figure so that they share a WebGL context.
        function renderFacets(host) {
          const W = host.clientWidth, cols = Math.max(1, Math.floor(W / 230)), rows = Math.ceil(META.contigs.length / cols);
          const ROW = 150, LABEL = 18, GAP = 8, AXIS_W = 24, H = rows * ROW;
          const traces = [], shapes = [], annotations = [];
          const layout = { height: H, margin: { l: 0, r: 0, t: 0, b: 0 }, hovermode: "closest", dragmode: false, shapes, annotations };
          META.contigs.forEach((ci, k) => {
            const s = k ? String(k + 1) : "", r = Math.floor(k / cols), c = k % cols;
            const a = ci.offset, b = ci.offset + ci.count, x0 = D.starts[a], x1 = Math.max(D.ends[b - 1], x0 + 1);
            const xDomain = [(c * W / cols + AXIS_W) / W, ((c + 1) * W / cols - GAP) / W];
            const yDomain = [1 - ((r + 1) * ROW - GAP) / H, 1 - (r * ROW + LABEL) / H];
            layout["xaxis" + s] = { domain: xDomain, anchor: "y" + s, range: [x0, x1], fixedrange: true, showticklabels: false,
                                    showgrid: false, zeroline: false, showline: true, linecolor: "#374151" };
            layout["yaxis" + s] = Object.assign({}, AXIS, { domain: yDomain, anchor: "x" + s, range: [0, OVERVIEW_YMAX], dtick: 1,
                                                            fixedrange: true, tickfont: { size: 9 }, ticklen: 3 });
            traces.push({ type: "scattergl", mode: "markers", xaxis: "x" + s, yaxis: "y" + s,
                          x: D.mids.subarray(a, b), y: D.caseCr.subarray(a, b), marker: { color: CASE_COLOR, size: 2, opacity: 0.2 },
                          hovertemplate: esc(ci.name) + ":%{x:,.0f}<br>Copy ratio: %{y:.3f}<extra></extra>" });
            for (const call of META.calls) {
              if (!call.pass || call.contig !== ci.name) continue;
              shapes.push({ type: "rect", layer: "below", xref: "x" + s, yref: "y" + s + " domain", y0: 0, y1: 1, line: { width: 0 },
                            x0: call.start, x1: Math.max(call.end, call.start + (x1 - x0) * 0.005),
                            fillcolor: "rgba(" + (CALL_RGB[call.alt] || FILTERED_RGB) + ",.3)" });
            }
            annotations.push({ text: esc(ci.name), xref: "x" + s + " domain", yref: "y" + s + " domain", x: 0.5, y: 1,
                               yanchor: "bottom", showarrow: false, captureevents: true, bgcolor: "#e5e7eb",
                               width: (xDomain[1] - xDomain[0]) * W, height: LABEL - 4, font: { size: 12, color: "#1f2937" } });
          });
          plot(host, traces, layout);
        }

        function renderGc(host) {
          plot(host, [{ type: "scattergl", mode: "markers", x: D.gc, y: D.caseAdj, hoverinfo: "skip",
                        marker: { color: CASE_COLOR, size: 2, opacity: 0.2 } }], {
            height: 380, margin: { l: 58, r: 12, t: 10, b: 44 }, hovermode: false,
            xaxis: Object.assign({}, AXIS, { range: [0, 1], dtick: 0.2, tickformat: ".1f", title: { text: "GC content" } }),
            yaxis: Object.assign({}, AXIS, { range: [0, OVERVIEW_YMAX], dtick: 1, title: { text: "Adjusted read counts" } }),
          });
        }

        function legendHtml(focusText) {
          const nCr = META.panel_cr_samples.length, nAdj = META.panel_adj_samples.length;
          const panel = (nCr === nAdj ? nCr + " samples" : nCr + " copy ratio / " + nAdj + " read count samples") +
                        (META.panel_excludes_case ? ", excluding the case" : "");
          const reps = META.representatives.map(esc).join(", ");
          return '<div class="legend">' +
            '<span><i class="sw" style="background:' + CASE_COLOR + '"></i>Case (' + esc(META.sample) + ')</span>' +
            (reps ? '<span><i class="sw" style="background:' + REP_COLOR + '"></i>Representative panel ' +
                    (META.representatives.length > 1 ? "samples" : "sample") + ' (' + reps + ')</span>' : "") +
            '<span><i class="sw" style="background:rgba(37,99,235,.45)"></i>Panel (' + panel + ')</span>' +
            '<span><i class="sw rect" style="background:rgba(' + CALL_RGB.DEL + ',.3)"></i>Passing DEL</span>' +
            '<span><i class="sw rect" style="background:rgba(' + CALL_RGB.DUP + ',.3)"></i>Passing DUP</span>' +
            '<span><i class="sw rect" style="background:rgba(' + FILTERED_RGB + ',.25);border:1px dashed #6b7280"></i>Filtered call</span>' +
            (focusText || "") + "</div>";
        }

        // A pair of plots (copy ratio and read counts) for one window. Zooming either plot calls onZoom with the new window
        // and double-clicking calls onReset; both are expected to redraw the pair.
        function makeRegionPair(parent, onZoom, onReset) {
          const onRelayout = ev => {
            const r = ev["xaxis.range"] || ("xaxis.range[0]" in ev ? [ev["xaxis.range[0]"], ev["xaxis.range[1]"]] : null);
            if (r) onZoom(Math.min(r[0], r[1]), Math.max(r[0], r[1]));
          };
          const pair = el("div", { className: "pair" });
          const hosts = TRACKS.map(t => {
            const cell = el("div"), host = el("div");
            host._handlers = { plotly_relayout: onRelayout };
            host.addEventListener("dblclick", onReset);
            cell.append(el("h4", {}, esc(t.title)), host);
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
        let overview = null;
        function renderOverview() {
          if (!overview) {
            overview = { facets: document.getElementById("facets"), gc: null };
            const openContig = k => openBrowser(META.contigs[k].name, "0");
            overview.facets._handlers = { plotly_clickannotation: e => openContig(e.index),
                                          plotly_click: e => openContig(e.points[0].curveNumber) };
            const gcHost = document.getElementById("gc-host");
            if (D.gc) {
              overview.gc = el("div");
              overview.gc.style.maxWidth = "760px";
              gcHost.append(overview.gc);
            } else {
              gcHost.innerHTML = "<div class='muted'>GC content is not available: gCNV was run without annotated intervals.</div>";
            }
          }
          renderFacets(overview.facets);
          if (overview.gc) renderGc(overview.gc);
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
          // Only events near the screen are plotted, to stay within the browser's limit on WebGL contexts.
          const observer = new IntersectionObserver(entries => {
            for (const entry of entries) {
              const box = entry.target;
              if (entry.isIntersecting && !box._rendered) { box._rendered = true; box._render(); }
              else if (!entry.isIntersecting && box._rendered) { box._rendered = false; box._hosts.forEach(releasePlot); }
            }
          }, { rootMargin: "300px 0px" });
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
            const initial = autoWindow(ci, c.start, c.end, CONTEXT.minCallWidth);
            let view = initial;
            box._render = () => {
              const shapes = callShapes(ci, view.ws, view.we, { call: c });
              box._hosts.forEach((h, k) => renderRegion(h, TRACKS[k], ci, view, shapes));
            };
            box._hosts = makeRegionPair(box, (ws, we) => { view = clampView(ws, we); box._render(); },
                                        () => { view = initial; box._render(); });
            observer.observe(box);
          });
        }
        function rerenderEvents() {
          document.querySelectorAll("#events-list .event").forEach(b => { if (b._rendered) b._render(); });
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
          const s = text.trim().replace(/[–—]/g, "-");
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
          browser.view = clampView(ws, we);
          renderBrowser();
        }
        function renderBrowser() {
          const q = browser.query, v = browser.view;
          document.querySelectorAll("#nav-controls button").forEach(b => { b.disabled = !q; });
          if (!q) return;
          const container = document.getElementById("browser-plots");
          if (!browser.hosts) {
            container.innerHTML = legendHtml('<span><i class="sw rect" style="background:rgba(' + QUERY_RGB + ',.2);border:1px dashed rgb(' + QUERY_RGB + ')"></i>Entered interval</span>');
            browser.hosts = makeRegionPair(container, setView, resetView);
          }
          const whole = q.start === D.starts[q.ci.offset] && q.end === D.ends[q.ci.offset + q.ci.count - 1];
          const shapes = callShapes(q.ci, v.ws, v.we, whole ? {} : { query: q });
          browser.hosts.forEach((h, k) => renderRegion(h, TRACKS[k], q.ci, v, shapes));
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
          if (name !== activeTab) document.querySelectorAll("#tab-" + activeTab + " .js-plotly-plot").forEach(releasePlot);
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
          // A representative sample missing from one of the panels is left out of that track.
          const reps = arrays => META.representatives.map((name, k) => ({ name, vals: arrays[k] })).filter(r => r.vals);
          TRACKS = [
            { key: "cr", title: "Denoised Copy Ratio", yLabel: "Denoised linear copy ratio", caseVals: D.caseCr, panel: D.panelCr, P: META.panel_cr_samples.length, reps: reps(D.repCr) },
            { key: "adj", title: "Adjusted Read Counts", yLabel: "Adjusted read counts", caseVals: D.caseAdj, panel: D.panelAdj, P: META.panel_adj_samples.length, reps: reps(D.repAdj) },
          ];
          const nPass = META.calls.filter(c => c.pass).length;
          document.getElementById("summary").textContent = "Sample " + META.sample + " · " + nPass + " passing and " +
            (META.calls.length - nPass) + " filtered CNV calls · " + fmt(META.n_intervals) + " gCNV intervals";
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
          // Plotly resizes plots to their containers itself; the overview is redrawn because its number of columns can change.
          let resizeTimer = null;
          window.addEventListener("resize", () => {
            clearTimeout(resizeTimer);
            resizeTimer = setTimeout(() => { if (activeTab === "overview") renderOverview(); }, 150);
          });

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
            --representative-samples-list ~{write_lines(representative_sample_ids)} \
            --title "~{output_prefix}" \
            --output ~{output_prefix}_cnv_event_report.html
    >>>

    runtime {
        # Built from Utilities/Dockers/Python-h5py/Dockerfile, which includes plotly.
        docker: "us.gcr.io/broad-dsde-methods/python-h5py:plotly-6.3.0"
        disks: "local-disk 100 HDD"
        memory: mem_gb + " GB"
    }

    output {
        File cnv_event_report = "~{output_prefix}_cnv_event_report.html"
    }
}
