import pandas as pd
import numpy as np
import pysam
from datetime import datetime

class BGEScorer():
    def __init__(self, ref_dict_path, prs_weights_path, output_basename, score_haploid_as_diploid):
        """
        Initialize the BGEScorer object.

        Args:
            ref_dict_path (str): The file path to the reference dictionary.
            prs_weights_path (str): The file path to the PRS weights.
            output_basename (str): The base name for the output files. Sites scored will be written during scoring, final scores and any_source_any_sample.sites_scored will be written when calling write_output.
            score_haploid_as_diploid (bool): Whether to score haploid genotypes as diploid.
        """
        self.ref_dict = self._read_ref_dict(ref_dict_path)
        self.prs_weights = self._read_prs_weights(prs_weights_path)
        self.output_basename = output_basename
        self.sample_names = None
        self.score_haploid_as_diploid = score_haploid_as_diploid

        self.gvcf_sample_score = None
        self.gvcf_sites_scored = None
        self.gvcf_sites_scored_set = None
        self.gvcf_low_quality_sites = None

        self.vcf_sample_score = None
        self.vcf_sites_scored = None
        self.vcf_num_sites_not_found = None

        self.any_source_any_sample_sites_scored = None

    def _read_ref_dict(self, ref_dict_path):
        with open(ref_dict_path) as ref_dict:
            return [line.split('\t')[1].replace('SN:', '') for line in ref_dict if line.startswith('@SQ')]

    # Required TSV format (including header line):
    # contig position ref alt effect_allele weight
    def _read_prs_weights(self, prs_weights_path):
        prs_weights = pd.read_table(prs_weights_path,
            dtype={'contig': str, 'position': np.int32, 'ref': str, 'alt': str, 'effect_allele': str, 'weight': np.float64})

        return prs_weights.sort_values(by='position').sort_values(kind='stable', by='contig', key=lambda contig: contig.map(self.ref_dict.index))
    
    def _compare_record_and_weight(self, record, weight):
        if record.contig != weight.contig:
            return self.ref_dict.index(record.contig) - self.ref_dict.index(weight.contig)
        if record.pos <= weight.position <= record.stop:
            return 0
        return record.pos - weight.position

    def _check_biallelic(self, record):
        if record is not None and len(record.alts) > 1:
            raise RuntimeError('Only biallelic variants are supported: ' + str(record))
    
    def _gvcf_score_site(self, record, weight, site_gq_threshold):
        samples_scored = []
        for sample_name in self.sample_names:
            if record.samples[sample_name]['GQ'] < site_gq_threshold or None in record.samples[sample_name]['GT']:
                self.gvcf_low_quality_sites[sample_name].append((weight.contig, weight.position, weight.ref, weight.alt))
                continue
            
            if weight.effect_allele == weight.ref:
                effect_allele_index = 0
            else:
                effect_allele_index = record.alts.index(weight.effect_allele) + 1 if weight.effect_allele in record.alts else None

            site_score = 0 if effect_allele_index is None else record.samples[sample_name]['GT'].count(effect_allele_index) * weight.weight

            if self.score_haploid_as_diploid and len(record.samples[sample_name]['GT']) == 1:
                site_score *= 2

            self.gvcf_sites_scored[sample_name].append((weight.contig, weight.position, weight.ref, weight.alt))
            self.gvcf_sample_score[sample_name] += site_score
            samples_scored.append(sample_name)
        return samples_scored
    
    def _vcf_score_dosage(self, record, weight):
        samples_scored = []
        for sample_name in self.sample_names:
            if (weight.contig, weight.position, weight.ref, weight.alt) in self.gvcf_sites_scored_set[sample_name]:
                continue

            if 'DS' not in record.samples[sample_name]:
                raise RuntimeError(f'VCF file does not contain dosage information for sample {sample_name} at site {weight.contig}:{weight.position}:{weight.ref}:{weight.alt}')

            dosage = record.samples[sample_name]['DS'] if record.alts[0] == weight.effect_allele else 2 - record.samples[sample_name]['DS']
            site_score = dosage * weight.weight

            self.vcf_sites_scored[sample_name].append((weight.contig, weight.position, weight.ref, weight.alt))
            self.vcf_sample_score[sample_name] += site_score
            samples_scored.append(sample_name)
        return samples_scored
    
    def _print_wes_gvcf_metrics(self):
        num_sites_scored = {sample_name: len(self.gvcf_sites_scored[sample_name]) for sample_name in self.sample_names}
        num_sites_scored_min_max = min(num_sites_scored.values()), max(num_sites_scored.values())
        num_low_quality_sites = {sample_name: len(self.gvcf_low_quality_sites[sample_name]) for sample_name in self.sample_names}
        num_low_quality_sites_min_max = min(num_low_quality_sites.values()), max(num_low_quality_sites.values())

        print('  Metrics:')
        print(f'    Sites scored: Min: {num_sites_scored_min_max[0]} Max: {num_sites_scored_min_max[1]}')
        print(f'    Low quality sites: Min: {num_low_quality_sites_min_max[0]} Max: {num_low_quality_sites_min_max[1]}')

    def _print_wgs_vcf_metrics(self):
        num_sites_scored = {sample_name: len(self.vcf_sites_scored[sample_name]) for sample_name in self.sample_names}
        num_sites_scored_min_max = min(num_sites_scored.values()), max(num_sites_scored.values())

        print('  Metrics:')
        print(f'    Sites scored: Min: {num_sites_scored_min_max[0]} Max: {num_sites_scored_min_max[1]}')
        print(f'    Sites not found: {self.vcf_num_sites_not_found}')

    def _print_wes_and_wgs_metrics(self):
        total_sites_scored = {sample_name: len(self.gvcf_sites_scored[sample_name]) + len(self.vcf_sites_scored[sample_name]) for sample_name in self.sample_names}
        sites_scored_min_max = min(total_sites_scored.values()), max(total_sites_scored.values())

        print('WES GVCF + WGS VCF Scoring:')
        print(f'    Total sites scored: Min: {sites_scored_min_max[0]} Max: {sites_scored_min_max[1]}')

    def _process_weight_wes(self, weight, gvcf, site_gq_threshold, out_sites_scored):
        site_records = gvcf.fetch(weight.contig, weight.position - 1, weight.position)
        record = next(site_records, None)
        # Skip all records that are before the current record
        while record is not None and self._compare_record_and_weight(record, weight) < 0:
            record = next(site_records, None)

        if record is None or self._compare_record_and_weight(record, weight) != 0:
            samples_scored = []
        else:
            samples_scored = self._gvcf_score_site(record, weight, site_gq_threshold)
                
        out_sites_scored.write(f'{weight.contig}:{weight.position}:{weight.ref}:{weight.alt}\t{",".join(samples_scored)}\n')
        if len(samples_scored) > 0:
            self.any_source_any_sample_sites_scored.add((weight.contig, weight.position, weight.ref, weight.alt))

    def score_wes_gvcf(self, gvcf_path, sample_names=None, site_gq_threshold=30):
        """
        Score variants in a Whole Exome Sequencing (WES) GVCF file if they meet the genotype quality threshold. No-calls and half-calls are treated as low quality.

        Args:
            gvcf_path (str): The file path to the WES GVCF file.
            sample_names (list, optional): A list of sample names to score. If None, all samples in the GVCF file will be scored. Default is None.
            variant_gq_threshold (int, optional): The genotype quality threshold for variants. Variants with GQ below this threshold will be considered low quality. Default is 30.
            ref_block_gq_threshold (int, optional): The genotype quality threshold for reference blocks. Reference blocks with GQ below this threshold will be considered low quality. Default is 30.
        """
        if self.gvcf_sample_score is not None or self.vcf_sample_score is not None:
            raise RuntimeError('Either WES GVCF or WGS VCF scoring has already been performed. If you want to re-score, create a new BGEScorer object.')

        print('WES GVCF Scoring:')
        with pysam.VariantFile(gvcf_path, "r", ignore_truncation=True) as gvcf:
            if sample_names is None:
                self.sample_names = list(gvcf.header.samples)
            else:
                if any(sample_name not in gvcf.header.samples for sample_name in sample_names):
                    raise RuntimeError('One or more sample_names are not present in the GVCF file.')
                self.sample_names = sample_names

            print(f'  Scoring {len(self.sample_names)} samples...')

            self.gvcf_sample_score = {sample_name: 0 for sample_name in self.sample_names}
            self.gvcf_sites_scored = {sample_name: [] for sample_name in self.sample_names}
            self.gvcf_low_quality_sites = {sample_name: [] for sample_name in self.sample_names}
            self.any_source_any_sample_sites_scored = set()

            start_time = datetime.now()
            step_time = start_time

            with open(self.output_basename + '.exome_gvcf.sites_scored', 'w') as out_sites_scored:
                out_sites_scored.write('site\tsamples_scored\n')
                for weight in self.prs_weights.itertuples():
                    i_weight = weight.Index
                    if i_weight % 100000 == 0:
                        print(f'Scored {i_weight:,} sites. Current locus: {weight.contig}:{weight.position}. Time elapsed: {(datetime.now() - start_time).total_seconds():,.0f}s. Time since last step: {(datetime.now() - step_time).total_seconds():,.0f}s')
                        step_time = datetime.now()
                    self._process_weight_wes(weight, gvcf, site_gq_threshold, out_sites_scored)

        # Save the scored sites as a set for faster lookup
        self.gvcf_sites_scored_set = {key: set(value) for key, value in self.gvcf_sites_scored.items()}
        self._print_wes_gvcf_metrics()

    def _process_weight_wgs(self, weight, vcf, out_sites_scored):

        site_records = vcf.fetch(weight.contig, weight.position - 1, weight.position)

        # Get first record
        record = next(site_records, None)
        self._check_biallelic(record)

        # Find the record that matches the weight
        while record is not None and (self._compare_record_and_weight(record, weight) != 0 or record.ref != weight.ref or weight.alt != record.alts[0]):
            record = next(site_records, None)
            self._check_biallelic(record)
        
        # Skip the weight if the record was not found
        if record is None:
            self.vcf_num_sites_not_found += 1
            print('  Site not found in VCF:', weight.contig, weight.position, weight.ref, weight.alt)
            samples_scored = []
        else:
            samples_scored = self._vcf_score_dosage(record, weight)

        out_sites_scored.write(f'{weight.contig}:{weight.position}:{weight.ref}:{weight.alt}\t{",".join(samples_scored)}\n')
        if len(samples_scored) > 0:
            self.any_source_any_sample_sites_scored.add((weight.contig, weight.position, weight.ref, weight.alt))

    def score_wgs_vcf(self, wgs_vcf_path, sample_names=None, allow_wgs_vcf_only=False):
        """
        Score variants in a Whole Genome Sequencing (WGS) VCF file.

        Args:
            wgs_vcf_path (str): The file path to the WGS VCF file.
            sample_names (list, optional): A list of sample names to score. If None, all samples in the VCF file will be scored. Must match the sample_names from WES GVCF scoring, if performed. Default is None.
            allow_wgs_vcf_only (bool, optional): Whether to allow scoring of the WGS VCF file only, without prior scoring of the WES GVCF file. If False and WES GVCF scoring was not performed, a RuntimeError will be raised. Default is False.
        """
        print('WGS VCF Scoring:')
        if not allow_wgs_vcf_only and self.gvcf_sites_scored_set is None:
            raise RuntimeError('WES GVCF scoring must be performed before WGS VCF scoring. If you want to score the WGS VCF only, set allow_wgs_vcf_only to True')

        with pysam.VariantFile(wgs_vcf_path, "r") as vcf:
            if sample_names is None:
                sample_names = list(vcf.header.samples)
            
            # Check if the sample names match the ones from the WES GVCF scoring, or if they are None (i.e. WES GVCF scoring was not performed), set them below
            if self.sample_names is not None and sample_names != self.sample_names:
                raise RuntimeError('Sample names in WES GVCF and WGS VCF do not match. If you want to score only a subset of the samples, provide the same sample_names argument to both score_wes_gvcf and score_wgs_vcf.\nWES GVCF sample names:' +str(self.sample_names) + '\nWGS VCF sample names:' + str(sample_names))
            self.sample_names = sample_names

            print(f'  Scoring {len(self.sample_names)} samples...')

            if self.gvcf_sites_scored_set is None:
                self.gvcf_sites_scored_set = {sample_name: set() for sample_name in self.sample_names}

            self.vcf_sample_score = {sample_name: 0 for sample_name in self.sample_names}
            self.vcf_sites_scored = {sample_name: [] for sample_name in self.sample_names}
            self.vcf_num_sites_not_found = 0
            if self.any_source_any_sample_sites_scored is None:
                self.any_source_any_sample_sites_scored = set()

            start_time = datetime.now()
            step_time = start_time

            with open(self.output_basename + '.imputed_wgs_vcf.sites_scored', 'w') as out_sites_scored:
                out_sites_scored.write('site\tsamples_scored\n')
                for weight in self.prs_weights.itertuples():
                    i_weight = weight.Index
                    if i_weight % 100000 == 0:
                        print(f'Scored {i_weight:,} sites. Current locus: {weight.contig}:{weight.position}. Time elapsed: {(datetime.now() - start_time).total_seconds():,.0f}s. Time since last step: {(datetime.now() - step_time).total_seconds():,.0f}s')
                        step_time = datetime.now()
                    self._process_weight_wgs(weight, vcf, out_sites_scored)
        
        self._print_wgs_vcf_metrics()
        if self.gvcf_sites_scored is not None:
            self._print_wes_and_wgs_metrics()

    def write_output(self, allow_single_source_scoring=False):
        """
        Write the scores to plink-like output files.
        """
        def write_gvcf_or_vcf_output(is_gvcf):
            score = self.gvcf_sample_score if is_gvcf else self.vcf_sample_score
            score_filename_suffix = '.exome_gvcf.score' if is_gvcf else '.imputed_wgs_vcf.score'

            with open(self.output_basename + score_filename_suffix, 'w') as out_score:
                out_score.write('#IID\tSCORE1_SUM\n')
                for sample_name in self.sample_names:
                    out_score.write(f'{sample_name}\t{score[sample_name]}\n')
        
        # GVCF
        if self.gvcf_sites_scored_set is None:
            if not allow_single_source_scoring:
                raise RuntimeError('No GVCF scoring performed. If you want to output scores from a single source, set allow_single_source_scoring to True.')
            else:
                print('No GVCF scoring performed. Skipping GVCF output.')
        else:
            write_gvcf_or_vcf_output(True)

        # VCF
        if self.vcf_sites_scored is None:
            if not allow_single_source_scoring:
                raise RuntimeError('No VCF scoring performed. If you want to output scores from a single source, set allow_single_source_scoring to True.')
            else:
                print('No VCF scoring performed. Skipping VCF output.')
        else:
            write_gvcf_or_vcf_output(False)
        
        # All sites (in plink-compatible format)
        with open(f'{self.output_basename}.any_source_any_sample.sites_scored', 'w') as out_sites:
            for site in sorted(self.any_source_any_sample_sites_scored, key=lambda site: (self.ref_dict.index(site[0]), site[1], site[2], site[3])):
                out_sites.write(f'{site[0]}:{site[1]}:{site[2]}:{site[3]}\n')

        # Sum
        with open(f'{self.output_basename}.score', 'w') as out_combined_score:
            out_combined_score.write('#IID\tSCORE1_SUM\n')
            for sample_name in self.sample_names:
                # We already checked above if single source scoring is allowed, so we don't need to check again
                vcf_score = 0 if self.vcf_sample_score is None else self.vcf_sample_score[sample_name]
                gvcf_score = 0 if self.gvcf_sample_score is None else self.gvcf_sample_score[sample_name]
                out_combined_score.write(f'{sample_name}\t{vcf_score + gvcf_score}\n')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Score BGE')
    parser.add_argument('--weights', type=str, help='Path to weights file', required=True)
    parser.add_argument('--gvcf', type=str, help='Path to WES GVCF', required=True)
    parser.add_argument('--vcf', type=str, help='Path to imputed WGS VCF', required=True)
    parser.add_argument('--ref-dict', type=str, help='Path to reference dict file', required=True)
    parser.add_argument('--basename', type=str, help='Path to reference dict file', required=True)
    parser.add_argument('--sample-names', type=str, nargs='+', help='Sample names to score', required=False, default=None)
    parser.add_argument('--score-haploid-as-diploid', action='store_true', help='Always score haploid GTs (such as on chrX) as diploid', required=False, default=False)
    args = parser.parse_args()

    bge_scorer = BGEScorer(args.ref_dict, args.weights, args.basename, args.score_haploid_as_diploid)
    bge_scorer.score_wes_gvcf(args.gvcf, sample_names=args.sample_names)
    bge_scorer.score_wgs_vcf(args.vcf, sample_names=args.sample_names)
    bge_scorer.write_output()
