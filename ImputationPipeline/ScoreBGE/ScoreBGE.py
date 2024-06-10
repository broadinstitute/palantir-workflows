import pandas as pd
import numpy as np
import pysam

class BGEScorer():
    def __init__(self, ref_dict_path, prs_weights_path):
        """
        Initialize the BGEScorer object.

        Args:
            ref_dict_path (str): The file path to the reference dictionary.
            prs_weights_path (str): The file path to the PRS weights.
        """
        self.ref_dict = self._read_ref_dict(ref_dict_path)
        self.prs_weights = self._read_prs_weights(prs_weights_path)
        self.sample_names = None

        self.gvcf_sample_score = None
        self.gvcf_sites_scored = None
        self.gvcf_sites_scored_set = None
        self.gvcf_low_quality_sites = None

        self.vcf_sample_score = None
        self.vcf_sites_scored = None
        self.vcf_num_sites_not_found = None

    def _read_ref_dict(self, ref_dict_path):
        with open(ref_dict_path) as ref_dict:
            return [line.split('\t')[1].replace('SN:', '') for line in ref_dict if line.startswith('@SQ')]

    def _read_prs_weights(self, prs_weights_path):
        prs_weights = pd.read_table(prs_weights_path, dtype={'effect_allele': str, 'weight': np.float64, 'chr': str, 'bp': np.int32, 'ref': str, 'alt': str})

        prs_weights[['contig', 'position', 'ref', 'alt']] = prs_weights['CHR:BP:REF:ALT'].str.split(':', expand=True)
        prs_weights = prs_weights.drop(['CHR:BP:REF:ALT'], axis=1)
        prs_weights['locus'] = prs_weights['contig'] + ':' + prs_weights['position']

        prs_weights = prs_weights.astype({'effect_allele': str, 'weight': np.float64, 'contig': str, 'position': np.int32, 'ref': str, 'alt': str})

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
        for sample_name in self.sample_names:
            if record.samples[sample_name]['GQ'] < site_gq_threshold:
                self.gvcf_low_quality_sites[sample_name].append((weight.locus, weight.ref, weight.alt))
                return
            
            #effect_allele_index = None if weight.effect_allele not in record.alts
            if weight.effect_allele == weight.ref:
                effect_allele_index = 0
            else:
                effect_allele_index = record.alts.index(weight.effect_allele) if weight.effect_allele in record.alts else None

            site_score = 0 if effect_allele_index is None else record.samples[sample_name]['GT'].count(effect_allele_index) * weight.weight

            self.gvcf_sites_scored[sample_name].append((weight.locus, weight.ref, weight.alt))
            self.gvcf_sample_score[sample_name] += site_score
    
    def _vcf_score_dosage(self, record, weight):
        for sample_name in self.sample_names:
            if (weight.locus, weight.ref, weight.alt) in self.gvcf_sites_scored_set[sample_name]:
                continue
            dosage = record.samples[sample_name]['DS'] if record.alts[0] == weight.effect_allele else 2 - record.samples[sample_name]['DS']
            site_score = dosage * weight.weight

            self.vcf_sites_scored[sample_name].append((weight.locus, weight.ref, weight.alt))
            self.vcf_sample_score[sample_name] += site_score
    
    def _print_wes_gvcf_metrics(self):
        num_sites_scored = {sample_name: len(self.gvcf_sites_scored[sample_name]) for sample_name in self.sample_names}
        num_sites_scored_min_max = min(num_sites_scored.values()), max(num_sites_scored.values())
        num_low_quality_sites = {sample_name: len(self.gvcf_low_quality_sites[sample_name]) for sample_name in self.sample_names}
        num_low_quality_sites_min_max = min(num_low_quality_sites.values()), max(num_low_quality_sites.values())

        print(f'  Metrics:')
        print(f'    Sites scored: Min: {num_sites_scored_min_max[0]} Max: {num_sites_scored_min_max[1]}')
        print(f'    Low quality sites: Min: {num_low_quality_sites_min_max[0]} Max: {num_low_quality_sites_min_max[1]}')

    def _print_wgs_vcf_metrics(self):
        num_sites_scored = {sample_name: len(self.vcf_sites_scored[sample_name]) for sample_name in self.sample_names}
        num_sites_scored_min_max = min(num_sites_scored.values()), max(num_sites_scored.values())

        print(f'  Metrics:')
        print(f'    Sites scored: Min: {num_sites_scored_min_max[0]} Max: {num_sites_scored_min_max[1]}')
        print(f'    Sites not found: {self.vcf_num_sites_not_found}')

    def _print_wes_and_wgs_metrics(self):
        total_sites_scored = {sample_name: len(self.gvcf_sites_scored[sample_name]) + len(self.vcf_sites_scored[sample_name]) for sample_name in self.sample_names}
        sites_scored_min_max = min(total_sites_scored.values()), max(total_sites_scored.values())

        print(f'WES GVCF + WGS VCF Scoring:')
        print(f'    Total sites scored: Min: {sites_scored_min_max[0]} Max: {sites_scored_min_max[1]}')

    def score_wes_gvcf(self, gvcf_path, sample_names=None, site_gq_threshold=30):
        """
        Score variants in a Whole Exome Sequencing (WES) GVCF file.

        Args:
            gvcf_path (str): The file path to the WES GVCF file.
            sample_names (list, optional): A list of sample names to score. If None, all samples in the GVCF file will be scored. Default is None.
            variant_gq_threshold (int, optional): The genotype quality threshold for variants. Variants with GQ below this threshold will be considered low quality. Default is 30.
            ref_block_gq_threshold (int, optional): The genotype quality threshold for reference blocks. Reference blocks with GQ below this threshold will be considered low quality. Default is 30.
        """
        print('WES GVCF Scoring:')
        gvcf = pysam.VariantFile(gvcf_path, "r")

        if sample_names is None:
            self.sample_names = list(gvcf.header.samples)
        else:
            if any(sample_name not in gvcf.header.samples for sample_name in sample_names):
                raise RuntimeError('One or more sample_names are not present in the GVCF file.')
            self.sample_names = sample_names

        self.gvcf_sample_score = {sample_name: 0 for sample_name in self.sample_names}
        self.gvcf_sites_scored = {sample_name: [] for sample_name in self.sample_names}
        self.gvcf_low_quality_sites = {sample_name: [] for sample_name in self.sample_names}

        current_record_iter = gvcf.fetch()
        current_record = next(current_record_iter)

        for _, weight in self.prs_weights.iterrows():
            # Skip all weights that are before the current record
            while self._compare_record_and_weight(current_record, weight) < 0:
                try:
                    current_record = next(current_record_iter)
                except StopIteration:
                    break
            
            if self._compare_record_and_weight(current_record, weight) == 0:
                self._gvcf_score_site(current_record, weight, site_gq_threshold)
        
        # Save the scored sites as a set for faster lookup
        self.gvcf_sites_scored_set = {key: set(value) for key, value in self.gvcf_sites_scored.items()}
        self._print_wes_gvcf_metrics()

    
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

        vcf = pysam.VariantFile(wgs_vcf_path, "r")

        if sample_names is None:
            sample_names = list(vcf.header.samples)
        
        # Check if the sample names match the ones from the WES GVCF scoring, or if they are None (i.e. WES GVCF scoring was not performed), set them below
        if self.sample_names is not None and sample_names != self.sample_names:
            raise RuntimeError('Sample names in WES GVCF and WGS VCF do not match. If you want to score only a subset of the samples, provide the same sample_names argument to both score_wes_gvcf and score_wgs_vcf.')
        self.sample_names = sample_names

        if self.gvcf_sites_scored_set is None:
            self.gvcf_sites_scored_set = {sample_name: set() for sample_name in self.sample_names}

        self.vcf_sample_score = {sample_name: 0 for sample_name in self.sample_names}
        self.vcf_sites_scored = {sample_name: [] for sample_name in self.sample_names}
        self.vcf_num_sites_not_found = 0

        for _, weight in self.prs_weights.iterrows():
            site_records = vcf.fetch(weight.contig, weight.position - 1, weight.position)

            # Get first record
            record = next(site_records, None)
            self._check_biallelic(record)

            # Find the record that matches the weight
            while record is None or record.pos != weight.position or record.ref != weight.ref or weight.alt != record.alts[0]:
                # If fetch yields no more records at this position, break here
                if record is None:
                    break
                record = next(site_records, None)
                self._check_biallelic(record)

            # Skip the weight if the record was not found
            if record is None:
                self.vcf_num_sites_not_found += 1
                print('  Site not found in VCF:', weight.locus, weight.ref, weight.alt)
                continue
            
            self._vcf_score_dosage(record, weight)
        self._print_wgs_vcf_metrics()
        if self.gvcf_sites_scored is not None:
            self._print_wes_and_wgs_metrics()

    def write_output(self, basename):
        """
        Write the scores to plink-like output files.
        """
        with open(f'{basename}.exome_gvcf.score', 'w') as out_exome_gvcf_score:
            out_exome_gvcf_score.write('#IID\tSCORE1_SUM\n')
            for sample_name in self.sample_names:
                out_exome_gvcf_score.write(f'{sample_name}\t{self.gvcf_sample_score[sample_name]}\n')
        with open(f'{basename}.imputed_wgs_vcf.score', 'w') as out_imputed_wgs_vcf_score:
            out_imputed_wgs_vcf_score.write('#IID\tSCORE1_SUM\n')
            for sample_name in self.sample_names:
                out_imputed_wgs_vcf_score.write(f'{sample_name}\t{self.vcf_sample_score[sample_name]}\n')
        with open(f'{basename}.score', 'w') as out_combined_score:
            out_combined_score.write('#IID	SCORE1_SUM\n')
            for sample_name in self.sample_names:
                out_combined_score.write(f'{sample_name}\t{self.vcf_sample_score[sample_name] + self.gvcf_sample_score[sample_name]}\n')
