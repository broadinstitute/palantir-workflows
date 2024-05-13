version 1.0

import "ScoringTasks.wdl"

workflow ScoreVcfWithPreferredGvcf {
    input {
        File exome_gvcf
        File exome_gvcf_index
        File imputed_wgs_vcf
        File imputed_wgs_vcf_index
        String basename
        File weights

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int preemptible = 1
    }

    call ScoreGvcfAndVcf {
        input:
            exome_gvcf = exome_gvcf,
            exome_gvcf_index = exome_gvcf_index,
            imputed_wgs_vcf = imputed_wgs_vcf,
            imputed_wgs_vcf_index = imputed_wgs_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            basename = basename,
            weights = weights,
            preemptible = preemptible
    }

    

    # call CombinePrimaryAndSecondaryScores {
    #     input:
    #         primary_score = PrimaryScore.score,
    #         secondary_score = SecondaryScore.score,
    #         primary_sites_scored = PrimaryScore.sites_scored,
    #         secondary_sites_scored = SecondaryScore.sites_scored,
    #         basename = basename,
    #         preemptible = preemptible
    # }

    output {
        # File score = CombinePrimaryAndSecondaryScores.score
        # File sites_scored = CombinePrimaryAndSecondaryScores.sites_scored
        File exome_gvcf_score = ScoreGvcfAndVcf.exome_gvcf_score
        File imputed_wgs_vcf_score = ScoreGvcfAndVcf.imputed_wgs_vcf_score
    }
}

task ScoreGvcfAndVcf {
    input {
        File exome_gvcf
        File exome_gvcf_index
        File imputed_wgs_vcf
        File imputed_wgs_vcf_index
        String basename
        File weights

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int preemptible = 1
        Int? disk_gb
        Int mem_gb = 4
        Int cpu = 4
    }

    Int disk_size_gb = select_first([disk_gb, ceil(size(exome_gvcf, "GiB") + size(imputed_wgs_vcf, "GiB") + size(weights, "GiB") + size(ref_fasta, "GiB") + 50)])

    command <<<
        set -xeuo pipefail
        
        cat <<'EOF' > script.py
class BGEScorer():
    def __init__(self, ref_dict_path, prs_weights_path):
        self.ref_dict = self._read_ref_dict(ref_dict_path)
        self.prs_weights = self._read_prs_weights(prs_weights_path)
        self.sample_name = None

        self.gvcf_sample_score = None
        self.gvcf_sites_scored = None
        self.gvcf_sites_scored_set = None
        self.gvcf_num_variants_scored = None
        self.gvcf_num_low_quality_variants = None
        self.gvcf_num_ref_blocks_scored = None
        self.gvcf_num_low_quality_ref_blocks = None

        self.vcf_sample_score = None
        self.vcf_sites_scored = None
        self.vcf_num_sites_scored = None
        self.vcf_num_sites_not_found = None

    def _read_ref_dict(self, ref_dict_path):
        with open(ref_dict_path) as ref_dict:
            return [line.split('\t')[1].replace('SN:', '') for line in ref_dict]

    def _read_prs_weights(self, prs_weights_path):
        prs_weights = pd.read_table(prs_weights_path, dtype={'effect_allele': str, 'weight': np.float64, 'chr': str, 'bp': np.int32, 'ref': str, 'alt': str})

        prs_weights[['contig', 'position', 'ref', 'alt']] = prs_weights['CHR:BP:REF:ALT'].str.split(':', expand=True)
        prs_weights = prs_weights.drop(['CHR:BP:REF:ALT'], axis=1)
        prs_weights['locus'] = prs_weights['contig'] + ':' + prs_weights['position']

        prs_weights = prs_weights.astype({'effect_allele': str, 'weight': np.float64, 'contig': str, 'position': np.int32, 'ref': str, 'alt': str})

        return prs_weights.sort_values(by='position').sort_values(kind='stable', by='contig', key=lambda contig: contig.map(self.ref_dict.index))

    def _compare_contigs_and_positions(self, lhs_contig, lhs_position, rhs_contig, rhs_position):
        if lhs_contig != rhs_contig:
            return self.ref_dict.index(lhs_contig) - self.ref_dict.index(rhs_contig)
        return lhs_position - rhs_position

    def _check_biallelic(self, record):
        if record is not None and len(record.alts) > 1:
            raise RuntimeError('Only biallelic variants are supported: ' + str(record))

    def _gvcf_check_quality_and_score_variant(self, record, weight, variant_qual_threshold):
        if record.qual < variant_qual_threshold:
            self.gvcf_num_low_quality_variants += 1
            return 0
        
        if weight.effect_allele not in record.alleles:
            raise RuntimeError('Effect allele not found in variant alleles: ' + str(record))
        
        effect_allele_index = record.alleles.index(weight.effect_allele)
        site_score = record.samples[self.sample_name]['GT'].count(effect_allele_index) * weight.weight

        self.gvcf_num_variants_scored += 1
        self.gvcf_sites_scored.append((weight.locus, weight.ref, weight.alt))
        self.gvcf_sample_score += site_score

    def _gvcf_score_ref_block(self, record, weight, ref_block_gq_threshold):
        if record.samples[self.sample_name]['GQ'] < ref_block_gq_threshold:
            self.gvcf_num_low_quality_ref_blocks += 1
            return 0
        
        site_score = 2 * weight.weight if weight.effect_allele == weight.ref else 0

        self.gvcf_num_ref_blocks_scored += 1
        self.gvcf_sites_scored.append((weight.locus, weight.ref, weight.alt))
        self.gvcf_sample_score += site_score
    
    def _vcf_score_dosage(self, record, weight):
        dosage = record.samples[self.sample_name]['DS'] if record.alts[0] == weight.effect_allele else 2 - record.samples[self.sample_name]['DS']
        site_score = dosage * weight.weight

        self.vcf_num_sites_scored += 1
        self.vcf_sites_scored.append((weight.locus, weight.ref, weight.alt))
        self.vcf_sample_score += site_score

    def score_wes_gvcf(self, gvcf_path, sample_name, variant_qual_threshold=30, ref_block_gq_threshold=30):
        gvcf = pysam.VariantFile(gvcf_path, "r")

        if sample_name not in gvcf.header.samples:
            raise RuntimeError('The sample_name is not present in the GVCF file.')
        self.sample_name = sample_name

        self.gvcf_sample_score = 0
        self.gvcf_sites_scored = []
        self.gvcf_num_variants_scored = 0
        self.gvcf_num_low_quality_variants = 0
        self.gvcf_num_ref_blocks_scored = 0
        self.gvcf_num_low_quality_ref_blocks = 0

        weights_iter = self.prs_weights.iterrows()
        _, current_weight = next(weights_iter)

        for record in gvcf.fetch():
            # Skip all weights that are before the current record
            while self._compare_contigs_and_positions(record.contig, record.pos, current_weight.contig, current_weight.position) > 0:
                try:
                    _, current_weight = next(weights_iter)
                except StopIteration:
                    break
            
            # If the record is a ref block
            if record.alleles_variant_types == ('REF', 'REF'):
                # If the ref block overlaps with the current weight
                if self._compare_contigs_and_positions(record.contig, record.pos, current_weight.contig, current_weight.position) <= 0 and self._compare_contigs_and_positions(record.contig, record.stop, current_weight.contig, current_weight.position) >= 0:
                    self._gvcf_score_ref_block(record, current_weight, ref_block_gq_threshold)
            else: # If the record is a variant
                if self._compare_contigs_and_positions(record.contig, record.pos, current_weight.contig, current_weight.position) == 0:
                    if current_weight.alt in record.alts:
                        self._gvcf_check_quality_and_score_variant(record, current_weight, variant_qual_threshold)
        
        # Save the scored sites as a set for faster lookup
        self.gvcf_sites_scored_set = set(self.gvcf_sites_scored)
        print(f'WES GVCF: {self.gvcf_num_variants_scored} variants scored, {self.gvcf_num_low_quality_variants} low quality variants, {self.gvcf_num_ref_blocks_scored} ref blocks scored, {self.gvcf_num_low_quality_ref_blocks} low quality ref blocks.')

    def score_wgs_vcf(self, wgs_vcf_path, sample_name, allow_wgs_vcf_only=False):
        if not allow_wgs_vcf_only and self.gvcf_sites_scored_set is None:
            raise RuntimeError('WES GVCF scoring must be performed before WGS VCF scoring. If you want to score the WGS VCF only, set allow_wgs_vcf_only to True')
        if self.gvcf_sites_scored_set is None:
            self.gvcf_sites_scored_set = set()

        print(f'{len(self.gvcf_sites_scored_set)} weights have already been scored by the WES GVCF. Scoring the remaining weights using the WGS VCF.')

        vcf = pysam.VariantFile(wgs_vcf_path, "r")
        if self.sample_name is not None and sample_name != self.sample_name:
            raise RuntimeError('Sample name in WES GVCF and WGS VCF do not match.')
        if sample_name not in vcf.header.samples:
            raise RuntimeError('The sample_name is not present in the VCF file.')
        # Set the sample_name in case it was not set by the WES GVCF scoring
        self.sample_name = sample_name

        self.vcf_sample_score = 0
        self.vcf_sites_scored = []
        self.vcf_num_sites_scored = 0
        self.vcf_num_sites_not_found = 0

        for _, weight in self.prs_weights.iterrows():
            if (weight.locus, weight.ref, weight.alt) in self.gvcf_sites_scored_set:
                continue

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
                print('Not found:', weight.locus, weight.ref, weight.alt)
                continue
            
            self._vcf_score_dosage(record, weight)
        print(f'WGS VCF:  {self.vcf_num_sites_scored} sites scored, {self.vcf_num_sites_not_found} sites not found.')

bge_scorer = BGEScorer('~{ref_dict}', '~{weights}')
bge_scorer.score_wes_gvcf('~{exome_gvcf}', '~{basename}')
bge_scorer.score_wgs_vcf('~{imputed_wgs_vcf}', '~{basename}')

with open('~{basename}.exome_gvcf.score', 'w') as out_exome_gvcf_score:
    out_exome_gvcf_score.write(f'{bge_scorer.gvcf_sample_score}\n')
with open('~{basename}.imputed_wgs_vcf.score', 'w') as out_imputed_wgs_vcf_score:
    out_imputed_wgs_vcf_score.write(f'{bge_scorer.vcf_sample_score}\n')
EOF
            python3 script.py
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/pysam:v1"
        memory: mem_gb + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible
    }

    output {
        File exome_gvcf_score = "~{basename}.exome_gvcf.score"
        File imputed_wgs_vcf_score = "~{basename}.imputed_wgs_vcf.score"
    }
}

task FilterSitesAndGetExtractedSiteIDs {
    input {
        File extracted_vcf
        String basename
        Int mem_gb = 4
        Int cpu = 4
        Int? disk_gb
        Int preemptible = 1
    }

    Int disk_size_gb = select_first([disk_gb, ceil(size(extracted_vcf, "GiB") + 50)])

    command <<<
        set -xeuo pipefail

        bcftools view -O z -o ~{basename}.extracted.filtered.vcf.gz -e 'QUAL<=20' ~{extracted_vcf}
        bcftools index ~{basename}.extracted.filtered.vcf.gz

        bcftools query -f '%CHROM:%POS:%REF:%ALT\n' ~{basename}.extracted.filtered.vcf.gz > ~{basename}.extracted.filtered.sites
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.3"
        memory: mem_gb + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible
    }

    output {
        File extracted_filtered_vcf = "~{basename}.extracted.filtered.vcf.gz"
        File extracted_filtered_vcf_index = "~{basename}.extracted.filtered.vcf.gz"
        File extracted_filtered_sites = "~{basename}.extracted.filtered.sites"
    }
}

task CombinePrimaryAndSecondaryScores {
    input {
        File primary_score
        File secondary_score
        File primary_sites_scored
        File secondary_sites_scored

        String basename

        Int disk_gb = 50
        Int mem_gb = 2
        Int cpu = 1
        Int preemptible = 1
    }

    command <<<
        set -xeuo pipefail

        cat ~{primary_sites_scored} ~{secondary_sites_scored} | sort | uniq > ~{basename}.sscore.vars

        cat <<EOF > script.py
import pandas as pd

primary_scores = pd.read_table('~{primary_score}')
secondary_scores = pd.read_table('~{secondary_score}')

scores = primary_scores.merge(secondary_scores, how='inner', on='#IID')
scores['SCORE1_SUM'] = scores['SCORE1_SUM_x'] + scores['SCORE1_SUM_y']
scores['NAMED_ALLELE_DOSAGE_SUM'] = scores['NAMED_ALLELE_DOSAGE_SUM_x'] + scores['NAMED_ALLELE_DOSAGE_SUM_y']
scores['SCORE1_AVG'] = 'NA'
scores[['#IID', 'NAMED_ALLELE_DOSAGE_SUM', 'SCORE1_AVG', 'SCORE1_SUM']].to_csv('~{basename}.sscore', sep='\t', index=None)
EOF
        python3 script.py
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        disks: "local-disk " + disk_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File score = "~{basename}.sscore"
        File sites_scored = "~{basename}.sscore.vars"
    }
}