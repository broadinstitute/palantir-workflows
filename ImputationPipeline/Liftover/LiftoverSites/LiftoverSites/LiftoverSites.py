import hail as hl
import pandas as pd
import subprocess
import collections.abc

def format_input(line, header):
    fields = line.split('\t')
    locus_alleles = fields[0].split(':')
    chrom = locus_alleles[0]
    pos = locus_alleles[1]
    a1 = locus_alleles[2]
    a2 = locus_alleles[3]
    effect_allele = fields[1]
    weight = fields[2]
    site_id = fields[0]
    return InputRow(chrom, pos, a1, a2, effect_allele, {'weight': weight, 'site_id': site_id})

def format_output(sites):
    #CHR:BP:A1:A2    A1      BETA    CHR     BP      A1      A2
    w = sites.key_by()
    w = w.annotate(CHR_BP_REF_ALT=w.locus.contig + ":" + hl.str(w.locus.position) + ":" + w.alleles[0] + ":" + w.alleles[1])
    w = w.select(CHR_BP_REF_ALT=w.CHR_BP_REF_ALT, effect_allele=w.effect, weight=w['info.original_weight'])
    w = w.rename({'CHR_BP_REF_ALT': 'CHR:BP:REF:ALT'})
    return w

class LiftoverSites:
    class LiftoverArguments:
        def __init__(self, run_java_picard_jar_instead_of_gatk:bool, path_to_gakt_or_picard_jar:str, chain_file_path:str, reference_path:str) -> None:
            self.run_java_picard_jar_instead_of_gatk = run_java_picard_jar_instead_of_gatk
            self.path_to_gakt_or_picard_jar = path_to_gakt_or_picard_jar
            self.chain_file_path = chain_file_path
            self.reference_path = reference_path
    
    class InputRow:
        def __init__(self, chrom, pos, a1, a2, effect_allele, custom_columns) -> None:
            if any(v is None for v in [chrom, pos, a1, a2]):
                raise RuntimeError('All of chrom, pos, a1, a2 have to be set.')
            self.chrom = chrom
            self.pos = pos
            self.a1 = a1
            self.a2 = a2
            self.effect_allele = effect_allele
            self.custom_columns = custom_columns
    
    def __init__(self, output_dir:str, input_function:collections.abc.Callable[[str, str], InputRow], output_function:collections.abc.Callable[[hl.Table], hl.Table], annotations_to_save:collections.abc.Iterable[str], liftover_arguments:LiftoverArguments, contains_effect_allele:bool, reference_panel_path:str, input_has_header:bool=True, input_encoding:str ='utf-8', confirm_continue_on_conflicts:bool=True) -> None:
        self.output_dir = output_dir
        self.input_function = input_function
        self.output_function = output_function
        self.annotations_to_save = annotations_to_save
        self.liftover_arguments = liftover_arguments
        self.contains_effect_allele = contains_effect_allele
        self.reference_panel_path = reference_panel_path
        self.input_has_header = input_has_header
        self.input_encoding = input_encoding
        self.reference_panel = None
        self.confirm_continue_on_conflicts = confirm_continue_on_conflicts
    
    @staticmethod
    def _print_rejected_site(df_group, file):
        # Get first site in group, because both entries will be the same site with a different allele order
        site = df_group.iloc[0]
        print(f"Position: {site['#CHROM']}:{site['POS']}, Alleles: {site['REF']}, {site['ALT']}. Reason: {site['FILTER']}", file=file)

    @staticmethod
    def _get_ambiguous_sites(sites):
        var_number_counts = sites.group_by(sites['rsid']).aggregate(n=hl.agg.count())
        ambiguous_var_numbers = var_number_counts.filter(var_number_counts.n > 1)
        ambiguous_sites = sites.key_by('rsid').semi_join(ambiguous_var_numbers).key_by('locus', 'alleles')
        return ambiguous_sites.cache()

    def prepare_liftover(self, input_table, basename):
        with open(input_table, encoding=self.input_encoding, newline='') as input_file:
            with open(f'{self.output_dir}/{basename}.hg19.vcf', 'w') as output_vcf:
                output_vcf.write('##fileformat=VCFv4.2\n')
                output_vcf.write('##INFO=<ID=effect_allele,Number=1,Type=String,Description="Effect allele, either ref or alt.">\n')
                for annotation in self.annotations_to_save:
                    output_vcf.write(f'##INFO=<ID=original_{annotation},Number=1,Type=String>\n')
                output_vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

                input_header_line = next(input_file).strip() if self.input_has_header else None

                for i_line, line in enumerate(input_file):
                    input_row = self.input_function(line.strip(), input_header_line)
                    original_columns = ';'.join([f'original_{annotation}={input_row.custom_columns[annotation]}' for annotation in self.annotations_to_save])
                    if not original_columns and not self.contains_effect_allele:
                        original_columns = '.'
                    
                    chrom = 'X' if input_row.chrom == '23' else input_row.chrom

                    if self.contains_effect_allele:
                        a1_is_effect_allele = None
                        if input_row.effect_allele == input_row.a1:
                            a1_is_effect_allele = True
                        elif input_row.effect_allele == input_row.a2:
                            a1_is_effect_allele = False
                        else:
                            raise RuntimeError(f'Error in input sites: Neither A1 nor A2 are the effect allele: {line}')
                        effect_allele_1 = f"effect_allele={'ref' if a1_is_effect_allele else 'alt'};"
                        effect_allele_2 = f"effect_allele={'alt' if a1_is_effect_allele else 'ref'};"

                    output_vcf.write(f"{chrom}\t{input_row.pos}\t{i_line}\t{input_row.a1}\t{input_row.a2}\t30\tPASS\t{effect_allele_1 if self.contains_effect_allele else ''}{original_columns}\n")
                    output_vcf.write(f"{chrom}\t{input_row.pos}\t{i_line}\t{input_row.a2}\t{input_row.a1}\t30\tPASS\t{effect_allele_2 if self.contains_effect_allele else ''}{original_columns}\n")

    def run_liftover(self, basename):
        java_or_not = ['java'] if self.liftover_arguments.run_java_picard_jar_instead_of_gatk else []
        process = subprocess.run(java_or_not + [self.liftover_arguments.path_to_gakt_or_picard_jar,
                    'LiftoverVcf',
                    '-I', f'{self.output_dir}/{basename}.hg19.vcf',
                    '-O', f'{self.output_dir}/{basename}.hg38.vcf',
                    '-CHAIN', self.liftover_arguments.chain_file_path,
                    '-REJECT', f'{self.output_dir}/{basename}.reject.vcf',
                    '-R', self.liftover_arguments.reference_path,
                    '-RECOVER_SWAPPED_REF_ALT', 'false'])
        if process.returncode != 0:
            raise RuntimeError(f'Picard liftover failed. Error {process.returncode}')

    def check_rejected_sites(self, basename):
        # Count header lines to skip when reading VCF
        skip_rows = 0
        with open(f'{self.output_dir}/{basename}.reject.vcf', 'r') as file:
            for line in file:
                if line.startswith('##'):
                    skip_rows += 1
                else:
                    break

        # Read reject VCF
        rejected_sites = pd.read_table(f'{self.output_dir}/{basename}.reject.vcf', skiprows=skip_rows)
        # Filter out all MismatchedRefAllele sites, those were sites where the REF allele wasn't correct
        # Because we submit each site with both allele orders to the liftover, we expect ~50% to fail for this reason
        rejected_sites = rejected_sites[rejected_sites.FILTER != 'MismatchedRefAllele']

        if len(rejected_sites) > 0:
            # Group sites by var_number (rsid)
            rejected_site_groups = rejected_sites.groupby('ID')
            print(f'{len(rejected_site_groups)} site(s) could not be lifted over.')
            print(f'If you want to fix these sites manually, download the following file, edit it, and re-upload it:\n{dir}/{basename}.hg38.vcf\n')
            print(f'Failed sites written to: {self.output_dir}/{basename}.reject.txt')
            with open(f'{self.output_dir}/{basename}.reject.txt', 'w') as file:
                rejected_site_groups.apply(LiftoverSites._print_rejected_site, file=file)
            if self.confirm_continue_on_conflicts:
                print('\nHit enter when you are ready to continue...')
                input()
                print('Continuing.')
        else:
            print('No rejected sites.')

    def disambiguate(self, basename):
        # Init hail
        hl.init(default_reference='GRCh38', idempotent=True)

        # Read lifted over VCF
        sites = hl.import_vcf(f'{self.output_dir}/{basename}.hg38.vcf').rows().flatten().cache()
        if self.reference_panel is None:
            self.reference_panel = hl.import_vcf(self.reference_panel_path, force_bgz=True).rows().cache()

        # Filter out sites that are not in the reference panel
        sites_in_panel = sites.key_by('locus', 'alleles').semi_join(self.reference_panel).cache()

        # Find ambiguous sites
        ambiguous_sites_in_panel = LiftoverSites._get_ambiguous_sites(sites_in_panel)
        # Annotate sites with AF from reference panel
        ambiguous_sites_in_panel = ambiguous_sites_in_panel.annotate(af=self.reference_panel[ambiguous_sites_in_panel.locus, ambiguous_sites_in_panel.alleles].info.AF[0]).cache()
        num_ambiguous_sites_in_panel = ambiguous_sites_in_panel.count()

        # Select sites based on AF filtering threshold
        ambiguous_sites_in_panel_selected = ambiguous_sites_in_panel.filter(ambiguous_sites_in_panel.af > 5e-2).cache()
        num_ambiguous_sites_in_panel_selected = ambiguous_sites_in_panel_selected.count()

        # This will catch if we couldn't disambiguate any sites
        still_ambiguous_sites = LiftoverSites._get_ambiguous_sites(ambiguous_sites_in_panel_selected)
        if still_ambiguous_sites.count() > 0:
            print('AF filtering criterion did not disambiguate all sites.')
            ambiguous_sites_in_panel_selected.export(f'{self.output_dir}/{basename}.still_ambiguous.tsv')
            print('Exported still ambiguous sites to file.')
            ambiguous_sites_in_panel_selected = ambiguous_sites_in_panel_selected.key_by('rsid').anti_join(still_ambiguous_sites.key_by('rsid')).key_by('locus', 'alleles')
        
        # If we passed the above then this will catch if we might have accidentally filtered out both candidates.
        # It may seem that the check below will catch both cases, but it could happen that we filter out both
        # candidates for one site and not disambiguate another site, which will result in the two terms below
        # being equal, so we need both checks
        if 2 * num_ambiguous_sites_in_panel_selected != num_ambiguous_sites_in_panel:
            print('AF filtering criterion filtered out both candidates for a site.')

        print(f'Successfully disambiguated {num_ambiguous_sites_in_panel_selected} sites.')

        # Remove all ambiguous sites and add only the selected sites
        sites_in_panel_without_ambiguous_sites = sites_in_panel.anti_join(ambiguous_sites_in_panel).cache()
        sites_disambiguated = sites_in_panel_without_ambiguous_sites.union(ambiguous_sites_in_panel_selected.drop('af')).cache()

        # If we use an effect_allele, write the correct allele into the "effect" column
        if self.contains_effect_allele:
            sites_disambiguated = sites_disambiguated.annotate(effect=hl.if_else(sites_disambiguated['info.effect_allele'] == 'ref', sites_disambiguated.alleles[0], sites_disambiguated.alleles[1])).cache()

        # Re-structure table for output
        sites_output = self.output_function(sites_disambiguated)
        sites_output.export(f'{self.output_dir}/{basename}.tsv')

    def run_all(self, input_table, basename):
        self.prepare_liftover(input_table, basename)
        self.run_liftover(basename)
        self.check_rejected_sites(basename)
        self.disambiguate(basename)

if __name__ == '__main__':
    liftover = LiftoverSites(
        'test/output', format_input, format_output, ['weight'],
        LiftoverSites.LiftoverArguments(False, '/Users/mgatzen/code/gatk/gatk', '/Users/mgatzen/liftover/b37ToHg38.over.chain', '/Users/mgatzen/reference/hg38/Homo_sapiens_assembly38.fasta'),
        True, 'gs://fc-f0032a5f-c108-4b69-aeb0-9d665405bfdc/_data/1000G_HGDP_no_singletons.sites.vcf.gz')

    liftover.run_all('test/input/hcl.txt', 'hcl2')
