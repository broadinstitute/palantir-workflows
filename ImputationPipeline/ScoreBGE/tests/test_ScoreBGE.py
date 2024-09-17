from unittest import TestCase
from collections import namedtuple
import tempfile

from ImputationPipeline.ScoreBGE import ScoreBGE

TestRecord = namedtuple('TestRecord', ['contig', 'pos', 'stop'])
TestWeight = namedtuple('TestWeight', ['contig', 'position'])


class TestScoreBGE(TestCase):
    def test_read_dict(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/test_weights_read.txt', temp_dir + '/test_output', False)
            self.assertEqual(s.ref_dict, ['1', '3', '2'])

    def test_read_weights(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/test_weights_read.txt', temp_dir + '/test_output', False)
            self.assertEqual(list(s.prs_weights.contig), ['1', '1', '3'])  # contig
            self.assertEqual(list(s.prs_weights.position), [1, 2, 1])  # position
            self.assertEqual(list(s.prs_weights.ref), ['A', 'A', 'A'])  # ref
            self.assertEqual(list(s.prs_weights.alt), ['C', 'G', 'T'])  # alt
            self.assertEqual(list(s.prs_weights.effect_allele), ['C', 'A', 'T'])  # effect_allele
            self.assertEqual(list(s.prs_weights.weight), [1.0, 2.0, 4.0])  # weight

    def test_compare_record_and_weight(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/wes_gvcf_ref_blocks/test_weights_ref_blocks.txt', temp_dir + '/test_output', False)

            self.assertEqual(s._compare_record_and_weight(TestRecord('1', 1, 1), TestWeight('1', 1)), 0)
            self.assertLess(s._compare_record_and_weight(TestRecord('1', 1, 1), TestWeight('1', 2)), 0)
            self.assertGreater(s._compare_record_and_weight(TestRecord('1', 2, 2), TestWeight('1', 1)), 0)
            self.assertLess(s._compare_record_and_weight(TestRecord('1', 1, 1), TestWeight('2', 1)), 0)
            self.assertGreater(s._compare_record_and_weight(TestRecord('2', 1, 1), TestWeight('1', 1)), 0)
            self.assertGreater(s._compare_record_and_weight(TestRecord('2', 1, 1), TestWeight('3', 1)), 0)

    def test_gvcf_score_ref_block(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/wes_gvcf_ref_blocks/test_weights_ref_blocks.txt', temp_dir + '/test_output', False)

            s.score_wes_gvcf('resources/wes_gvcf_ref_blocks/test_ref_blocks.gvcf.gz', ['testsample'], 30)
            expected_sites_scored = [
                #                               ('1', 100, 'A', 'C'), no call low gq
                ('1', 200, 'A', 'C'),  # score 0
                ('1', 300, 'A', 'C'),  # score 4
                ('1', 310, 'A', 'AC'),  # score 0
                ('1', 320, 'AC', 'A'),  # score 16
                #                               ('1', 400, 'A', 'C'), negative GQ
                ('1', 500, 'A', 'C'),  # score 64
                ('1', 700, 'G', 'C'),  # score 0
                #                               ('1', 800, 'A', 'C'), gap
                ('1', 900, 'A', 'C'),  # score 1024
            ]
            total_expected_score = (0 + 4 + 0 + 16 + 64 + 256 + 1024) * 2

            self.assertEqual(s.gvcf_sites_scored['testsample'], expected_sites_scored)
            self.assertEqual(s.gvcf_sample_score['testsample'], total_expected_score)

    def test_gvcf_score_variants(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            for score_haploid_as_diploid in [False, True]:
                s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/wes_gvcf_variants/test_weights_gvcf_variants.txt', temp_dir + '/test_output', score_haploid_as_diploid)

                s.score_wes_gvcf('resources/wes_gvcf_variants/test_gvcf_variants.gvcf.gz', ['testsample'], 30)
                expected_sites_scored = [
                    ('3', 100, 'A', 'C'),  # score 2*1
                    ('3', 100, 'A', 'C'),  # score 0
                    ('3', 110, 'A', 'C'),  # score 4
                    ('3', 110, 'A', 'C'),  # score 8
                    ('3', 120, 'A', 'C'),  # score 0
                    ('3', 120, 'A', 'C'),  # score 2*32
                    # ('3', 130, 'A', 'C'), # score 0
                    # ('1', 140, 'A', 'C'), # score 0
                    # ('3', 150, 'A', 'C'), # score 0
                    ('3', 200, 'A', 'C'),  # score 512
                    ('3', 300, 'A', 'AC'),  # score 1024
                    ('3', 400, 'A', 'AC'),  # score 2048
                    ('3', 500, 'A', 'G'),  # score 4096
                    ('3', 600, 'A', 'C'),  # score 8192 if chrX is scored as haploid, otherwise double
                    ('3', 700, 'A', 'G'),  # score 0
                ]
                total_expected_score = (2 * 1 + 4 + 8 + 2 * 32 + 512 + 1024 + 2048 + 4096 + 8192 * (2 if score_haploid_as_diploid else 1))

                self.assertEqual(expected_sites_scored, s.gvcf_sites_scored['testsample'], f'{score_haploid_as_diploid=}')
                self.assertEqual(total_expected_score, s.gvcf_sample_score['testsample'], f'{score_haploid_as_diploid=}')

    def test_score_wgs_vcf(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/wgs_vcf/test_weights_wgs_vcf.txt', temp_dir + '/test_output', False)

            s.score_wgs_vcf('resources/wgs_vcf/test_wgs.vcf.gz', ['testsample'], allow_wgs_vcf_only=True)
            expected_sites_scored = [
                ('3', 100, 'A', 'C'),  # score 2*1
                ('3', 100, 'A', 'C'),  # score 0
                ('3', 110, 'A', 'C'),  # score 1*4
                ('3', 110, 'A', 'C'),  # score 1*8
                ('3', 120, 'A', 'C'),  # score 0
                ('3', 120, 'A', 'C'),  # score 2*32
                ('3', 130, 'A', 'C'),  # score 2*64
                ('3', 140, 'A', 'C'),  # score 0.1*128
                ('3', 150, 'A', 'C'),  # score 1*256
                ('3', 300, 'A', 'AC'),  # score 1*512
                ('3', 400, 'AC', 'A'),  # score 1*1024
                ('3', 600, 'A', 'C'),  # score 1*2048
            ]
            total_expected_score = (2*1 + 1*4 + 1*8 + 2*32 + 2*64 + 0.1*128 + 1*256 + 1*512 + 1*1024 + 1*2048)

            self.assertEqual(s.vcf_sites_scored['testsample'], expected_sites_scored)
            self.assertAlmostEqual(s.vcf_sample_score['testsample'], total_expected_score, 4)

    def test_combined(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/combined/test_weights_combined.txt', temp_dir + '/test_output', False)

            s.score_wes_gvcf('resources/combined/test_combined.wes.gvcf.gz', site_gq_threshold=30)
            s.score_wgs_vcf('resources/combined/test_combined.wgs.vcf.gz')

            s.write_output()

            files_to_compare = [
                'exome_gvcf.score',
                'exome_gvcf.sites_scored',
                'imputed_wgs_vcf.score',
                'imputed_wgs_vcf.sites_scored',
                'any_source_any_sample.sites_scored'
            ]

            for file_to_compare in files_to_compare:
                test_filename = f'{temp_dir}/test_output.{file_to_compare}'
                expected_filename = f'resources/combined/expected_output/expected_output.{file_to_compare}'
                with open(test_filename) as test_output:
                    with open(expected_filename) as expected_output:
                        self.assertListEqual(list(test_output), list(expected_output), file_to_compare)

            # Expected score: 152
            # locus score source
            # 1', 100     0   GVCF
            # 1', 200   2*4    VCF
            # 1', 300  1*16   GVCF
            # 1', 400  2*64    VCF
            # 1', 500     0   none
            with open(f'{temp_dir}/test_output.score') as test_output:
                with open('resources/combined/expected_output/expected_output.score') as expected_output:
                    self.assertListEqual(list(test_output), list(expected_output))

