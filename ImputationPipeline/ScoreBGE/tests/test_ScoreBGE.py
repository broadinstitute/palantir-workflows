from unittest import TestCase
from collections import namedtuple

from ImputationPipeline.ScoreBGE import ScoreBGE

TestRecord = namedtuple('TestRecord', ['contig', 'pos', 'stop'])
TestWeight = namedtuple('TestWeight', ['contig', 'position'])


class TestScoreBGE(TestCase):
    def test_read_dict(self):
        s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/test_weights_read.txt')
        self.assertEqual(s.ref_dict, ['1', '3', '2'])

    def test_read_weights(self):
        s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/test_weights_read.txt')
        self.assertEqual(list(s.prs_weights.contig), ['1', '1', '3'])  # contig
        self.assertEqual(list(s.prs_weights.position), [1, 2, 1])  # position
        self.assertEqual(list(s.prs_weights.locus), ['1:1', '1:2', '3:1'])  # locus
        self.assertEqual(list(s.prs_weights.ref), ['A', 'A', 'A'])  # ref
        self.assertEqual(list(s.prs_weights.alt), ['C', 'G', 'T'])  # alt
        self.assertEqual(list(s.prs_weights.effect_allele), ['C', 'A', 'T'])  # effect_allele
        self.assertEqual(list(s.prs_weights.weight), [1.0, 2.0, 4.0])  # weight

    def test_compare_record_and_weight(self):
        s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/test_weights_ref_blocks.txt')

        self.assertEqual(s._compare_record_and_weight(TestRecord('1', 1, 1), TestWeight('1', 1)), 0)
        self.assertLess(s._compare_record_and_weight(TestRecord('1', 1, 1), TestWeight('1', 2)), 0)
        self.assertGreater(s._compare_record_and_weight(TestRecord('1', 2, 2), TestWeight('1', 1)), 0)
        self.assertLess(s._compare_record_and_weight(TestRecord('1', 1, 1), TestWeight('2', 1)), 0)
        self.assertGreater(s._compare_record_and_weight(TestRecord('2', 1, 1), TestWeight('1', 1)), 0)
        self.assertGreater(s._compare_record_and_weight(TestRecord('2', 1, 1), TestWeight('3', 1)), 0)

    def test_gvcf_score_ref_block(self):
        s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/test_weights_ref_blocks.txt')

        s.score_wes_gvcf('resources/test_ref_blocks.gvcf.gz', ['testsample'], 30)
        expected_sites_scored = [
            #                               ('1:100', 'A', 'C'), no call low gq
            ('1:200', 'A', 'C'),  # score 0
            ('1:300', 'A', 'C'),  # score 4
            ('1:310', 'A', 'AC'),  # score 0
            ('1:320', 'AC', 'A'),  # score 16
            #                               ('1:400', 'A', 'C'), negative GQ
            ('1:500', 'A', 'C'),  # score 64
            ('1:700', 'G', 'C'),  # score 0
            #                               ('1:800', 'A', 'C'), gap
            ('1:900', 'A', 'C'),  # score 1024
        ]
        total_expected_score = (0 + 4 + 0 + 16 + 64 + 256 + 1024) * 2

        self.assertEqual(s.gvcf_sites_scored['testsample'], expected_sites_scored)
        self.assertEqual(s.gvcf_sample_score['testsample'], total_expected_score)

    def test_gvcf_score_variants(self):
        s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/test_weights_gvcf_variants.txt')

        s.score_wes_gvcf('resources/test_gvcf_variants.gvcf.gz', ['testsample'], 30)
        expected_sites_scored = [
            ('3:100', 'A', 'C'),  # score 2*1
            ('3:100', 'A', 'C'),  # score 0
            ('3:110', 'A', 'C'),  # score 4
            ('3:110', 'A', 'C'),  # score 8
            ('3:120', 'A', 'C'),  # score 0
            ('3:120', 'A', 'C'),  # score 2*32
            # ('3:130', 'A', 'C'), # score 0
            # ('1:140', 'A', 'C'), # score 0
            # ('3:150', 'A', 'C'), # score 0
            ('3:200', 'A', 'C'),  # score 512
            ('3:300', 'A', 'AC'),  # score 1024
            ('3:400', 'A', 'AC'),  # score 2048
            ('3:500', 'A', 'G'),  # score 4096
            ('3:600', 'A', 'C'),  # score 8192
        ]
        total_expected_score = (2 * 1 + 4 + 8 + 2 * 32 + 512 + 1024 + 2048 + 4096 + 8192)

        self.assertEqual(s.gvcf_sites_scored['testsample'], expected_sites_scored)
        self.assertEqual(s.gvcf_sample_score['testsample'], total_expected_score)
