import sys
sys.path.insert(0, '/Users/mgatzen/code/palantir-workflows')
import ImputationPipeline.ScoreBGE.ScoreBGE as ScoreBGE
from collections import namedtuple

import pysam

TestRecord = namedtuple('TestRecord', ['contig', 'pos', 'stop'])
TestWeight = namedtuple('TestWeight', ['contig', 'position'])

def test_read_dict():
    s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/test_weights_read.txt')
    assert s.ref_dict == ['1', '3', '2']

def test_read_weights():
    s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/test_weights_read.txt')
    assert list(s.prs_weights.contig) == ['1', '1', '3'] # contig
    assert list(s.prs_weights.position) == [1, 2, 1] # position
    assert list(s.prs_weights.locus) == ['1:1', '1:2', '3:1'] # locus
    assert list(s.prs_weights.ref) == ['A', 'A', 'A'] # ref
    assert list(s.prs_weights.alt) == ['C', 'G', 'T'] # alt
    assert list(s.prs_weights.effect_allele) == ['C', 'A', 'T'] # effect_allele
    assert list(s.prs_weights.weight) == [1.0, 2.0, 4.0] # weight

def test_compare_record_and_weight():
    s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/test_weights_ref_blocks.txt')


    assert s._compare_record_and_weight(TestRecord('1', 1, 1), TestWeight('1', 1)) == 0
    assert s._compare_record_and_weight(TestRecord('1', 1, 1), TestWeight('1', 2)) < 0
    assert s._compare_record_and_weight(TestRecord('1', 2, 2), TestWeight('1', 1)) > 0
    assert s._compare_record_and_weight(TestRecord('1', 1, 1), TestWeight('2', 1)) < 0
    assert s._compare_record_and_weight(TestRecord('2', 1, 1), TestWeight('1', 1)) > 0
    assert s._compare_record_and_weight(TestRecord('2', 1, 1), TestWeight('3', 1)) > 0

def test_gvcf_score_ref_block():
    s = ScoreBGE.BGEScorer('resources/ref.dict', 'resources/test_weights_ref_blocks.txt')

    s.score_wes_gvcf('resources/test_ref_blocks.gvcf.gz', ['testsample'], 30)
    expected_sites_scored = [
        #                               ('1:100', 'A', 'C'), no call low gq
        ('1:200', 'A', 'C'), # score 0
        ('1:300', 'A', 'C'), # score 4
        ('1:310', 'A', 'AC'), #score 0
        ('1:320', 'AC', 'A'), #score 16
        #                               ('1:400', 'A', 'C'), negative GQ
        ('1:500', 'A', 'C'), # score 64
        ('1:700', 'G', 'C'), # score 0
        #                               ('1:800', 'A', 'C'), gap
        ('1:900', 'A', 'C'), # score 1024
    ]
    total_expected_score = (0 + 4 + 0 + 16 + 64 + 256 + 1024) * 2

    assert s.gvcf_sites_scored['testsample'] == expected_sites_scored
    assert s.gvcf_sample_score['testsample'] == total_expected_score

    s.write_output('test_output', allow_single_source_scoring=True)


if __name__ == '__main__':
    test_read_dict()
    test_read_weights()
    test_compare_record_and_weight()
    test_gvcf_score_ref_block()
    print('All tests passed!')