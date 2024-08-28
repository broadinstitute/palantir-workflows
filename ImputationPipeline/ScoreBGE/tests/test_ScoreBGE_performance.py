from unittest import TestCase
from collections import namedtuple
import tempfile
import os
import subprocess

from ImputationPipeline.ScoreBGE import ScoreBGE

TestRecord = namedtuple('TestRecord', ['contig', 'pos', 'stop'])
TestWeight = namedtuple('TestWeight', ['contig', 'position'])


class TestScoreBGEPerformance(TestCase):
    def test_performance(self):
        #os.environ['GCS_OAUTH_TOKEN'] = subprocess.check_output('gcloud auth application-default print-access-token',
        #                                                        shell=True).decode()
        s = ScoreBGE.BGEScorer('resources/performance/Homo_sapiens_assembly38.dict',
                               "resources/performance/t2d_PRS_weights.hg38.txt.tsv", score_haploid_as_diploid=True, use_emerge_weight_format=True)

        s.score_wes_gvcf('resources/performance/201741830-001_1216091246_blood_blended.wes.hard-filtered.gvcf.gz', site_gq_threshold=30)
        s.score_wgs_vcf("resources/performance/imputation60.imputed.vcf.gz", sample_names=['201741830-001_1216091246_blood_blended'])
        for weight in s.prs_weights
        s.write_output('test_output_write', allow_single_source_scoring=True)
