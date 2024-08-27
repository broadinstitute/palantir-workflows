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
        os.environ['GCS_OAUTH_TOKEN'] = subprocess.check_output('gcloud auth application-default print-access-token',
                                                                shell=True).decode()
        s = ScoreBGE.BGEScorer('resources/performance/Homo_sapiens_assembly38.dict',
                               "resources/performance/t2d_PRS_weights.hg38.txt.tsv", score_haploid_as_diploid=True, )

        s.score_wes_gvcf("gs://fc-f0032a5f-c108-4b69-aeb0-9d665405bfdc/PDO-28184_hg38_hgdp_1kg/option2/201741830-001_1216091246_blood_blended/201741830-001_1216091246_blood_blended.wes.hard-filtered.gvcf.gz", site_gq_threshold=30)
        s.score_wgs_vcf("gs://fc-f0032a5f-c108-4b69-aeb0-9d665405bfdc/submissions/90650b26-3d08-4088-add0-88cf15e505ed/Glimpse2Imputation/85476e81-c1bc-4f55-879a-417ee32da793/call-GlimpseLigate/imputation60.imputed.vcf.gz")
