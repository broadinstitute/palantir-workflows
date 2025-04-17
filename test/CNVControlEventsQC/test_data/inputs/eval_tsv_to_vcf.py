header = '''##fileformat=VCFv4.2
##contig=<ID=chr1,length=?>
##contig=<ID=chr2,length=?>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
'''

with open('eval.vcf', 'w') as vcf:
    vcf.write(header)
    with open('eval.tsv') as tsv:
        next(tsv)
        for line in tsv:
            contig, start, end, svtype, test_type, overlap = line.strip().split('\t')
            vcf.write(f"{contig}\t{start}\t{test_type.replace(' ', '_')}\tN\t<{svtype}>\t100\tPASS\tEND={end};SVTYPE={svtype};OVERLAP_WITH_TRUTH={overlap}\tGT:CN\t0/1:3\n")
