import sys

num_called_methylated_cpg_sites = 0
num_called_unmethylated_cpg_sites = 0

with open(sys.argv[1]) as infile:
	for line in infile:
		line = line.rstrip()
		columns = line.split()
		if columns[5] == "CG":
			num_called_methylated_cpg_sites += int(columns[3])
			num_called_unmethylated_cpg_sites += int(columns[4])

called_cpg = num_called_methylated_cpg_sites / (num_called_methylated_cpg_sites + num_called_unmethylated_cpg_sites) * 100
print(called_cpg)
