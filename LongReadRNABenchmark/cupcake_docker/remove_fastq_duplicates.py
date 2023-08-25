import sys

if len(sys.argv) != 2:
	print("Usage: " + sys.argv[0] + " <in.fastq>")
	sys.exit(0)

entries = {}

name = ""
seq = ""
qual = ""

line_n = 0
with open(sys.argv[1]) as infile:
	for line in infile:
		if line_n == 0:
			name = line
		elif line_n == 1:
			seq = line
		elif line_n == 3:
			qual = line
		
		line_n += 1
		if line_n >= 4:
			if name not in entries:
				entries[name] = (seq, qual)
			line_n = 0

outfile = open("out.fastq", 'w')

for key, value in entries.items():
	outfile.write(key)
	outfile.write(value[0])
	outfile.write("+\n")
	outfile.write(value[1])

outfile.close()
