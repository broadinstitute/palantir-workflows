#notes on running this:
#assumes that input file is a vcf file - normalized? yes
#[python] flip_haplotype_vcf_dill.py inputfile.vcf outputfile.vcf
#or
#[python] flip_haplotype_vcf_dill.py inputfile.vcf
#(output file called inputfile_flipped.vcf)

import sys

print("Arg list: (python)", " ".join(sys.argv))

if len(sys.argv) != 3 and len(sys.argv) != 2:
    print("ERROR: please check command; should have input vcf file and (optional) output vcf file")
    exit()
elif len(sys.argv[1]) < 5 or sys.argv[1][ int(len(sys.argv[1]))-4 : int(len(sys.argv[1])) ] != ".vcf":
    print("ERROR: check input file name?")
    exit()
elif len(sys.argv) == 3 and (len(sys.argv[2]) < 5 or sys.argv[2][ int(len(sys.argv[2]))-4 : int(len(sys.argv[2])) ] != ".vcf"):
    print("ERROR: check output file name?")
    exit()
else:
    print("input command structure looks fine!")

if len(sys.argv) == 3:
    outputfile = sys.argv[2]
else:
    outputfile = sys.argv[1][0:len(sys.argv)-6] + "_flipped.vcf"

print("-----------------------------------------------------------")
print("...processing vcf file...")
try:
    with open (str(sys.argv[1]), "r") as file:
        data = file.readlines()
except:
    print("ERROR: file ["+sys.argv[1]+"] does not appear to exist! please check your files?")
    print("===========================================================")
    exit()

with open (outputfile, "w") as f:
    errorchar = []
    errorproc = []
    processed = 0
    passed_contigs = 0
    passed_chr22 = 0
    proc_variant = 0
    for line in data:
        processed = processed+1
        if line[0] == "#":
            passed_contigs = passed_contigs + 1
            f.write(line)
        elif line[0:5] == "chr22":
            print('chr22: RMVD')
            passed_chr22 = passed_chr22 + 1
            pass
            #skipping this chromosome because it is not within the confidence regions for what i am doing
        else:
            line_arr = line.split('\t')
            linesegment = line_arr[9].split(':')
            if len(linesegment[0]) != 3 or (linesegment[0][0] != "0" and linesegment[0][0] != "1") or (linesegment[0][1] != "/" and linesegment[0][1] != "|") or (linesegment[0][2] != "0" and linesegment[0][2] != "1"):
                new_hap = "?/?"
                errorchar.append("["+linesegment[0]+'] replaced with '+new_hap+' in '+line_arr[0]+" "+line_arr[1]+" "+line_arr[9]+line)
            else:
                new_hap = linesegment[0][2]+linesegment[0][1]+linesegment[0][0]
            
            print(linesegment[0],'-->', new_hap)
            if ((linesegment[0] == "0/0" or linesegment[0] == "1/1") and new_hap != linesegment[0]) or ((linesegment[0] == "0/1" and new_hap != "1/0") or (linesegment[0] == "1/0" and new_hap != "0/1")):
                errorproc.append("orig haplotype:["+linesegment[0]+'], new haplotype ['+new_hap+'] at '+line_arr[0]+" "+line_arr[1]+" "+line_arr[9]+line)
                new_hap = "?/?"

            #copying changed haplotype into the line, to write it to the new file
            linesegment[0] = new_hap
            linesegment = ':'.join(linesegment)
            line_arr[9] = linesegment
            line = '\t'.join(line_arr)
            
            f.write(line)
            proc_variant = proc_variant + 1
print("-----------------------------------------------------------")
print("lines in file:", processed)
print("contig lines:", passed_contigs)
print("chr22 lines removed:", passed_chr22)
print("variants outputted:",proc_variant)
print("all variants processed:",processed == passed_contigs+passed_chr22+proc_variant)
print("-----------------------------------------------------------")
#show incorrect-character errors
y = len(errorchar)
print(y, 'haplotype-character error(s) out of',processed,"lines processed")
if y != 0 :
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    x = 0
    while x<(min(6, y)):
        print(errorchar[x])
        x = x + 1
print("-----------------------------------------------------------")
y = len(errorproc)
print(y, 'processing error(s) out of',processed,"lines processed")
if y != 0 :
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    x = 0
    while x<(min(6, y)):
        print(errorchar[x])
        x = x + 1
print("===========================================================")
