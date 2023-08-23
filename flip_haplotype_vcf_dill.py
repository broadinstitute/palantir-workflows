#notes on running this:
#assumes that input file is a vcf file
#[python] flip_haplotype_vcf_dill.py inputfile.vcf outputfile.vcf
#or
#[python] flip_haplotype_vcf_dill.py inputfile.vcf
#(output file called inputfile_flipped.vcf)

import sys
print("Arg list: (python)", str(sys.argv))

if len(sys.argv) != 3 and len(sys.argv) != 2:
    print("please check command; should have input vcf file and (optional) output vcf file")
    exit()
elif len(sys.argv[1]) < 5 or sys.argv[1][ int(len(sys.argv[1]))-4 : int(len(sys.argv[1])) ] != ".vcf":
    print("check input file name?")
    exit()
elif len(sys.argv) == 3 and (len(sys.argv[2]) < 5 or sys.argv[2][ int(len(sys.argv[2]))-4 : int(len(sys.argv[2])) ] != ".vcf"):
    print("check output file name?")
    exit()
else:
    print 

if len(sys.argv) == 3:
    outputfile = sys.argv[2]
else:
    outputfile = sys.argv[1][0:len(sys.argv)-6] + "_flipped.vcf"


#exit() #include for testing
with open (str(sys.argv[1]), "r") as file:
    data = file.readlines()

with open (outputfile, "w") as f:
    errorchar = []
    for line in data:
        if line[0] == "#":
            f.write(line) #remove for testing
            #pass #include for testing
        elif line[0:5] == "chr22":
            pass
            #skipping this chromosome because it is not within the confidence regions for what i am doing
        else:
            line_arr = line.split('\t')
            linesegment = line_arr[9].split(':')
            new_hap = ""
            for char in linesegment[0]:
                if char == '0':
                    new_hap = new_hap + '1'
                elif char == '1':
                    new_hap = new_hap +  '0'
                elif char == '/' or char == '|':
                    new_hap = new_hap +  char
                else:
                    errorchar.append("["+char+'] in '+line_arr[0]+" "+line_arr[1]+" "+line_arr[9]+line)
            
            #error checking for if the flipping does not work
            if new_hap == linesegment[0]:
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNEXPECTED ERROR?")
                print(line)
                print("original haplotype:",linesegment[0])
                print("flipped haplotype:",new_hap)
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNEXPECTED ERROR?")
                exit()
                
            #copying changed haplotype into the line, to write it to the new file
            linesegment[0] = new_hap
            linesegment = ':'.join(linesegment)
            line_arr[9] = linesegment
            line = '\t'.join(line_arr)
            
            f.write(line) #remove for testing
            #break #include for testing
            #inp = input() #include for testing

#show incorrect-character errors
y = len(errorchar)
print(y, 'unexpected-haplotype-character errors')
if y != 0 :
    x = 0
    while x<(min(5, y)):
        print(errorchar[x])
        x = x + 1
