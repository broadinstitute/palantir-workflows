#NOT DONE
import sys
# [python] confirm_haplotype_flip_dill.py [unflipped file name] [flipped file name]
print("Arg list: (python)", str(sys.argv))

if len(sys.argv) != 3:
    print("please check command; should have two input vcf files")
    exit()
elif len(sys.argv[1]) < 5 or sys.argv[1][ int(len(sys.argv[1]))-4 : int(len(sys.argv[1])) ] != ".vcf":
    print("check unflipped input file name?")
    exit()
elif len(sys.argv) == 3 and (len(sys.argv[2]) < 5 or sys.argv[2][ int(len(sys.argv[2]))-4 : int(len(sys.argv[2])) ] != ".vcf"):
    print("check flipped input file name?")
    exit()
else:
    print("input command looks fine!")
print("-----------------------------------------------------------")
print("...data being read from files...")
with open (str(sys.argv[1]), "r") as file:
    data1 = file.readlines()

with open (str(sys.argv[2]), "r") as file:
    data2 = file.readlines()
print("-----------------------------------------------------------")
data1c = 0
data2c = 0
print("length of unflipped file:", len(data1))
print("length of flipped file:", len(data2))
verify_order = len(data1)>=len(data2)
if verify_order:
    print("are files in correct order (unflipped file cannot be smaller than flipped file):", verify_order)
else:
    print("ERROR: size of ["+sys.argv[1]+"] < size of ["+sys.argv[2]+"]! please check your files?")
    exit()
print("-----------------------------------------------------------")
print("...comparing the files...")
while data1c < len(data1):
    if data1[data1c][0] == "#":
        if (data1[data1c]) != (data2[data2c]):
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print("a line may be missing?")
            print("on line "+str(data1c)+" in "+sys.argv[1]+":")
            print(data1[data1c])
            print("on line "+str(data2c)+" in "+sys.argv[2]+":")
            print(data2[data2c])
            print("===========================================================")
            exit()
        else:
            data1c = data1c + 1
            data2c = data2c + 1
    else:
        line1_arr = data1[data1c].split('\t')
        line1_segment = line1_arr[9].split(':')
        if line1_arr[0] == "chr22":
            data1c = data1c + 1
            pass
        else:
            line2_arr = data2[data2c].split('\t')
            line2_segment = line2_arr[9].split(':')
            if  (line1_arr[0:1] == line2_arr[0:1]) and (line1_segment[0] == "0/0" or line1_segment[0] == "1/1") and (line1_segment[0] != line2_segment[0]) or \
                (line1_arr[0:1] == line2_arr[0:1]) and ((line1_segment[0] == "0/1" and line2_segment[0] != "1/0") or (line1_segment[0] == "1/0" and line2_segment[0] != "0/1")) or \
                (line2_segment[0] == "?/?") or (line2_segment[0] == "?/?") or \
                (line1_arr[0:1] == line2_arr[0:1]) and (line1_segment[0] == "0|0" or line1_segment[0] == "1|1") and (line1_segment[0] != line2_segment[0]) or \
                (line1_arr[0:1] == line2_arr[0:1]) and ((line1_segment[0] == "0|1" and line2_segment[0] != "1|0") or (line1_segment[0] == "1|0" and line2_segment[0] != "0|1")) or \
                (line2_segment[0] == "?|?") or (line2_segment[0] == "?|?"):
                print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                print("on line "+str(data1c+1)+" in "+sys.argv[1]+":")
                print(data1[data1c])
                print("on line "+str(data2c+1)+" in "+sys.argv[2]+":")
                print(data2[data2c])
                print("===========================================================")
                exit()
            else:
                data1c = data1c + 1
                data2c = data2c + 1

print("-----------------------------------------------------------")
print("haplotypes flipped correctly: yes")
print("===========================================================")