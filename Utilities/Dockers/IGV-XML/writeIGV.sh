#!/bin/bash

ref=$1
echo "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>"
echo "<Session genome=\""$ref"\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"All\" nextAutoscaleGroup=\"4\" path=\"\" version=\"8\">"
echo -e "\t<Resources>"

NAMES=""
FILES=""

while [ "$2" != "" ]; do
	case "$2" in 
		-n ) shift
			 NAMES="$NAMES $2"
			 ;;
		* )  FILES="$FILES $2";;
	esac
	shift
done

names=($NAMES)
files=($FILES)

#if names is empty
if [ -z "$names" ]; then
    for ((i=0;i<${#files[@]};++i)); do
    	echo -e "\t\t<Resource path=\""${files[i]}"\"/>" 
	done    
else
    for ((i=0;i<${#files[@]};++i)); do
    	echo -e "\t\t<Resource path=\""${files[i]}"\" name=\""${names[i]}"\"/>" 
	done  
fi

echo -e "\t</Resources>"
echo "</Session>"

