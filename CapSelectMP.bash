#!/bin/bash

########################################################################
##                                                                    ##
## This script fully automates CapSelect ranging from the generation  ##
## of the necessary input files (ICM), the execution of the           ##
## algorithm, and the generation of frags_for_enum.sdf for enumeration##
##                                                                    ##
## Version 02082022  -  Antonina L. Nazarova   nazarova@usc.edu       ## 
##                                                                    ##
########################################################################

# Make CapSelect files directory - All input and output files of CapSelect are stored here
mkdir CapSelect_files

# Execute icm_generate_CapSelect_files.icm using ICM (icm-3.9-2b)
/usr/icm-3.9-2b/icm64 icm_generate_CapSelect_files.icm > output_icm_generate_CapSelect_files.log

# Fallback options if no inputs are passed along with the bash script
nproc=16                # Default is 16 CPU cores
sdfPr1=CapSelect_files/protein1.sdf   # Default names & path of input files 
sdfPr2=CapSelect_files/protein2.sdf
sdfFrag=CapSelect_files/fragments.sdf

while getopts f:P:p:c: flag
do
    case "${flag}" in
        f) sdfFrag=${OPTARG};;
        P) sdfPr1=${OPTARG};;
        p) sdfPr2=${OPTARG};;
        c) nproc=${OPTARG};;
    esac
done

echo "NOTE: Starting CapSelectMP ..."
# Check if input files exist
if [ ! -e "$sdfFrag" ]; then
    echo "Please provide the name/path for the fragment sdf using -f FRAGMENTS.sdf"
    exit
fi

[ ! -e "$sdfPr1" ] && echo "No protein1.sdf file" && exit
[ ! -e "$sdfPr2" ] && echo "No protein2.sdf file, running CapSelect for 3D Model."

# In case fragments.sdf was generated using Windows
sed -i -e 's/\r$//' $sdfFrag

nof_frags=$(grep -c \$\$\$\$ $sdfFrag)
chunksize=$(((nof_frags+nproc-1)/nproc))

echo "Number of Processors: $nproc"
echo "Number of Fragments: $nof_frags"
echo "Chunk Size: $chunksize"

echo "Splitting fragments into $nproc Chunks ..."
# Splits fragment.sdf into chunks 
awk -v var=$chunksize '/^\$\$\$\$$/ { if(++delim % var == 0) { next } } { file = sprintf("CapSelect_files/fragment_%s.sdf", int(delim / var )); print > file; }' < $sdfFrag
echo "Successfully created $nproc chunks."

mkdir -p CapSelect_files/tmp
for ((i=0; i<$nproc; i++));
do
  chmod +x CapSelect_files/fragment_${i}.sdf
  [ $i -lt $(($nproc-1)) ] && printf "\$\$\$\$\n" >> CapSelect_files/fragment_${i}.sdf

  mkdir -p CapSelect_files/tmp/CapSelect_${i}
  mv CapSelect_files/fragment_${i}.sdf CapSelect_files/tmp/CapSelect_${i}/fragments.sdf
  cp $sdfPr1 CapSelect_files/tmp/CapSelect_${i}/protein1.sdf
  [ -e "$sdfPr2" ] && cp $sdfPr2 CapSelect_files/tmp/CapSelect_${i}/protein2.sdf

  # Check for CapSelect.py or compiled CapSelect and execute accordingly
  if [ -e "CapSelect.py" ]; then
    cp CapSelect.py CapSelect_files/tmp/CapSelect_${i}/
    cd CapSelect_files/tmp/CapSelect_${i}
    python CapSelect.py > output_${i}.log &
  elif [ -e "CapSelect" ]; then
    cp CapSelect CapSelect_files/tmp/CapSelect_${i}/
    cd CapSelect_files/tmp/CapSelect_${i}
    ./CapSelect > output_${i}.log &
  else
    echo "CapSelect.py or CapSelect executable not found."
    exit 1
  fi

  let PID${i}=${!}
  cd ../../../
done

echo "Started $nproc CapSelect instances."
echo "Please wait ..."
wait
echo "NOTE: All threads are completed"

find . -name CapSelect.sdf |

