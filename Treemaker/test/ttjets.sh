#!/bin/bash

#Weight for MC 
#Powheg <700    : 1.00/21675970*245.8*19700
#Powheg 700-1000: 0.074/3082812*245.8*19700    : 0.074 * 1.571
#Powheg 1000+   : 0.014/1249111*245.8*19700    : 0.014 * 3.877

weight=0.22339
weight7=0.11623
weight10=0.05427

num_jobs=20

for i in `seq $num_jobs`;
do

python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/full_powheg/nom/' --outfile=output/ttjets_powhegLT7_sec_${i} --sec=${i} --totalSec=$num_jobs --isMC=1 --unfoldWeight=$weight --invMassCut=700 &
python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/full_powheg/nom/' --outfile=output/ttjets_powhegNC_sec_${i}  --sec=${i} --totalSec=$num_jobs --isMC=1 --unfoldWeight=$weight &
python run_top_xs_TreeMaker.py --dirs='/uscms_data/d3/dfehling/NTUPLES/ttjets7/new_truth/'            --outfile=output/ttjets_powheg7_sec_${i}   --sec=${i} --totalSec=$num_jobs --isMC=1 --unfoldWeight=$weight7 &
python run_top_xs_TreeMaker.py --dirs='/uscms_data/d3/dfehling/NTUPLES/ttjets10/new_truth/'           --outfile=output/ttjets_powheg10_sec_${i}  --sec=${i} --totalSec=$num_jobs --isMC=1 --unfoldWeight=$weight10 &

done
exit

