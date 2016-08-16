#!/bin/bash

weight=0.22339
weight7=0.11623
weight10=0.05427

num_jobs=20

for i in `seq $num_jobs`;
do
for s in "JERup" "JERdown" "JESup" "JESdown"
do

python run_top_xs_TreeMaker.py --dirs=/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/full_powheg/${s}/  --outfile=output/ttjets_powhegLT7_${s}_sec_${i} --sec=${i} --totalSec=$num_jobs --isMC=1 --unfoldWeight=$weight   --doUnfold=0 --invMassCut=700 &
python run_top_xs_TreeMaker.py --dirs=/uscms_data/d3/dfehling/NTUPLES/ttjets7/${s}/                   --outfile=output/ttjets_powheg7_${s}_sec_${i}   --sec=${i} --totalSec=$num_jobs --isMC=1 --unfoldWeight=$weight7  --doUnfold=0 &
python run_top_xs_TreeMaker.py --dirs=/uscms_data/d3/dfehling/NTUPLES/ttjets10/${s}/                  --outfile=output/ttjets_powheg10_${s}_sec_${i}  --sec=${i} --totalSec=$num_jobs --isMC=1 --unfoldWeight=$weight10 --doUnfold=0 &

done
done
exit
