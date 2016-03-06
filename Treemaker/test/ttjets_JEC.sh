#!/bin/bash

num_jobs=10

for i in `seq $num_jobs`;
do
for s in "JERup" "JERdown" "JESup" "JESdown"
do

python run_top_xs_TreeMaker.py --dirs=/uscms_data/d3/dfehling/NTUPLES/ttjets7/${s}/  --outfile=output/ttjets_powheg7_${s}_sec_${i}  --sec=${i} --totalSec=$num_jobs &
python run_top_xs_TreeMaker.py --dirs=/uscms_data/d3/dfehling/NTUPLES/ttjets10/${s}/ --outfile=output/ttjets_powheg10_${s}_sec_${i} --sec=${i} --totalSec=$num_jobs &

done
done
exit
