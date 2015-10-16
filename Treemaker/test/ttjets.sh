#!/bin/bash

num_jobs=10

for i in `seq $num_jobs`;
do

python run_top_xs_TreeMaker.py --dirs='/uscms_data/d3/dfehling/NTUPLES/June/ttjets7/'  --outfile=output/ttjets_powheg7_sec_${i}  --sec=${i} --totalSec=$num_jobs &
python run_top_xs_TreeMaker.py --dirs='/uscms_data/d3/dfehling/NTUPLES/June/ttjets10/' --outfile=output/ttjets_powheg10_sec_${i} --sec=${i} --totalSec=$num_jobs &

done
exit
