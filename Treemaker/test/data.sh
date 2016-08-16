#!/bin/bash

num_jobs=20

for i in `seq $num_jobs`;
do 

python run_top_xs_TreeMaker.py --dirs='/uscms_data/d3/dfehling/NTUPLES/dataA/' --outfile=output/data_A_sec_${i} --sec=${i} --totalSec=$num_jobs &
python run_top_xs_TreeMaker.py --dirs='/uscms_data/d3/dfehling/NTUPLES/dataB/' --outfile=output/data_B_sec_${i} --sec=${i} --totalSec=$num_jobs &
python run_top_xs_TreeMaker.py --dirs='/uscms_data/d3/dfehling/NTUPLES/dataC/' --outfile=output/data_C_sec_${i} --sec=${i} --totalSec=$num_jobs &
python run_top_xs_TreeMaker.py --dirs='/uscms_data/d3/dfehling/NTUPLES/dataD/' --outfile=output/data_D_sec_${i} --sec=${i} --totalSec=$num_jobs &

done
exit
