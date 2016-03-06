#!/bin/bash

num_jobs=10

for i in `seq $num_jobs`;
do

python run_top_xs_TreeMaker.py --dirs='/eos/uscms/store/user/dfehling/TT_Mtt-700to1000_CT10_TuneZ2star_8TeV-powheg-tauola/TTBar_Powheg7_v4/150618_161425/0000/'  --outfile=output/ttjets_powheg7_sec_${i}  --sec=${i} --totalSec=$num_jobs &
python run_top_xs_TreeMaker.py --dirs='/eos/uscms/store/user/dfehling/TT_Mtt-1000toInf_CT10_TuneZ2star_8TeV-powheg-tauola/TTBar_Powheg10_v4/150618_161443/0000/' --outfile=output/ttjets_powheg10_sec_${i} --sec=${i} --totalSec=$num_jobs &

done
exit
