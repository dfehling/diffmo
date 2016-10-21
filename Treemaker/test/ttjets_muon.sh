#!/bin/bash

num_jobs=10

for i in `seq $num_jobs`;
do

# python run_top_xs_TreeMaker.py --dirs='/uscms_data/d3/dfehling/NTUPLES/ttjets7/new_truth/'  --config='muon' --outfile=output/ttjets_powheg7_muon_sec_${i}  --sec=${i} --isMC=1 --includeTrigger=1 --totalSec=$num_jobs &
# python run_top_xs_TreeMaker.py --dirs='/uscms_data/d3/dfehling/NTUPLES/ttjets10/new_truth/' --config='muon' --outfile=output/ttjets_powheg10_muon_sec_${i} --sec=${i} --isMC=1 --includeTrigger=1 --totalSec=$num_jobs &
python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/full_powheg/nom/' --outfile=output/ttjets_powhegNC_muon_sec_${i} --config='muon' --sec=${i} --isMC=1 --includeTrigger=1 --totalSec=$num_jobs --includePileup=1 --pileupFile="muon_pu_100" &
# python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/single_top/s/'       --config='muon' --outfile=output/single_top_t_s_muon_sec_${i}     --sec=${i} --isMC=1 --includeTrigger=1 --totalSec=$num_jobs &
# python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/single_top/t/'       --config='muon' --outfile=output/single_top_t_t_muon_sec_${i}     --sec=${i} --isMC=1 --includeTrigger=1 --totalSec=$num_jobs &
# python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/single_top/tW/'      --config='muon' --outfile=output/single_top_t_tW_muon_sec_${i}    --sec=${i} --isMC=1 --includeTrigger=1 --totalSec=$num_jobs &
# python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/single_top/tbar_s/'  --config='muon' --outfile=output/single_top_tbar_s_muon_sec_${i}  --sec=${i} --isMC=1 --includeTrigger=1 --totalSec=$num_jobs &
# python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/single_top/tbar_t/'  --config='muon' --outfile=output/single_top_tbar_t_muon_sec_${i}  --sec=${i} --isMC=1 --includeTrigger=1 --totalSec=$num_jobs &
# python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/single_top/tbar_tW/' --config='muon' --outfile=output/single_top_tbar_tW_muon_sec_${i} --sec=${i} --isMC=1 --includeTrigger=1 --totalSec=$num_jobs &

done
exit
