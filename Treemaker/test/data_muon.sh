#!/bin/bash

num_jobs=20

for i in `seq $num_jobs`;
do 

python run_top_xs_TreeMaker_muon.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/SingleMu/Data_Run_A_v1/150819_221942/0000/' --outfile=output/data2_A_muon_sec_${i} --sec=${i} --totalSec=$num_jobs &
python run_top_xs_TreeMaker_muon.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/SingleMu/Data_Run_B_v1/150819_222250/0000/' --outfile=output/data2_B_muon_sec_${i} --sec=${i} --totalSec=$num_jobs &
python run_top_xs_TreeMaker_muon.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/SingleMu/Data_Run_C_v1/150819_222437/0000/' --outfile=output/data2_C_muon_sec_${i} --sec=${i} --totalSec=$num_jobs &
python run_top_xs_TreeMaker_muon.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/SingleMu/Data_Run_D_v1/150819_223011/0000/' --outfile=output/data2_D_muon_sec_${i} --sec=${i} --totalSec=$num_jobs &

done
