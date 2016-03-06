#!/bin/bash

num_jobs=20

for i in `seq $num_jobs`;
do 

python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/muon/runA/' --config='muon' --outfile=output/data_A_muon_sec_${i} --sec=${i} --totalSec=$num_jobs &
python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/muon/runB/' --config='muon' --outfile=output/data_B_muon_sec_${i} --sec=${i} --totalSec=$num_jobs &
python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/muon/runC/' --config='muon' --outfile=output/data_C_muon_sec_${i} --sec=${i} --totalSec=$num_jobs &
python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/muon/runD/' --config='muon' --outfile=output/data_D_muon_sec_${i} --sec=${i} --totalSec=$num_jobs &

done