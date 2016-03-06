#!/bin/bash

num_jobs=20

for i in `seq $num_jobs`;
do 

python run_top_xs_TreeMaker.py --dirs='/eos/uscms/store/user/dfehling/Jet/Data_Run_A_v5/150618_161327/0000/'   --outfile=output/data_A_sec_${i} --sec=${i} --totalSec=$num_jobs &
python run_top_xs_TreeMaker.py --dirs='/eos/uscms/store/user/dfehling/JetHT/Data_Run_B_v5/150618_155820/0000/' --outfile=output/data_B_sec_${i} --sec=${i} --totalSec=$num_jobs &
python run_top_xs_TreeMaker.py --dirs='/eos/uscms/store/user/dfehling/JetHT/Data_Run_C_v5/150618_160451/0000/' --outfile=output/data_C_sec_${i} --sec=${i} --totalSec=$num_jobs &
python run_top_xs_TreeMaker.py --dirs='/eos/uscms/store/user/dfehling/JetHT/Data_Run_D_v5/150618_161308/0000/' --outfile=output/data_D_sec_${i} --sec=${i} --totalSec=$num_jobs &

done
exit
