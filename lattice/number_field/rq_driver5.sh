#!/bin/bash
cd /home/ubuntu/code/math-modular-sources/lattice/number_field
./rq_scan 2 rq_fields3.txt V 2000 200 rq_out_R2_normV.txt -0.35 2 > rq_log_R2_normV.txt 2>&1
echo V2 >> rq_stage5.txt
./rq_scan 3 rq_fields3.txt W 500 150 rq_out_R3_normW.txt -0.35 3 > rq_log_R3_normW.txt 2>&1
echo W3 >> rq_stage5.txt
./rq_scan 2 rq_fields3.txt W 400 100 rq_out_R2_normW.txt -0.35 2 > rq_log_R2_normW.txt 2>&1
echo W2 >> rq_stage5.txt
./rq_scan 3 rq_fields3.txt T 500 150 rq_out_R3_normT.txt -0.35 3 > rq_log_R3_normT.txt 2>&1
echo T3 >> rq_stage5.txt
./rq_scan 2 rq_fields3.txt T 400 100 rq_out_R2_normT.txt -0.35 2 > rq_log_R2_normT.txt 2>&1
echo T2 >> rq_stage5.txt
echo DONE > rq_pass6.done
