#!/bin/bash
cd /home/ubuntu/code/math-modular-sources/lattice/number_field
./rq_scan 3 rq_fields2.txt U 2000 400 rq_out_R3_normU.txt -0.35 3 > rq_log_R3_normU.txt 2>&1
echo T1 >> rq_stage3.txt
./rq_scan 2 rq_fields2.txt U 2000 400 rq_out_R2_normU.txt -0.35 2 > rq_log_R2_normU.txt 2>&1
echo T2 >> rq_stage3.txt
./rq_scan 3 rq_fields2.txt T 800 200 rq_out_R3_normT.txt -0.35 3 > rq_log_R3_normT.txt 2>&1
echo T3 >> rq_stage3.txt
./rq_scan 2 rq_fields2.txt T 800 200 rq_out_R2_normT.txt -0.35 2 > rq_log_R2_normT.txt 2>&1
echo T4 >> rq_stage3.txt
echo DONE > rq_pass4.done
