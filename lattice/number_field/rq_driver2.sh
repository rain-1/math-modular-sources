#!/bin/bash
cd /home/ubuntu/code/math-modular-sources/lattice/number_field
./rq_scan 3 rq_fields2.txt U 3000 1000 rq_out_R3_normU.txt -0.35 3 > rq_log_R3_normU.txt 2>&1
echo S1 >> rq_stage.txt
./rq_scan 2 rq_fields2.txt U 3000 1000 rq_out_R2_normU.txt -0.35 2 > rq_log_R2_normU.txt 2>&1
echo S2 >> rq_stage.txt
./rq_scan 3 rq_fields2.txt T 1000 300 rq_out_R3_normT.txt -0.35 3 > rq_log_R3_normT.txt 2>&1
echo S3 >> rq_stage.txt
./rq_scan 2 rq_fields2.txt T 1000 300 rq_out_R2_normT.txt -0.35 2 > rq_log_R2_normT.txt 2>&1
echo S4 >> rq_stage.txt
echo DONE > rq_pass3.done
