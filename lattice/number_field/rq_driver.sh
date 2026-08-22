#!/bin/bash
cd /home/ubuntu/code/math-modular-sources/lattice/number_field
set -x
./rq_scan 3 rq_fields.txt S 200 60 rq_out_R3_smallnorm.txt > rq_log_R3_smallnorm.txt 2>&1
./rq_scan 2 rq_fields.txt S 200 60 rq_out_R2_smallnorm.txt > rq_log_R2_smallnorm.txt 2>&1
./rq_scan 3 rq_fields.txt Q  60 40 rq_out_R3_ratval.txt   > rq_log_R3_ratval.txt   2>&1
./rq_scan 2 rq_fields.txt Q  60 40 rq_out_R2_ratval.txt   > rq_log_R2_ratval.txt   2>&1
./rq_scan 3 rq_fields.txt C 1500 600 rq_out_R3_wideMI.txt -0.35 3 > rq_log_R3_wideMI.txt 2>&1
./rq_scan 2 rq_fields.txt C 1500 600 rq_out_R2_wideMI.txt -0.35 2 > rq_log_R2_wideMI.txt 2>&1
./rq_scan 3 rq_fields.txt S  600 300 rq_out_R3_wideMIs.txt -0.35 3 > rq_log_R3_wideMIs.txt 2>&1
./rq_scan 2 rq_fields.txt S  600 300 rq_out_R2_wideMIs.txt -0.35 2 > rq_log_R2_wideMIs.txt 2>&1
echo DONE > rq_pass2.done
