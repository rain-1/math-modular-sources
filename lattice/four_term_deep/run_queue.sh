#!/bin/bash
# ./run_queue.sh JOBFILE NJOBS   -- runs the job list with NJOBS in parallel, in order
cd "$(dirname "$0")"; mkdir -p out logs
JF=$1; NJ=${2:-8}
xargs -a "$JF" -d '\n' -P "$NJ" -I{} bash -c '{}'
echo "QUEUE DONE $JF" >> logs/queue.done
