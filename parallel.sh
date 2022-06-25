#!/bin/bash

# This is the concurrency limit
MAX_POOL_SIZE=$(grep -c ^processor /proc/cpuinfo)

# Output file
OUTPUT=output.txt

# This is used within the program. Do not change.
CURRENT_POOL_SIZE=0

#Job inputs
Lx="100"
Ly="100"
dl="0.5"
tf="200"
dt="1"
corr_len="2"
V0="0.1"
tot_frames="100"

#Total jobs
TOTAL_JOBS=120

# This is a just a function to print the output as a log with timestamp
_log() {
        echo " $(date +'[%F %T]') - $1"
}

# This is the custom function to process each job read from the file
process_job() {
  # customize your job function as required
  # in our example, we just "ping" each hostname read from the file
  ./RandPot.o $Lx $Ly $dl $tf $dt $corr_len $V0 $1 $tot_frames $2
}

computeMean(){
  python3 statAnalysis.py $1 $2 $3 $4 $5
}

# ------ This is the main program code --------
# Starting the timer
T1=$(date +%s)

i="1"
while [ $i -le $TOTAL_JOBS ]; do
  
  # This is the blocking loop where it makes the program to wait if the job pool is full
  while [ $CURRENT_POOL_SIZE -ge $MAX_POOL_SIZE ]; do
    CURRENT_POOL_SIZE=$(jobs | wc -l)
  done
  
  
  # This is a custom function to process each job
  rand=$(tr -cd "[:digit:]" < /dev/urandom | head -c 8)
  _log "Starting job $i with seed $rand"
  process_job $rand $i &
  
  # When a new job is created, the program updates the $CURRENT_POOL_SIZE variable before next iteration
  CURRENT_POOL_SIZE=$(jobs | wc -l)
  #_log "Current pool size = $CURRENT_POOL_SIZE"
  ((i++))
done  # this is where we feed the $JOB_LIST file for the read operation

# wait for all background jobs (forks) to exit before exiting the parent process
wait

# Ending the timer
T2=$(date +%s)

_log "All jobs completed in $((T2-T1)) seconds. Parent process exiting."
_log "Computing mean..."
computeMean $Lx $Ly $dl $V0 $corr_len
_log "Plotting mean..."
python3 wavefunc_ev.py $Lx $Ly $dl $V0 $corr_len
_log "Final output is written in $OUTPUT"
exit 0