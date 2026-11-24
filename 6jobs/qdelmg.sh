#!/bin/bash -l
for j in $(qstat -u tgibbons | grep F9Mg_ | awk '{print $1}'); do qdel $j; done
sleep 5
bash /lcrc/project/eic/tgibbons/f9submit3.sh | wc -l
