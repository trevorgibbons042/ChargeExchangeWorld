# D_track only: F8_TrackerFits.cpp prints pure / all reconstructed tracks, no pdfs. Submit: qsub -v RUN=<run folder in build/root>,SCRIPT=run_dtrack.sh Analysis.pbs
source ./run_lib.sh
export D_TRACK_ONLY=1
run 'F8_TrackerFits.cpp'
