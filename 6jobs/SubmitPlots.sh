#!/bin/bash
# Submits the plotting chain for one run folder as PBS jobs (Plot.pbs), the lanes of 5rootMacros/plot.sh in parallel, instead of running plot.sh in one go.
# Quick redo of one cut: CUT_ONLY=megacut3 bash SubmitPlots.sh <run> "" CutPlots 6_cuts/9_megaCut3 rewrites only the cut_megacut3* files (CutPlots.cpp) and then replots only that folder.
# Usage: bash SubmitPlots.sh <folder in build/root> [job ids CutPlots waits for, colon separated] [first job] [plots/ folder to redo], e.g. bash SubmitPlots.sh Physics_xs1e10
# The Pdfs job is only submitted with PDFS=1 in the environment or with Pdfs as the first job (the examples below that mention it need one of them).
# first job: CutPlots (default), Cuts, Fom, CEblock, SecCE or Pdfs. The jobs before it are not submitted (their files are already in cuts/), so e.g. "Cuts" runs Cuts, Fom, CEblock, SecCE and Pdfs: bash SubmitPlots.sh Physics_801 "" Cuts
# plots/ folder to redo: one folder of plots/ (a leading plots/ is dropped), e.g. bash SubmitPlots.sh Physics_801 "" CutPlots 6_cuts/3_ResCuts. It is deleted, then only the jobs that write it run (none of CutPlots), then Pdfs.
#   1_position, 2_energy          Plot3D part 1
#   6_cuts, 6_cuts/<cut folder>   Plot3D parts 2, 3, 4 of that cut (all cuts for 6_cuts), part 5 (rank csv files), QualityStudy for 3_ResCuts, 4_nParticles, 5_nHits
#   3_ReconCalcVals, 5_fomCalc    the Fom lane      4_CEblock  the CEblock lane      7_secondaryCE  the SecCE lane
# Needs cut_*.root and settings.csv already in the run folder (CutEvents.cpp, ReconAngleAll.cpp, CalcAPandT.cpp, CountCE.cpp done), or the job ids of the jobs that write them (SubmitEvents.sh).
RUN=$1
ORDER="CutPlots Cuts Fom CEblock SecCE Pdfs"
# the Pdfs job (tidy_plots.py, cut_side_box.py, crop_pdfs.py) is not submitted unless it is asked for: PDFS=1 in the environment, or Pdfs as the first job
[ "${3:-}" = Pdfs ] && PDFS=1
FROM=${3:-CutPlots}
TARGET=$4
[ -n "$RUN" ] || { echo "usage: bash SubmitPlots.sh <run folder> [job ids to wait for] [first job: $ORDER] [plots/ folder to redo]"; exit 1; }
case " $ORDER " in *" $FROM "*) ;; *) echo "unknown first job $FROM (one of: $ORDER)"; exit 1;; esac
cd "$(dirname "$0")" || exit 1
source ./submit_lib.sh

# cut_indexes <folder under 6_cuts, empty for all>: the P3_CUT numbers of Plot3D.cpp (0: 0_cuts, 1_hitZDC and 2_KECUT80, then sampleCuts in order). Keep in step with Plot3D.cpp and Plot.pbs.
cut_indexes() {
    local n=${1#*/}
    case $1 in
        "") echo 0 1 2 3 4 5 6 7 8 9 10 11 ;;
        0_cuts|0_cuts/*|1_hitZDC|1_hitZDC/*|2_KECUT80|2_KECUT80/*) echo 0 ;;
        3_ResCuts|3_ResCuts/*) echo 1 ;;
        4_nParticles) echo 2 3 ;;
        4_nParticles/1_*) echo 2 ;;
        4_nParticles/2_*) echo 3 ;;
        5_nHits|5_nHits/*) echo 4 ;;
        6_megaCut) echo 5 6 ;;
        6_megaCut/[1-2]_*) echo $((4 + ${n:0:1})) ;;
        7_megaCut2) echo 7 8 ;;
        7_megaCut2/[1-2]_*) echo $((6 + ${n:0:1})) ;;
        8_nZdcParticles|8_nZdcParticles/*) echo 9 ;;
        9_megaCut3) echo 10 11 ;;
        9_megaCut3/[1-2]_*) echo $((9 + ${n:0:1})) ;;
        *) return 1 ;;
    esac
}

PARTS= LANE= QS=
if [ -n "$TARGET" ]; then
    TARGET=${TARGET#plots/}; TARGET=${TARGET%/}
    case "$TARGET" in ""|.*|*..*|/*) echo "bad folder: $4"; exit 1;; esac
    top=${TARGET%%/*}; rest=; [ "$top" != "$TARGET" ] && rest=${TARGET#*/}
    case $top in
        1_position|2_energy) PARTS=1::0+1::1+1::2+1::3 ;;
        6_cuts)
            idx=$(cut_indexes "$rest") || { echo "unknown cut folder: 6_cuts/$rest"; exit 1; }
            for i in $idx; do
                if [ "$i" = 0 ]; then PARTS+="2:0+3:0:0+3:0:1+3:0:2+3:0:3+4:0:0+4:0:1+4:0:2+4:0:3+" # cut 0 holds four folders (P3_SUB), a process each
                else PARTS+="2:$i+3:$i+4:$i+"; fi
            done
            [ -z "$rest" ] && PARTS+="4:12+4:13+4:14+" # all of 6_cuts: also the KE scan sets of Plot3D.cpp part 4 (cuts 12 to 14)
            case " $idx " in *" 0 "*) PARTS+="1::3+" ;; esac # 1_hitZDC, 0_cuts and 2_KECUT80 also hold plots of part 1
            PARTS=${PARTS%+}
            case $rest in ""|3_ResCuts*|4_nParticles*|5_nHits*) QS=1 ;; esac ;;
        3_ReconCalcVals|5_fomCalc) LANE=fom ;;
        4_CEblock) LANE=ceblock ;;
        7_secondaryCE) LANE=secce ;;
        *) echo "unknown plots folder $top (one of: 1_position 2_energy 3_ReconCalcVals 4_CEblock 5_fomCalc 6_cuts 7_secondaryCE)"; exit 1 ;;
    esac
    [ -d "../build/root/$RUN/plots/$TARGET" ] || { echo "no folder ../build/root/$RUN/plots/$TARGET"; exit 1; }
    rm -rf "../build/root/$RUN/plots/$TARGET" || exit 1 # the old plots go, so nothing stale stays beside the new ones
fi

rsync -a --exclude __pycache__ ../5rootMacros/ "../build/root/$RUN/" || exit 1 # the jobs run the run folder's copy of the macros
# a run from the start (no folder, first job CutPlots) rebuilds every plot: the old plots/ and the Plot3D logs go first, so nothing of an older macro version stays beside the new files. cuts/ stays.
if [ -z "$TARGET" ] && [ "$FROM" = CutPlots ]; then
    rm -rf "../build/root/$RUN/plots"; rm -f ../build/root/$RUN/Plot3D_part*.log
fi

if [ -n "$TARGET" ]; then # one folder: its jobs only (the cuts/ files are already there), then the python steps
    ids=
    [ -n "$PARTS" ] && { p3=$(PARTS=$PARTS submit Redo plot3d "$2" 128); ids=$p3; }
    [ -n "$QS" ] && ids="$ids $(PARTS= DEPKIND=afterany submit RedoQS quality "$p3" 16)" # PARTS only goes to the plot3d job
    [ -n "$LANE" ] && ids=$(PARTS= submit Redo $LANE "$2")
    if [ -n "$PDFS" ]; then PARTS= DEPKIND=afterany submit Pdfs pdfs "$(echo $ids | tr ' ' ':')"; fi
    exit
fi

# on <job>: true when the job is FROM or comes after it in ORDER
on() { [ "$1" = "$FROM" ] || case " $ORDER " in *" $FROM "*"$1 "*) return 0;; *) return 1;; esac; }

# CutPlots (KE cut lists, ZDC proton cut files, sample cut files) is on the critical path, ahead of every lane: all the cores
on CutPlots && cut=$(submit CutPlots cuttrees "$2" 128)
# the lanes start after CutPlots has written the cuts/ files; without CutPlots they wait for the ids in $2, if any
w=${cut:-$2}
on Cuts && cuts=$(submit Cuts cuts "$w" 128) # Plot3D runs nine processes side by side (Plot.pbs): a whole node
on Fom && fom=$(submit Fom fom "$w")
on CEblock && ce=$(submit CEblock ceblock "$w")
on SecCE && sec=$(submit SecCE secce "$w")
# tidy, side box and crop once every lane has ended, so the pdfs of the lanes that worked are still finished when one failed
ids=$(echo $cuts $fom $ce $sec | tr ' ' ':')
if [ -n "$PDFS" ]; then DEPKIND=afterany submit Pdfs pdfs "${ids:-$w}"; fi
