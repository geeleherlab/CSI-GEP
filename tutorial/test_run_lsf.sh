#!/bin/bash
# Test harness for TUTORIAL.md on an LSF cluster: runs every tutorial step as its own job, with the same
# commands as the tutorial, so the numbers quoted in the tutorial come from a real run.
# Run from the repo root on a login node:  bash tutorial/test_run_lsf.sh
# Out: tutorial_run/ (data, CSI-GEP output, logs/, test_summary.txt with a STATUS line)
set -euo pipefail
REPO=$(pwd)
RUN=$REPO/tutorial_run/                      # trailing slash required by the R scripts
NAME=demo
GPU_Q=${GPU_Q:-gpu}                           # GPU queue for NMF; GPU_Q=none runs NMF on CPU (CPU_Q, 8 cores)
CPU_Q=${CPU_Q:-compbio}
TAG=csigt$$                                   # unique job-name prefix for dependencies
mkdir -p $RUN/data $RUN/logs $RUN/results $RUN/GEP_results $RUN/containers
L=$RUN/logs
SING="module load singularity/4.2.0 &&"
PY="singularity exec --nv $RUN/containers/csi-gep_py_latest.sif"
RR="singularity exec $RUN/containers/csi-gep_r_latest.sif"
KS=$(cat $REPO/tutorial/K.txt)
KCSV=$(paste -sd, $REPO/tutorial/K.txt)

sub() {  # sub <name> <queue> <mem_GB> <cores> <dependency or ""> <extra bsub args> <command>
  local dep=(); [ -n "$5" ] && dep=(-w "$5")
  bsub -P csigep_tutorial -J "${TAG}_$1" -q "$2" -n "$4" -R "rusage[mem=${3}GB]" -M ${3}GB ${dep[@]+"${dep[@]}"} $6 \
       -o "$L/$1_%J.out" -e "$L/$1_%J.err" "bash -c '$7'" > /dev/null
}

# 1-2. containers (skipped when already present)
sub get_py $CPU_Q 8 1 "" "" "$SING cd $RUN/containers && { [ -s csi-gep_py_latest.sif ] || singularity pull docker://ghcr.io/geeleherlab/csi-gep_py:latest; }"
sub get_r  $CPU_Q 8 1 "" "" "$SING cd $RUN/containers && { [ -s csi-gep_r_latest.sif ] || singularity pull docker://ghcr.io/geeleherlab/csi-gep_r:latest; }"
# 3. simulate the demo data (no download needed)
sub simulate $CPU_Q 32 1 "done(${TAG}_get_py)" "" \
  "$SING $PY python $REPO/tutorial/simulate_demo.py $RUN/data"
# 4. prepare
sub prepare $CPU_Q 32 1 "done(${TAG}_simulate)" "" \
  "$SING $PY python /app/prepare.py --numgenes 2000 -c $RUN/data/demo_counts.h5ad --output-dir $RUN"
# 5. factorize + consensus, one job per K
for k in $KS; do
  if [ "$GPU_Q" == "none" ]; then q=$CPU_Q; n=8; g=""; else q=$GPU_Q; n=1; g="-gpu num=1"; fi
  sub nmf_$k $q 16 $n "done(${TAG}_prepare)" "$g" \
    "$SING $PY python /app/cNMF.py --output-dir $RUN --components $k --n-iter 100"
done
# 6. Jaccard, one job per K
for k in $KS; do
  sub jac_$k $CPU_Q 32 16 "done(${TAG}_nmf_$k) && done(${TAG}_get_r)" "" \
    "$SING $RR Rscript --vanilla /app/Jaccard.R $k $RUN"
done
# 7. GEP analysis (rescue off = 0)
sub gep $CPU_Q 16 1 "done(${TAG}_jac_*)" "" \
  "$SING $RR Rscript --vanilla /app/GEP_analysis.R $RUN $NAME ${RUN}GEP_results/ $KCSV 0"
# 8. ground-truth check + summary
sub check $CPU_Q 8 1 "done(${TAG}_gep)" "" \
  "$SING $RR Rscript --vanilla $REPO/tutorial/check_ground_truth.R ${RUN}GEP_results/${NAME}_community_geps.RDS $RUN/data/ground_truth_programs.tsv ${RUN}GEP_results/ground_truth_check.tsv > $RUN/check.txt && { cat ${RUN}GEP_results/output_summary.txt; cat $RUN/check.txt; } > $RUN/test_summary.txt && echo \"STATUS: OK | \$(grep -h \"GEPs found\" $RUN/check.txt)\" >> $RUN/test_summary.txt"
echo "Submitted ${TAG}_* jobs ($(echo $KS | wc -w) K values). Watch: bjobs -J '${TAG}_*'   Result: $RUN/test_summary.txt"
