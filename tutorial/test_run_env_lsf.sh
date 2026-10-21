#!/bin/bash
# Test harness for TUTORIAL.md Option B (conda environments + the repo's scripts) on an LSF cluster.
# Input is the Docker-route run's simulated data (tutorial_run/data). Job 1 builds both environments from CSI-GEP/environment_{py,R}.yml (+ the two GitHub R packages), checks the
# imports, and records whether the repo scripts are identical to those in the containers; then the demo runs
# with `python CSI-GEP/*.py` / `Rscript CSI-GEP/*.R`, one job per step and per K, as in test_run_lsf.sh.
# Run from the repo root on a login node:  GPU_Q=gpu_priority bash tutorial/test_run_env_lsf.sh
# Out: tutorial_run_envB/ (envs/, data, CSI-GEP output, logs/, env_check.txt, test_summary.txt with STATUS line)
set -euo pipefail
REPO=$(pwd)
RUN=$REPO/tutorial_run_envB/                 # trailing slash required by the R scripts
ENV=$RUN/envs
CONT=$REPO/tutorial_run/containers           # only used to compare scripts with the container copies
NAME=demo
GPU_Q=${GPU_Q:-gpu}; CPU_Q=${CPU_Q:-compbio}
TAG=csigB$$
mkdir -p $RUN/data $RUN/logs $RUN/results $RUN/GEP_results $ENV
L=$RUN/logs
PY="$ENV/csigep_py/bin/python"
RR="$ENV/csigep_r/bin/Rscript"
KS=$(cat $REPO/tutorial/K.txt)
KCSV=$(paste -sd, $REPO/tutorial/K.txt)

sub() {  # sub <name> <queue> <mem_GB> <cores> <dependency or ""> <extra bsub args> <command>
  local dep=(); [ -n "$5" ] && dep=(-w "$5")
  bsub -P csigep_tutorial -J "${TAG}_$1" -q "$2" -n "$4" -R "rusage[mem=${3}GB]" -M ${3}GB ${dep[@]+"${dep[@]}"} $6 \
       -o "$L/$1_%J.out" -e "$L/$1_%J.err" "bash -c '$7'" > /dev/null
}

# 0. environments (skipped when already built)
sub envs $CPU_Q 16 4 "" "-W 240" "module load conda3/202402 && cd $REPO && \
 { [ -x $PY ] || conda env create -q -p $ENV/csigep_py -f CSI-GEP/environment_py.yml; } && \
 { [ -x $RR ] || conda env create -q -p $ENV/csigep_r -f CSI-GEP/environment_R.yml; } && \
 PATH=$ENV/csigep_r/bin:\$PATH $RR -e \"if (!all(c(\\\"jaccard\\\",\\\"elbow\\\") %in% rownames(installed.packages()))) { remotes::install_github(\\\"ncchung/jaccard@ad1d8fc\\\", upgrade=\\\"never\\\"); remotes::install_github(\\\"ahasverus/elbow@cf58784\\\", upgrade=\\\"never\\\") }\" && \
 $PY -c \"import torch, torchnmf, scanpy, anndata, fastcluster, sklearn; print(\\\"python ok\\\", torch.__version__)\" > $RUN/env_check.txt && \
 $RR -e \"suppressMessages({library(jaccard); library(igraph); library(vroom); library(elbow); library(foreach); library(doParallel); library(dplyr); library(ggplot2)}); cat(\\\"R ok\\\", R.version.string, \\\"\\\\n\\\")\" >> $RUN/env_check.txt && \
 module load singularity/4.2.0 && for f in prepare.py cNMF.py; do singularity exec $CONT/csi-gep_py_latest.sif cat /app/\$f | cmp -s - CSI-GEP/\$f && echo \"\$f identical to container\" || echo \"\$f DIFFERS from container\"; done >> $RUN/env_check.txt && \
 for f in Jaccard.R GEP_analysis.R; do singularity exec $CONT/csi-gep_r_latest.sif cat /app/\$f | cmp -s - CSI-GEP/\$f && echo \"\$f identical to container\" || echo \"\$f DIFFERS from container\"; done >> $RUN/env_check.txt"
# 1. input: the same simulated data as the Docker-route run (tutorial_run/data), so both routes get identical input
[ -s $REPO/tutorial_run/data/demo_counts.h5ad ] || { echo "run tutorial/test_run_lsf.sh first (needs tutorial_run/data)"; exit 1; }
ln -sf $REPO/tutorial_run/data/demo_counts.h5ad $REPO/tutorial_run/data/ground_truth_programs.tsv $RUN/data/
# 2. prepare
sub prepare $CPU_Q 32 1 "done(${TAG}_envs)" "" \
  "cd $REPO && $PY CSI-GEP/prepare.py --numgenes 2000 -c $RUN/data/demo_counts.h5ad --output-dir $RUN"
# 3. factorize + consensus, one job per K
for k in $KS; do
  sub nmf_$k $GPU_Q 16 1 "done(${TAG}_prepare)" "-gpu num=1" \
    "cd $REPO && $PY CSI-GEP/cNMF.py --output-dir $RUN --components $k --n-iter 100"
done
# 4. Jaccard, one job per K
for k in $KS; do
  sub jac_$k $CPU_Q 32 16 "done(${TAG}_nmf_$k)" "" "cd $REPO && $RR --vanilla CSI-GEP/Jaccard.R $k $RUN"
done
# 5. GEP analysis (rescue off = 0)
sub gep $CPU_Q 16 1 "done(${TAG}_jac_*)" "" \
  "cd $REPO && $RR --vanilla CSI-GEP/GEP_analysis.R $RUN $NAME ${RUN}GEP_results/ $KCSV 0"
# 6. ground-truth check + summary
sub check $CPU_Q 8 1 "done(${TAG}_gep)" "" \
  "$RR --vanilla $REPO/tutorial/check_ground_truth.R ${RUN}GEP_results/${NAME}_community_geps.RDS $RUN/data/ground_truth_programs.tsv ${RUN}GEP_results/ground_truth_check.tsv > $RUN/check.txt && { cat $RUN/env_check.txt ${RUN}GEP_results/output_summary.txt $RUN/check.txt; } > $RUN/test_summary.txt && echo \"STATUS: OK | \$(grep -h \"GEPs found\" $RUN/check.txt)\" >> $RUN/test_summary.txt"
echo "Submitted ${TAG}_* jobs ($(echo $KS | wc -w) K values). Watch: bjobs -J '${TAG}_*'   Result: $RUN/test_summary.txt"
