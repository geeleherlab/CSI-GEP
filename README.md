# CSI-GEP
Consensus and Scalable Interpretations of Gene Expression Programs in Single Cells.\
We provide a bash file for submitting jobs on HPC and all scripts package dependencies can be found at CSI-GEP/.

## Getting started
New users: follow the **[end-to-end tutorial](TUTORIAL.md)**, from raw counts to final GEPs, with a simulated example
that runs in under an hour on a cluster.

## Docker
We provide Docker images that minimize complexity of configuration, and a step-by-step general instruction of running it on HPC.
1. On your local host, pull Docker images:
   ```
   docker pull ghcr.io/geeleherlab/csi-gep_py:2026.11
   docker pull ghcr.io/geeleherlab/csi-gep_r:2026.11
   ```
2. Below is an example on St. Jude HPC that utilizes singularity. On an HPC node:
   ```
   ## load module
   module load singularity/4.1.1
   ## remote login
   singularity registry login -u <my_userID> docker://svlprhpcreg01.stjude.org
   ```
3. On your local host:
   ```
   ## docker login
   docker login -u <my_userID> svlprhpcreg01.stjude.org
   ## tag images
   docker tag ghcr.io/geeleherlab/csi-gep_py:2026.11 svlprhpcreg01.stjude.org/hpcf/csi-gep_py:latest
   docker tag ghcr.io/geeleherlab/csi-gep_r:2026.11 svlprhpcreg01.stjude.org/hpcf/csi-gep_r:latest
   ## push images
   docker image push svlprhpcreg01.stjude.org/hpcf/csi-gep_py:latest
   docker image push svlprhpcreg01.stjude.org/hpcf/csi-gep_r:latest
   ```
4. From HPC, pull images and submit jobs:
   ```
   ## pull
   singularity pull docker://svlprhpcreg01.stjude.org/hpcf/csi-gep_py:latest
   singularity pull docker://svlprhpcreg01.stjude.org/hpcf/csi-gep_r:latest
   ```
An example code of executing the Dockerized CSI-GEP *csigep_submit.bsub* can be found at Docker/.

## Changelog                                                                  
                                                                                
  ### Security updates to Python dependencies (2026-11-23)

  **Changed:** `CSI-GEP/requirements_py.txt` and `CSI-GEP/environment_py.yml` now pin Python 3.12 and patched versions of
  the Python packages (torch 2.13.0, numpy 1.26.4, scikit-learn 1.5.2, scanpy 1.10.3, anndata 0.10.9), resolving the
  Dependabot alerts for the old pins (Python 3.7, torch 1.13.1). The torch wheels bundle the CUDA 13.0 runtime and need
  a recent NVIDIA driver. The tutorial demo recovers all 21 simulated programs with the new environment. The container
  images were rebuilt with the same versions from `Docker/Dockerfile_py` and `Docker/Dockerfile_R` (R 4.4.3) by the
  `docker-images` GitHub Actions workflow in geeleherlab/CSI-GEP: `ghcr.io/geeleherlab/csi-gep_py:2026.11` and `csi-gep_r:2026.11` (also `:latest`).
  The previous images remain available by digest.

  ### Tutorial and conda environments (2026-10-22)

  **Added:** an end-to-end [tutorial](TUTORIAL.md) with a simulated example (`tutorial/`), and conda environment
  files (`CSI-GEP/environment_py.yml`, `CSI-GEP/environment_R.yml`) for running the scripts without containers.

  **Fixed:** `CSI-GEP/csigep_submit.bsub` passed the rescue argument as `FALSE`, which `GEP_analysis.R` converts to
  `NA`, stopping with "missing value where TRUE/FALSE needed". It now passes `rescue=0` (set to `1` to enable rescue),
  as `Docker/csigep_submit.bsub` does. The Docker images are unchanged.

  ### Docker Image Update (2026-02-10)                                          
                                                                                
  **Fixed:** Boolean argument parsing bug in `GEP_analysis.R`                   
                                                                                
  The `check_for_rescue` parameter (argument 5) was not being properly converted
   from string to logical, causing the rescue logic to fail silently.           
                                                                                
  **Fix:** Changed argument parsing to explicitly convert:                      
  ```r                                                                          
  # Before (buggy)                                                              
  check_for_rescue = args[5]                                                    
                                                                                
  # After (fixed)                                                               
  check_for_rescue = as.logical(as.numeric(args[5]))                            
   ```                                                                             
  Docker image tags:                                                            
  - ghcr.io/geeleherlab/csi-gep_r:latest - Fixed version (recommended)          
  - ghcr.io/geeleherlab/csi-gep_r:v1-deprecated - Original version with bug


## Benchmarking 
We show a complete example of benchmarking analysis on a simulated dataset with 250,000 cells and 25,000 genes. The simulated data is stored at simulated_data/.\
Insert figure here.\
Details can be found in Readme.md at benchmarking/.


