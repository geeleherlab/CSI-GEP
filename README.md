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
   docker pull ghcr.io/geeleherlab/csi-gep_py:latest
   docker pull ghcr.io/geeleherlab/csi-gep_r:latest
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
   docker tag ghcr.io/geeleherlab/csi-gep_py:latest svlprhpcreg01.stjude.org/hpcf/csi-gep_py:latest
   docker tag ghcr.io/geeleherlab/csi-gep_r:latest svlprhpcreg01.stjude.org/hpcf/csi-gep_r:latest
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


