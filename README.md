# HRD Testing Project

Important links:
 - Tracker issue in genetics/projects repo: https://gitlab.labmed.uw.edu/genetics/projects/issues/100
 - Current research pipeline: https://bitbucket.org/nithishak/loh-pipeline/src/master/


## Pipeline installation

Install the nextflow binary in this directory
  
```bash
wget -qO- https://get.nextflow.io | bash
```
  
Execute locally, using the singularity images
  
```bash
./nextflow run main.nf -profile singularity
```
