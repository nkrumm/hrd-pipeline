# hrd_nextflow
Install the nextflow binary in this directory
  
wget -qO- https://get.nextflow.io | bash
  
Execute locally, using docker images (must be available via docker locally)
  
./nextflow run main.nf -profile docker
