# saberfish mcmicro example


## Installation on O2

```sh
# Load Java
module load java/jdk-1.8u112
# Install NextFlow
curl -s https://get.nextflow.io | bash
# Move nextflow executable to a directory in your PATH
mv nextflow ~/software/bin/nextflow
```

Check that the installed tools work:

```sh
which java
which nextflow
java -version
nextflow -version
```

## Run with example data

Use the nextflow pipeline to download the example files:

```sh
mkdir -p ~/scratch/saberfish-mcmicro
nextflow run labsyspharm/mcmicro/exemplar.nf --name exemplar-001 --path ~/scratch/saberfish-mcmicro
```

Run the pipeline:
```sh
nextflow run labsyspharm/mcmicro --in ~/scratch/saberfish-mcmicro/exemplar-001 -w ~/scratch/saberfish-mcmicro/work -c ~/research/saberfish-mcmicro/mk596.config -with-report "~/scratch/saberfish-mcmicro/reports/$USER-$(date -Is).html"
```

## Troubleshoot singularity on O2

When using the pipeline on O2, you may run into this error: `singularity image is not in an allowed configured path`

To fix, locate the singularity image in `~/scratch/saberfish-mcmicro/singularity/` and submit it for [automated testing](https://wiki.rc.hms.harvard.edu/display/O2/Running+Singularity+Containers+in+O2):

```sh
# Change to a singularity image extension recognized by O2
mv ~/scratch/saberfish-mcmicro/singularity/labsyspharm-ashlar-1.13.0.img ~/scratch/saberfish-mcmicro/singularity/labsyspharm-ashlar-1.13.0.sif
# Submit the image
module load csubmitter/latest
csubmitter --name labsyspharm-ashlar-1.13.0 --image-path /n/scratch3/users/m/mk596/saberfish-mcmicro/singularity/labsyspharm-ashlar-1.13.0.sif
# Check the status
csubmitter --status
# When done processing, locate the verified image
ls /n/app/singularity/containers/mk596
```

Now change the container path in the Nextflow config

- ashlar
```
# from
container = "docker://labsyspharm/ashlar:${params.ashlarVersion}"
# to
container = "file:///n/app/singularity/containers/mk596/labsyspharm-ashlar-1.13.0.sif"
```

- unmicst
```
# from
container = "docker://labsyspharm/unmicst:${params.unmicstVersion}"
# to
container = "file:///n/app/singularity/containers/mk596/labsyspharm-unmicst-2.6.10.sif"
```
