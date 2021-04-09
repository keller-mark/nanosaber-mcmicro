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

## MCMICRO example data

Use the nextflow pipeline to download the example files:

```sh
mkdir -p ~/scratch/saberfish-mcmicro
nextflow run labsyspharm/mcmicro/exemplar.nf --name exemplar-001 --path ~/scratch/saberfish-mcmicro
```

Run the pipeline:
```sh
nextflow run labsyspharm/mcmicro --in ~/scratch/saberfish-mcmicro/exemplar-001 -w ~/scratch/saberfish-mcmicro/work -c ~/scratch/saberfish-mcmicro/mk596.config -with-report "~/scratch/saberfish-mcmicro/reports/$USER-$(date -Is).html"
```

