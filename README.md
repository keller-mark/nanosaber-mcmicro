# nanosaber mcmicro example


## Installation on Mac

Install Java with JDK: https://adoptopenjdk.net/

```sh
# Install Nextflow
curl -s https://get.nextflow.io | bash
```

Move the Nextflow executable to a directory in your PATH.

Install bioformats: download the command-line tools from http://www.openmicroscopy.org/bio-formats/downloads/

Move the `bftools` directory to `~/software/bftools`

Install the Python environment:

```sh
conda env create -f environment.yml
conda activate nanosaber-mcmicro-env
```

Set the path to FIJI at the top of the Snakefile.

## Run the pipeline

```sh
snakemake -j 1
jupyter lab ./fc_all.ipynb
```

### Prepare for visualization with Vitessce

Open the registered 12-channel TIFF in ImageJ to convert the Z-stack to a channel stack.
- Image -> Hyperstacks -> Stack to Hyperstack
- 12 channel slices, 1 z slice, Composite display mode

Convert TIFF files to OME-TIFF using bioformats

```sh
bash ~/software/bftools/bfconvert ./data/lung_1_1/registration/lung_1_1.tif ./data/lung_1_1/registration/lung_1_1.ome.tiff

bash ~/software/bftools/bfconvert ./data/lung_1_1/segmentation/unmicst-lung_1_1/cellMask.tif ./data/lung_1_1/segmentation/unmicst-lung_1_1/cellMask.ome.tiff

tiffcomment ./data/lung_1_1/registration/lung_1_1.ome.tiff > ./data/lung_1_1/registration/lung_1_1.in.ome.xml
python add_channel_names.py -s lung_1_1
tiffcomment -set './data/lung_1_1/registration/lung_1_1.out.ome.xml' ./data/lung_1_1/registration/lung_1_1.ome.tiff
```

```sh
bash ~/software/bftools/bfconvert ./data/lung_2_1/registration/lung_2_1.ome.tif ./data/lung_2_1/registration/lung_2_1.ome.tiff

bash ~/software/bftools/bfconvert ./data/lung_2_1/segmentation/unmicst-lung_2_1/cellMask.tif ./data/lung_2_1/segmentation/unmicst-lung_2_1/cellMask.ome.tiff

tiffcomment ./data/lung_2_1/registration/lung_2_1.ome.tiff > ./data/lung_2_1/registration/lung_2_1.in.ome.xml
python add_channel_names.py -s lung_2_1
tiffcomment -set './data/lung_2_1/registration/lung_2_1.out.ome.xml' ./data/lung_2_1/registration/lung_2_1.ome.tiff
```

```sh
bash ~/software/bftools/bfconvert ./data/lung_2_2/registration/lung_2_2.tif ./data/lung_2_2/registration/lung_2_2.ome.tiff

bash ~/software/bftools/bfconvert ./data/lung_2_2/segmentation/unmicst-lung_2_2/cellMask.tif ./data/lung_2_2/segmentation/unmicst-lung_2_2/cellMask.ome.tiff

tiffcomment ./data/lung_2_2/registration/lung_2_2.ome.tiff > ./data/lung_2_2/registration/lung_2_2.in.ome.xml
python add_channel_names.py -s lung_2_2
tiffcomment -set './data/lung_2_2/registration/lung_2_2.out.ome.xml' ./data/lung_2_2/registration/lung_2_2.ome.tiff
```

## Vitessce instances

### Local

http://localhost:3000/?dataset=nanosaber-lung_1_1

http://localhost:3000/?dataset=nanosaber-lung_2_1

http://localhost:3000/?dataset=nanosaber-lung_2_2

### Cloud

Lung 1.1

```
https://s3.amazonaws.com/vitessce-data/demos/2021-05-24/a1802c9/index.html?url=data:,{"name":"lung_1_1","version":"1.0.1","description":"","public":true,"datasets":[{"uid":"lung_1_1","name":"lung_1_1","description":"lung_1_1","files":[{"type":"raster","fileType":"raster.json","options":{"schemaVersion":"0.0.2","usePhysicalSizeScaling":false,"images":[{"name":"Mask","type":"ome-tiff","url":"https://storage.googleapis.com/vitessce-demo-data/nanosaber-mcmicro/lung_1_1/segmentation/unmicst-lung_1_1/cellMask.ome.tiff","metadata":{"isBitmask":true,"transform":{"matrix":[1.0331,0,0,0,0,1.0331,0,0,0,0,1,0,0,0,0,1]}}},{"name":"Image","type":"ome-tiff","url":"https://storage.googleapis.com/vitessce-demo-data/nanosaber-mcmicro/lung_1_1/registration/lung_1_1.ome.tiff","metadata":{"isBitmask":false}}],"renderLayers":["Image","Mask"]}},{"type":"cell-sets","fileType":"cell-sets.json","url":"https://storage.googleapis.com/vitessce-demo-data/nanosaber-mcmicro/lung_1_1/flowcore/lung_1_1.cell-sets.json"},{"type":"cells","fileType":"cells.json","url":"https://storage.googleapis.com/vitessce-demo-data/nanosaber-mcmicro/lung_1_1/flowcore/lung_1_1.cells.json"}]}],"initStrategy":"auto","coordinationSpace":{},"layout":[{"component":"description","x":0,"y":0,"w":2,"h":1},{"component":"layerController","x":10,"y":0,"w":2,"h":4},{"component":"status","x":0,"y":1,"w":2,"h":1},{"component":"spatial","coordinationScopes":{},"x":2,"y":0,"w":8,"h":4},{"component":"cellSets","x":0,"y":2,"w":2,"h":2}]}
```

Lung 2.1

```
https://s3.amazonaws.com/vitessce-data/demos/2021-05-24/a1802c9/index.html?url=data:,{"name":"lung_2_1","version":"1.0.1","description":"","public":true,"datasets":[{"uid":"lung_2_1","name":"lung_2_1","description":"lung_2_1","files":[{"type":"raster","fileType":"raster.json","options":{"schemaVersion":"0.0.2","usePhysicalSizeScaling":false,"images":[{"name":"Mask","type":"ome-tiff","url":"https://storage.googleapis.com/vitessce-demo-data/nanosaber-mcmicro/lung_2_1/segmentation/unmicst-lung_2_1/cellMask.ome.tiff","metadata":{"isBitmask":true,"transform":{"matrix":[1.0234,0,0,0,0,1.0234,0,0,0,0,1,0,0,0,0,1]}}},{"name":"Image","type":"ome-tiff","url":"https://storage.googleapis.com/vitessce-demo-data/nanosaber-mcmicro/lung_2_1/registration/lung_2_1.ome.tiff","metadata":{"isBitmask":false}}],"renderLayers":["Image","Mask"]}},{"type":"cell-sets","fileType":"cell-sets.json","url":"https://storage.googleapis.com/vitessce-demo-data/nanosaber-mcmicro/lung_2_1/flowcore/lung_2_1.cell-sets.json"},{"type":"cells","fileType":"cells.json","url":"https://storage.googleapis.com/vitessce-demo-data/nanosaber-mcmicro/lung_2_1/flowcore/lung_2_1.cells.json"}]}],"initStrategy":"auto","coordinationSpace":{},"layout":[{"component":"description","x":0,"y":0,"w":2,"h":1},{"component":"layerController","x":10,"y":0,"w":2,"h":4},{"component":"status","x":0,"y":1,"w":2,"h":1},{"component":"spatial","coordinationScopes":{},"x":2,"y":0,"w":8,"h":4},{"component":"cellSets","x":0,"y":2,"w":2,"h":2}]}
```

Lung 2.2

```
https://s3.amazonaws.com/vitessce-data/demos/2021-05-24/a1802c9/index.html?url=data:,{"name":"lung_2_2","version":"1.0.1","description":"","public":true,"datasets":[{"uid":"lung_2_2","name":"lung_2_2","description":"lung_2_2","files":[{"type":"raster","fileType":"raster.json","options":{"schemaVersion":"0.0.2","usePhysicalSizeScaling":false,"images":[{"name":"Mask","type":"ome-tiff","url":"https://storage.googleapis.com/vitessce-demo-data/nanosaber-mcmicro/lung_2_2/segmentation/unmicst-lung_2_2/cellMask.ome.tiff","metadata":{"isBitmask":true,"transform":{"matrix":[1.167,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]}}},{"name":"Image","type":"ome-tiff","url":"https://storage.googleapis.com/vitessce-demo-data/nanosaber-mcmicro/lung_2_2/registration/lung_2_2.ome.tiff","metadata":{"isBitmask":false}}],"renderLayers":["Image","Mask"]}},{"type":"cell-sets","fileType":"cell-sets.json","url":"https://storage.googleapis.com/vitessce-demo-data/nanosaber-mcmicro/lung_2_2/flowcore/lung_2_2.cell-sets.json"},{"type":"cells","fileType":"cells.json","url":"https://storage.googleapis.com/vitessce-demo-data/nanosaber-mcmicro/lung_2_2/flowcore/lung_2_2.cells.json"}]}],"initStrategy":"auto","coordinationSpace":{},"layout":[{"component":"description","x":0,"y":0,"w":2,"h":1},{"component":"layerController","x":10,"y":0,"w":2,"h":4},{"component":"status","x":0,"y":1,"w":2,"h":1},{"component":"spatial","coordinationScopes":{},"x":2,"y":0,"w":8,"h":4},{"component":"cellSets","x":0,"y":2,"w":2,"h":2}]}
```


<!--

## Installation on O2

```sh
# Load Java
module load java/jdk-1.8u112
# Install Nextflow
curl -s https://get.nextflow.io | bash
# Move Nextflow executable to a directory in your PATH
mv nextflow ~/software/bin/nextflow
```

Check that the installed tools work:

```sh
which java
which nextflow
java -version
nextflow -version
```

## Run with example data on O2

Use the nextflow pipeline to download the example files:

```sh
mkdir -p ~/scratch/nanosaber-mcmicro
nextflow run labsyspharm/mcmicro/exemplar.nf --name exemplar-001 --path ~/scratch/nanosaber-mcmicro
```

Run the pipeline:
```sh
nextflow run labsyspharm/mcmicro --in ~/scratch/nanosaber-mcmicro/exemplar-001 -w ~/scratch/nanosaber-mcmicro/work -c ~/research/nanosaber-mcmicro/mk596.config -with-report "~/scratch/nanosaber-mcmicro/reports/$USER-$(date -Is).html"
```

## Troubleshoot singularity on O2

When using the pipeline on O2, you may run into this error: `singularity image is not in an allowed configured path`

To fix, locate the singularity image in `~/scratch/nanosaber-mcmicro/singularity/` and submit it for [automated testing](https://wiki.rc.hms.harvard.edu/display/O2/Running+Singularity+Containers+in+O2):

```sh
# Change to a singularity image extension recognized by O2
mv ~/scratch/nanosaber-mcmicro/singularity/labsyspharm-ashlar-1.13.0.img ~/scratch/nanosaber-mcmicro/singularity/labsyspharm-ashlar-1.13.0.sif
# Submit the image
module load csubmitter/latest
csubmitter --name labsyspharm-ashlar-1.13.0 --image-path /n/scratch3/users/m/mk596/nanosaber-mcmicro/singularity/labsyspharm-ashlar-1.13.0.sif
# Check the status
csubmitter --status
# When done processing, locate the verified image
ls /n/app/singularity/containers/mk596
```

To pull all of the other containers and submit:


```sh
cd ~/scratch/nanosaber-mcmicro/singularity

# illumination
singularity pull docker://labsyspharm/basic-illumination:1.0.1
mv basic-illumination_1.0.1.sif labsyspharm-basic-illumination-1.0.1.sif
csubmitter --name labsyspharm-basic-illumination-1.0.1 --image-path labsyspharm-basic-illumination-1.0.1.sif

# ashlar
singularity pull docker://labsyspharm/ashlar:1.13.0
mv ashlar_1.0.1.sif labsyspharm-ashlar-1.13.0.sif
csubmitter --name labsyspharm-ashlar-1.13.0 --image-path labsyspharm-ashlar-1.13.0.sif

# coreograph
singularity pull docker://labsyspharm/unetcoreograph:2.2.2
mv unetcoreograph_2.2.2.sif labsyspharm-unetcoreograph-2.2.2.sif
csubmitter --name labsyspharm-unetcoreograph-2.2.2 --image-path labsyspharm-unetcoreograph-2.2.2.sif
# or
singularity build labsyspharm-unetcoreograph-2.2.2.img docker://labsyspharm/unetcoreograph:2.2.2
mv labsyspharm-unetcoreograph-2.2.2.img labsyspharm-unetcoreograph-2.2.2.sif
csubmitter --name labsyspharm-unetcoreograph-2.2.2 --image-path labsyspharm-unetcoreograph-2.2.2.sif

# unmicst
singularity pull docker://labsyspharm/unmicst:2.6.10
mv unmicst_2.6.10.sif labsyspharm-unmicst-2.6.10.sif
csubmitter --name labsyspharm-unmicst-2.6.10 --image-path labsyspharm-unmicst-2.6.10.sif
# or
singularity build labsyspharm-unmicst-2.6.10.img docker://labsyspharm/unmicst:2.6.10
mv labsyspharm-unmicst-2.6.10.img labsyspharm-unmicst-2.6.10.sif
csubmitter --name labsyspharm-unmicst-2.6.10 --image-path labsyspharm-unmicst-2.6.10.sif

# cypository
singularity pull docker://labsyspharm/cypository:1.0.11
mv cypository_1.0.11.sif labsyspharm-cypository-1.0.11.sif
csubmitter --name labsyspharm-cypository-1.0.11 --image-path labsyspharm-cypository-1.0.11.sif

# ilastik
singularity pull docker://labsyspharm/mcmicro-ilastik:1.4.2
mv mcmicro-ilastik_1.4.2.sif labsyspharm-mcmicro-ilastik-1.4.2.sif
csubmitter --name labsyspharm-mcmicro-ilastik-1.4.2 --image-path labsyspharm-mcmicro-ilastik-1.4.2.sif

# s3seg
singularity pull docker://labsyspharm/s3segmenter:1.2.5
mv s3segmenter_1.2.5.sif labsyspharm-s3segmenter-1.2.5.sif
csubmitter --name labsyspharm-s3segmenter-1.2.5 --image-path labsyspharm-s3segmenter-1.2.5.sif

# quantification
singularity pull docker://labsyspharm/quantification:1.3.2
mv quantification_1.3.2.sif labsyspharm-quantification-1.3.2.sif
csubmitter --name labsyspharm-quantification-1.3.2 --image-path labsyspharm-quantification-1.3.2.sif

# naivestates
singularity pull docker://labsyspharm/naivestates:1.6.2
mv naivestates_1.6.2.sif labsyspharm-naivestates-1.6.2.sif
csubmitter --name naivestates-1.6.2 --image-uri docker://labsyspharm/naivestates:1.6.2
```


-->
