sampleId = getArgument();

IJ.log("Sample ID: " + sampleId);

dataDir = "/Users/mkeller/research/dbmi/nanosaber-mcmicro/data"

cycle01 = IJ.openImage(dataDir + "/" + sampleId + "/raw/" + sampleId + "-cycle-01.tif");
cycle02 = IJ.openImage(dataDir + "/" + sampleId + "/raw/" + sampleId + "-cycle-02.tif");
cycle03 = IJ.openImage(dataDir + "/" + sampleId + "/raw/" + sampleId + "-cycle-03.tif");

imp = Concatenator.run(cycle01, cycle02, cycle03);

numChannels = (sampleId === "lung_1_1" ? 11 : 12);

imp2 = HyperStackConverter.toHyperStack(imp, numChannels, 1, 1, "Composite");

IJ.saveAs(imp2, "Tiff", dataDir + "/" + sampleId + "/registration/" + sampleId + ".ome.tif");