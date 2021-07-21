
arg = getArgument();
args = arg.split("__");

sampleId = args[0];
cycleId = args[1];

IJ.log("Sample ID: " + sampleId);
IJ.log("Cycle ID: " + cycleId);

dataDir = "/Users/mkeller/research/dbmi/nanosaber-mcmicro/data"

registeredSuffix = (cycleId === "cycle-02" ? "" : "-registered");

imp1 = IJ.openImage(dataDir + "/" + sampleId + "/pre-raw/" + sampleId + "-" + cycleId + "-channel-01" + registeredSuffix + ".tif");
imp2 = IJ.openImage(dataDir + "/" + sampleId + "/pre-raw/" + sampleId + "-" + cycleId + "-channel-02" + registeredSuffix + ".tif");
imp3 = IJ.openImage(dataDir + "/" + sampleId + "/pre-raw/" + sampleId + "-" + cycleId + "-channel-03" + registeredSuffix + ".tif");

arrayOfImages = [imp1, imp2, imp3];

// Cycle 1 of Sample 1.1 only has 3 channels available.
if(!(sampleId === "lung_1_1" && cycleId === "cycle-01")) {
  imp4 = IJ.openImage(dataDir + "/" + sampleId + "/pre-raw/" + sampleId + "-" + cycleId + "-channel-04" + registeredSuffix + ".tif");
  arrayOfImages.push(imp4);
}

for (var i = 0; i < arrayOfImages.length; i++){
  IJ.run(arrayOfImages[i], "32-bit", "");
  IJ.run(arrayOfImages[i], "16-bit", "");
}


imp = ImagesToStack.run(arrayOfImages);
IJ.run(imp, "16-bit", "");
IJ.saveAs(imp, "Tiff", dataDir + "/" + sampleId + "/raw/" + sampleId + "-" + cycleId + ".tif");

IJ.log("Done stacking for " + sampleId + " " + cycleId);
