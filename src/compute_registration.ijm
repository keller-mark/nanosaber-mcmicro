sampleId = getArgument;

dataDir = "/Users/mkeller/research/dbmi/nanosaber-mcmicro/data";

open(dataDir + "/" + sampleId + "/pre-raw/" + sampleId + "-cycle-01-channel-01.tif");
run("16-bit");
open(dataDir + "/" + sampleId + "/pre-raw/" + sampleId + "-cycle-02-channel-01.tif");
run("16-bit");
open(dataDir + "/" + sampleId + "/pre-raw/" + sampleId + "-cycle-03-channel-01.tif");
run("16-bit");

targetImage = sampleId + "-cycle-02-channel-01.tif";



sourceImage = sampleId + "-cycle-03-channel-01.tif";
transFile = dataDir + "/" + sampleId + "/bunwarpj/" + sampleId + "-cycle-03-channel-01_direct_transf.txt";

registrationCommand = "source_image=" + sourceImage + " target_image=" + targetImage + " registration=Mono image_subsample_factor=0 initial_deformation=[Very Coarse] final_deformation=Fine divergence_weight=0 curl_weight=0 landmark_weight=0 image_weight=1 consistency_weight=10 stop_threshold=0.01 verbose save_transformations save_direct_transformation=" + transFile;

run("bUnwarpJ", registrationCommand);



sourceImage = sampleId + "-cycle-01-channel-01.tif";
transFile = dataDir + "/" + sampleId + "/bunwarpj/" + sampleId + "-cycle-01-channel-01_direct_transf.txt";

registrationCommand = "source_image=" + sourceImage + " target_image=" + targetImage + " registration=Mono image_subsample_factor=0 initial_deformation=[Very Coarse] final_deformation=Fine divergence_weight=0 curl_weight=0 landmark_weight=0 image_weight=1 consistency_weight=10 stop_threshold=0.01 verbose save_transformations save_direct_transformation=" + transFile;

run("bUnwarpJ", registrationCommand);