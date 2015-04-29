/* c1 = phase; c2 = fluor; c3 = nucl marker; c1_pha = threshold
 *  
 * Manual steps: 
 * 0) Place raw image stacks in data folder; preprocess.jim in folder root
 * 1) Measure angle to be rotated (pos is CW, neg is CCW), set rot var
 * 2) Choose range of frames to keep, set start, end vars
 * 3) If images need to be flipped vertically so traps face "up", set reverse = 1
 * 
 * Macro function: 
 * 1) Make substack of frames in desired range for all 3 channel stacks
 * 2) Measure angle in phase, apply transformation to all 3 channels
 * 3) Create threshold image from phase
 * 4) Substact background from 2x flu channel stacks
 * 5) Invert so the traps are facing "up" if necessary
 * 6) Export 4 stacks; export each stack as individual .tif images
 * 
 * (optionally) Delete the artifacts above and below traps in the threshold stack
 */


 //TO DO: STREAMLINE BETTER AS FUNCTIONS; would allow for setting of variables
 //for every position and processing all of the image stack positions at once


// Set variables
xy = 26;
rot = -6.016; 
start = 1;
end = 798;
reverse = 1;
directory = "/Users/yuanz/Desktop/xy25-36/" //working directory


//Phase
open(directory+"data/xy"+xy+"c1.tif"); //open file
selectWindow("xy"+xy+"c1.tif");
run("Make Substack...", " slices="+start+"-"+end); //make substack of desired frame range
run("Rotate... ", "angle=rot grid=1 interpolation=Bilinear stack"); // rotate substack
if (reverse == 1) {
	run("Flip Vertically", "stack"); //if needed, flip the substack vertically
}
saveAs("Tiff", directory+"/xy"+xy+"/xy"+xy+"c1_t.tif"); //save stack to xy_pos directory
//export sequence of individual images to xy_pos/channel_type/raw for matlab code processing
run("Image Sequence... ", "format=TIFF save=["+directory+"/xy"+xy+"/phase/raw/]"); 

//Create threshold from phase
//run("Threshold...");
setOption("BlackBackground", true);
run("Convert to Mask", "method=Default background=Dark calculate");
run("Invert", "stack");
saveAs("Tiff", directory+"/xy"+xy+"/xy"+xy+"c1_pha_t.tif");
run("Image Sequence... ", "format=TIFF save=["+directory+"/xy"+xy+"/thresh/raw/]"); 


//Fluorescence
open(directory+"data/xy"+xy+"c2.tif");
selectWindow("xy"+xy+"c2.tif");
run("Make Substack...", " slices="+start+"-"+end);
run("Rotate... ", "angle=rot grid=1 interpolation=Bilinear stack");
if (reverse == 1) {
	run("Flip Vertically", "stack");
}
run("Subtract Background...", "rolling=50 stack"); //substract background
saveAs("Tiff", directory+"/xy"+xy+"/xy"+xy+"c2_t.tif");
run("Image Sequence... ", "format=TIFF save=["+directory+"/xy"+xy+"/flu/raw/]"); 


//Nuclear
open(directory+"data/xy"+xy+"c3.tif");
selectWindow("xy"+xy+"c3.tif");
run("Make Substack...", " slices="+start+"-"+end);
run("Rotate... ", "angle=rot grid=1 interpolation=Bilinear stack");
if (reverse == 1) {
	run("Flip Vertically", "stack");
}
run("Subtract Background...", "rolling=50 stack");
saveAs("Tiff", directory+"/xy"+xy+"/xy"+xy+"c3_t.tif");
run("Image Sequence... ", "format=TIFF save=["+directory+"/xy"+xy+"/nuc/raw/]"); 


//Close all windows
run("Close All");