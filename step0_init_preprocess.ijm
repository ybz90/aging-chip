/* c1 = phase; c2 = fluor; c3 = nucl marker; c1_pha = threshold
 *  
 * Manual steps: 
 * 0) Place raw image stacks in data folder; preprocess.jim in folder root
 * 1) MANUALLY set the following variables in 1D array
 *    a) xy = positions to process
 *    b) rot = angle to be rotated, measured in phase (+ CW, - CCW)
 *    c) start, end = range of 'good' frames to keep; must be consistient across all positions, so 
 *       choose the lowest of all the positions
 *    d) reverse = 1 if stack needs to be vertically flipped so straps face "up", = 0 if no flip needed
 *    e) directory = working directory path for folder root
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



//////////
// Set variables
xyarray = newArray("27");
rotarray = newArray("-6.02"); 
startarray = newArray("1");
endarray = newArray("796");
reversearray = newArray("1");

directory = "/Volumes/Data HD/Workspace/xy25-36/" //working directory


//////////
// Initialize empty xy_pos/channel/raw folders
for (i = 0; i<xyarray.length; i++) { //for every position in 1D xyarray
	xy = xyarray[i]; //print(xy);
	myDir1 = directory+"xy"+xy; //make xy_pos folders
	File.makeDirectory(myDir1);
	channels = newArray("phase","thresh","flu","nuc");
	for (j = 0; j<channels.length; j++) { 
		subf = channels[j]; //print(subf); 
		//make a subfolder for pha, thrsh, flu, and nuc
		myDir2 = directory+"xy"+xy+"/"+subf;
		File.makeDirectory(myDir2);
		//make a /raw/ subfolder in each of the four channels
		myDir3 = directory+"xy"+xy+"/"+subf+"/raw"; //print(myDir3); 
		File.makeDirectory(myDir3);
		//if (!File.exists(myDir))
		//	exit("Unable to create directory");
		//print("");
		//print(myDir);
	}
}


//////////
function processStacks(xy,rot,start,end,reverse,directory) {
	//Process images for each of the four channel/raw subfolders, and export as stacks as well
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
}


// Run processStacks(xy,rot,start,end,reverse,directory) for all positions in xyarray
for (i = 0; i<xyarray.length; i++) { //for every position in 1D xyarray
	xy = xyarray[i]; 
	rot = parseFloat(rotarray[i]);
	start = parseInt(startarray[i]);
	end = parseInt(endarray[i]);
	reverse = parseInt(reversearray[i]);
	processStacks(xy,rot,start,end,reverse,directory);
}