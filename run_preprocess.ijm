/*
 * ImageJ Macro for automatic image stack processing
 * Yuan Zhao 05/04/2015
 *
 * Manual steps:
 * 0) Place input image stacks in directory/data/ folder
 * 1) MANUALLY set the following variables as rows in config.txt
 *    a) xy; positions to process (ie 01,02,03,...)
 *    b) rot; angle to be rotated for each xy, measured from phase channel where + is CW and - is CCW (ie -0.74,-0.6,...)
 *    c) start, end; range of 'good' frames to keep; must be consistent across all positions (ie 1,796)
 *    d) reverse; 1 if stack needs to be vertically flipped so traps face "up" (ie 1,1,0,1,...)
 *    e) thrsh1, thrsh2; lower and upper threshold limits (ie 500,1500)
 *    f) directory; working directory path for folder root (ie /Volumes/Data/experiment)
 *    g) phase_ch_prefix; phase channel prefix (ie c1)
 *    h) flu_ch_prefix; all fluorescent channel prefixes (ie c2,c3,...)
 *
 * Fluorescent channel macro:
 * 1) Duplicate raw stack, make substack based on start/end frames
 * 2) Rotate image stack, flip vertically if necessary
 * 3) Substract background
 * 4) Save fluorescent stack as .tif, export frames individually
 *
 * Phase channel / threshold macro:
 * 1) Duplicate raw stack, make substack based on start/end frames
 * 2) Rotate image stack, flip vertically if necessary
 * 3) Save phase stack as .tif, export frames individually
 * 4) Threshold based on thrsh1, thrsh2 parameters; apply mask
 * 5) Analyze particles to find bounding rectangle(s) of chip/trap features; use the upper bound & height of
 *    bounding rectangle to clear the area above/below, and fill the holes in rectangles above/below traps
 * 6) Save threshold stack as .tif, export frames individually
 *
 */


//Import manually set variable configuration
config = File.openAsString("/Users/yuanz/Git/aging-chip/run_preprocess_config.txt");

config = split(config,"\n"); //split each line of config
xyarray = split(config[0],","); //split string of xy positions into xyarray
rotarray = split(config[1],","); //split string of rot degrees into rotarray
startend = split(config[2],",");
start = startend[0]; end = startend[1]; //define start and end frames
reversearray = split(config[3],","); //split string of reversal flags into reversearray
thrshall = split(config[4],",");
thrsh1 = thrshall[0]; thrsh2 = thrshall[1]; //define lower and upper thresolds
directory = config[5]; //define working directory
phase_ch_prefix = config[6]; //define phase channel prefix
flu_ch_prefix = split(config[7],",");


//Function for processing and subtracting background from fluorescent channels
//cname prefix (ie c1, c1_thr, c2...)
function processFlu(xy,rot,start,end,reverse,directory,cname) {
	open(directory+"data/xy"+xy+cname+".tif");
	//run("Duplicate...", "duplicate");
	run("Make Substack...", " slices="+start+"-"+end);
	run("Rotate... ", "angle=rot grid=1 interpolation=Bilinear stack");
	if (reverse == 1) {
		run("Flip Vertically", "stack");
	}
	run("Subtract Background...", "rolling=100 stack"); //substract background
	saveAs("Tiff", directory+"/xy"+xy+"/xy"+xy+"_"+cname+"_t.tif"); //save stack to xy_pos directory
	run("Image Sequence... ", "format=TIFF start=1 save=["+directory+"/xy"+xy+"/"+cname+"/]"); //export sequence of individual images to xy_pos/flu_prefix
	run("Close All"); //Close all windows
}


//Function for processing phase images and creating/cleaning threshold files
function processPhTh(xy,rot,start,end,reverse,thrsh1,thrsh2,directory,cname) {
	open(directory+"data/xy"+xy+cname+".tif"); //open file
	//run("Duplicate...", "duplicate");
	run("Make Substack...", " slices="+start+"-"+end); //make substack of desired frame range
	run("Rotate... ", "angle=rot grid=1 interpolation=Bilinear stack"); // rotate substack
	if (reverse == 1) {
		run("Flip Vertically", "stack"); //if needed, flip the substack vertically
	}
	saveAs("Tiff", directory+"/xy"+xy+"/xy"+xy+"_"+cname+"_t.tif");
	run("Image Sequence... ", "format=TIFF start=1 save=["+directory+"/xy"+xy+"/"+cname+"/]");

	//Create threshold from phase
	setAutoThreshold("Default");
	setThreshold(thrsh1, thrsh2);
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=Default background=Light");
	//Measure perimeter bounding to automatically clear debris artifacts / clogged cells, based on first frame and applied to entire stack
	setAutoThreshold("Default dark"); //threshold existing mask at 100% to find trap features
	run("Set Measurements...", "bounding redirect=None decimal=3"); //set perimeter bounding measurement to yes
	run("Analyze Particles...", "size=5000-Infinity display clear include"); //only consider particles area>5000 to exclude artifacts
	topheight = 18; //set the height of the top trap feature to 18px
	bottomheight = 32; //set the height of the bottom trap feature to 32px
	if (nResults == 1) { //if there is 1 particle detected
		BY = getResult("BY",0); //upper Y of bounding rectangle
		height = getResult("Height",0); //height of bounding rectangle
		//clearing black space
		makeRectangle(0,0,512,BY); //(start x, start y, width going right, height going down)
		run("Clear", "stack"); //clear the upper black space of entire stack
		makeRectangle(0,BY+height,512,512);
		run("Clear", "stack"); //clear the lower black space
		//filling white space
		makeRectangle(0,BY,512,topheight);
		run("Fill", "stack");
		makeRectangle(0,BY+height-bottomheight,512,bottomheight);
		run("Fill", "stack");
	}
	if (nResults == 2) { //if there are 2 particles
		BY1 = getResult("BY",0); //upper Y of first particle's bounding rectangle
		BY2 = getResult("BY",1); //upper Y of second particle's bound rect
		height = getResult("Height",1); //height of the second particle's bound rect
		//clearing black space
		makeRectangle(0,0,512,BY1); //(start x, start y, width going right, height going down)
		run("Clear", "stack"); //clear the upper black space of entire stack
		makeRectangle(0,BY2+height,512,512);
		run("Clear", "stack"); //clear the lower black space
		//filling white space
		makeRectangle(0,BY1,512,topheight);
		run("Fill", "stack");
		makeRectangle(0,BY2+height-bottomheight,512,bottomheight);
		run("Fill", "stack");
	}
	saveAs("Tiff", directory+"/xy"+xy+"/xy"+xy+"_"+cname+"_thr_t.tif");
	run("Image Sequence... ", "format=TIFF start=1 save=["+directory+"/xy"+xy+"/"+cname+"_thr/]");
	selectWindow("Results");
	run("Close All"); //Close all windows
}



///////////////////
// Initialize empty xy_pos/channel/raw folders
ph_thrsh_prefix = newArray(phase_ch_prefix,phase_ch_prefix+"_thr"); // add _thr to phase channel prefix for threshold prefix
folders = Array.concat(ph_thrsh_prefix,flu_ch_prefix); //concat ph, thrsh prefixes with flu channel prefixes; use this array for folder names
//Array.print(folders);
for (i = 0; i<xyarray.length; i++) { //for every position in 1D xyarray
	xy = xyarray[i];
	myDir1 = directory+"xy"+xy; //make xy_pos folders
	File.makeDirectory(myDir1);
	for (j = 0; j<folders.length; j++) { //for each folder listed in folders
		subf = folders[j];
		//make a subfolder
		myDir2 = directory+"xy"+xy+"/"+subf;
		File.makeDirectory(myDir2);
	}
}

// For every position in xy array, run the PhTh and Flu codes for all the channels
for (i = 0; i<xyarray.length; i++) { //for every position in 1D xyarray
	xy = xyarray[i];
	rot = parseFloat(rotarray[i]);
	reverse = parseInt(reversearray[i]);

	// Phase channel; exports processed phase and threshold
	processPhTh(xy,rot,start,end,reverse,thrsh1,thrsh2,directory,phase_ch_prefix);

	// Fluorescent channels
	for (j = 0; j<flu_ch_prefix.length; j++) { //for every flu channel in flu_ch_prefix array
		curr_flu = flu_ch_prefix[j];
		processFlu(xy,rot,start,end,reverse,directory,curr_flu);
	}
	run("Close All"); //Close all windows
}
///////////////////
