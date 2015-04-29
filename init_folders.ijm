directory = "/Volumes/Data HD/Workspace/xy25-36/"

channels = newArray("phase","thresh","flu","nuc");
positions = newArray("26","27","28");


for (i = 0; i<positions.length; i++) { //for every position in 1D positions array
	xy = positions[i]; //print(xy);
	myDir1 = directory+"xy"+xy; //make xy_pos folders
	File.makeDirectory(myDir1);
	for (j = 0; j<channels.length; j++) { 
		subf = channels[j]; //print(subf); 
		//make a subfolder for pha, thrsh, flu, and nuc
		myDir2 = directory+"xy"+xy+"/"+subf;
		File.makeDirectory(myDir2);
		//make a /raw/ subfolder in each of the four channels
		myDir3 = directory+"xy"+xy+"/"+subf+"/raw"; print(myDir3); 
		File.makeDirectory(myDir3);
		//if (!File.exists(myDir))
		//	exit("Unable to create directory");
		//print("");
		//print(myDir);
	}
}