// Script to flat-field correct a channel in multiple folders with tile scans.

#@ File (label="Select directory with folders containing the image sequences", style="directory") DirPath
#@ String (label="Pattern in target directories (Optional): ", value="") pattern
#@ String (label="Channels to process (comma-separated; e.g. 'C00,C01')") channels
#@ Integer (label="First file index i: ") index_start
#@ Integer (label="Index length (e.g. 4 for '0000'): ") index_length
//#@ String (label="Histogram settings for each channel (In this syntax: 'Min1,Max1;None;auto'") minmax
#@ String (label="Convert to following bit type before correction:", choices={"8-bit","16-bit","RGB Color"}, style="listBox") bit_type
#@ String (label="Blurring radius for pseudo flatfield correction. '0' if no correction. Comma-separated for each channel (e.g. 50,0,50): ") radius

setBatchMode(true); // run all processes in background

macro  ApplyFunctionToMultipleFolders
{

	// get directories
	sep = File.separator;
	DirPath = DirPath + sep;
	DirList = getFileList(DirPath);
	DirList = GetDirectories(DirList, DirPath);
	DirList = BackslashOnly(DirList, DirPath); // substitue all forward slashes with backslashes

	BasePath = NavigateUp(DirPath, 1, sep);

	// check for pattern
	if(pattern.length>0){
		DirList = GetFilesWithStr(DirList, pattern);
	}

	channel_ar = split(channels, ",");
	radius_ar = split(radius, ",");

	// create output directory
	//OutDir = BasePath + "Processed" + sep;

	// create save directory in case it does not exist yet
	//if(!File.exists(OutDir)){
	//	File.makeDirectory(OutDir);
	//	}

	// save settings in file
	outfile=File.open(DirPath + "settings.txt");
	print(outfile, "Settings for Image Projection in Directory " + DirPath); // title
	print(outfile, "Pattern in input folders: " + pattern);
	print(outfile, "Channels: " + channels);
	//print(outfile, "Pseudo flatfield correction yes/no: " + correction);
	print(outfile, "Pseudo flatfield correction radius: " + radius);
	print(outfile, "Output format: " + bit_type);
	//print(outfile, "Projection method: " + projection_method);
	
	for(Cpt=0; Cpt<DirList.length; Cpt++){
		
		Folder = DirList[Cpt];
		FolderName = split(Folder, sep);
		FolderName = FolderName[0];

	
		// create save directory
		//SaveDir = OutDir + Folder;

		// create save directory in case it does not exist yet
		//if(!File.exists(SaveDir)){
		//	File.makeDirectory(SaveDir);
		//}

		//SaveName = RootName + "_C0" + channel - 1 + "_stitched.tif";

		FolderPath = DirPath + Folder;

		// here the save directory is the input directory
		SaveDir = FolderPath;
		//SavePath = NameExist(SaveDir + SaveName, sep=sep);

		FlatfieldCorrection();
		
	}
	print("Finished processing of all directories.");
}

function FlatfieldCorrection()
{

	print("Load images...");
	for(i=0; i<channel_ar.length; i++){
		// load image set 1
		run("Image Sequence...", "open=[" + FolderPath + "] file=" + channel_ar[i] + " sort");
		
		// rename image using channel name
		
		rename(channel_ar[i]);
		
		// change bit type as specified
		run(bit_type);

		if( radius_ar[i] > 0 ){
			// run pseudo flat field correction
			print("Run pseudo flat field correction for channel " + channel_ar[i]);
			run("Pseudo flat field correction", "blurring=" + radius_ar[i] + " hide stack");
		}
	
		// save files
		print("Save to " + SaveDir);
		SaveName = FolderName + "_" + channel_ar[i] + "ffcm";
	
		SaveDirList = getFileList(SaveDir);
		SaveDirList = GetFilesWithStr(SaveDirList,SaveName);
	
		if (SaveDirList.length > 0){
			print("Save directory contains files with output name. Delete files...");
			for (i=0; i<SaveDirList.length; i++){
			file = SaveDirList[i];
			delete_return = File.delete(SaveDir + file);
			}
			print("Files deleted.");
		}
		
		print(SaveName);
		print("Save image sequence...");
		// run("Image Sequence... ", "format=TIFF start=" + index_start + " digits=" + index_length + " name=" + SaveName + " save=[" + SaveDir + "]"); // olde ImageJ version
		run("Image Sequence... ", "format=TIFF start=" + index_start + " digits=" + index_length + " name=[" + SaveName + "] dir=[" + SaveDir + "]");

		close();
	}

}


// Functions


function NavigateUp(directory, n, sep){
	dirsplit = split(directory, sep);
	outdir = Array.trim(dirsplit, dirsplit.length-n);
	outdir = String.join(outdir, sep) + sep;

	return outdir;
}

function getMinMaxStack(){
	stack_mins = newArray();
	stack_maxs = newArray();
	for (n=1; n<=nSlices; n++){
		setSlice(n);
		getMinAndMax(min, max);
		stack_mins = Array.concat(stack_mins,min);
		stack_maxs = Array.concat(stack_maxs,max);
	}
	Array.getStatistics(stack_mins, stack_min, max, mean, stdDev);
	Array.getStatistics(stack_maxs, min, stack_max, mean, stdDev);
	stack_minmax = newArray(stack_min, stack_max);
	return stack_minmax;
}

function GetRootName(FileName,Seperator)
{
	Idx = GetStringIdx(FileName,Seperator);
	RootName = substring(FileName,0,Idx[Idx.length-1]);
	return RootName;
}

function GetStringIdx(Str,Pattern)
{
	Idx = -1;
	Idxs = newArray;
	do{
		Idx = indexOf(Str,Pattern,Idx+1);
		if(Idx!=-1){
			Idxs = Array.concat(Idxs,Idx);
		}
	} while(Idx!=-1)

	return Idxs;
}

function NameExist(Path, sep="\\") {
	// Checks whether the file exists and if yes it adds a running number.
	
	if( endsWith(Path, sep) ){
		RootName = substring(Path, 0, lengthOf(Path)-1);
		Extension = sep;
	}
	else {
		IdxP = GetStringIdx(Path,".");
		RootName = substring(Path, 0, IdxP[IdxP.length-1]);
		Extension = substring(Path, IdxP[IdxP.length-1],lengthOf(Path));
	}

	NewPath = Path;
	Cpt = 0;
	while( File.exists(NewPath) ){
		Cpt = Cpt+1;
		NewPath = RootName + "(" +  d2s(Cpt,0) + ")" + Extension;
	}

	return NewPath;
}

function GetFilesWithExt(FileList,Extension)
{
	NewFileList = newArray;
	for (Cpt = 0; Cpt < FileList.length; Cpt++) {
		if(endsWith(FileList[Cpt], Extension)){
			NewFileList = Array.concat(NewFileList,FileList[Cpt]);
		}
	}
	
	return NewFileList;
}

function CloseAllImagesExpect(KeepID)
{

IDs = newArray(nImages);
for(Cpt=1; Cpt<nImages+1; Cpt++){
	selectImage(Cpt);
	IDs[Cpt-1] = getImageID();
}
for(Cpt=0; Cpt<IDs.length; Cpt++){
	if(IDs[Cpt]!=KeepID){
		selectImage(IDs[Cpt]);
		close();
	}
}

}

function GetFilesWithStr(FileList,Pattern)
{
	NewFileList = newArray;
	for (Cpt = 0; Cpt < FileList.length; Cpt++) {
		Idx = lastIndexOf(FileList[Cpt], Pattern);
		if (Idx !=-1) {
			NewFileList = Array.concat(NewFileList,FileList[Cpt]);
		}
	}
	
	return NewFileList;
}

function GetDirectories(FileList, DirPath)
{
	NewFileList = newArray;
	for (Cpt = 0; Cpt < FileList.length; Cpt++) {
		if (File.isDirectory(DirPath + FileList[Cpt])) {
			currentdir = FileList[Cpt];
			NewFileList = Array.concat(NewFileList,currentdir);
		}
	}
	
	return NewFileList;
}

function BackslashOnly(FileList, DirPath)
{
	NewFileList = newArray;
	for (Cpt = 0; Cpt < FileList.length; Cpt++) {
		currentdir = FileList[Cpt];
		currentdir = replace(currentdir, "/", "\\");
		NewFileList = Array.concat(NewFileList,currentdir);
	}
	
	return NewFileList;
}

