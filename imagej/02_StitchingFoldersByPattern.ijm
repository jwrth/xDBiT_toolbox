#@ File (label="Select directory with folders containing the image sequences", style="directory") DirName
#@ Integer (label="Grid size x: ") GridX
#@ Integer (label="Grid size y: ") GridY
#@ Integer (label="Overlap (%): ", style="slider", min=0, max=100, value=20) overlap
#@ String (label="Type: ", choices={"Grid: row-by-row", "Grid: column-by-column", "Grid: snake by rows", "Grid: snake by columns", "Positions from file"}, style="listBox") type
#@ String (label="Order: ", choices={"Right & Down                ", "Left & Down", "Right & Up", "Left & Up"}, style="listBox") order
#@ Integer (label="First file index i: ") index_start
#@ Integer (label="Index length (e.g. 4 for '0000'): ") index_length
#@ String (label="Index prefix (e.g. M for 'M0000'): ") index_prefix
#@ String (label="Input file format: ", choices={".tif", ".png", ".jpg"}, style="listBox") file_ending
#@ String (label="Channels to process (comma-separated; e.g. 'C00,C01')") channels
#@ String (label="Pattern in target directories (Optional): ", value="") pattern_dir
#@ String (label="Convert to following bit type after stitching:", choices={"8-bit","16-bit","RGB Color"}, style="listBox") bit_type
#@ String (label="Flip? ", choices={"Vertically", "Horizontally", "Do not flip"}, style="listBox") flip
#@ Boolean (label="Compute overlap? ") compute
#@ Boolean (label="Overwrite output file if it exists? ") overwrite

setBatchMode(true); // run all processes in background

macro  ApplyStitchingToMultipleFolders
{
	sep = File.separator;
	DirName = DirName + sep;
	//DirName = getDirectory("Choose Directory with files");
	//prefix = getString("Prefix", "tile");
	//channel = getString("Channel to process", "1");
	DirList = getFileList(DirName);
	DirList = GetDirectories(DirList, DirName);
	DirList = BackslashOnly(DirList, DirName); // substitue all forward slashes with backslashes

	// check for pattern_files
	channels_ar = split(channels, ",");
	
	// check for pattern_dir
	if(pattern_dir.length>0){
		DirList = GetFilesWithStr(DirList, pattern_dir);
	}
	
	for(Cpt=0; Cpt<DirList.length; Cpt++){
		
		Folder = DirList[Cpt];
		FolderName = split(Folder, sep);
		FolderName = FolderName[0];
		FolderPath = DirName + Folder;
		SaveDir = DirName + "Stitched" + sep;

		if(!File.exists(SaveDir)){
			File.makeDirectory(SaveDir);
		}

		// determine structure of index in input files
		first_index = repeatStr("0", index_length - 1);
		first_index = index_prefix + first_index + toString(index_start);
		index_sub = "{" + repeatStr("i", index_length) + "}";
		


		for(i=0; i<channels_ar.length; i++){
			pattern_file = channels_ar[i];
			SaveName = FolderName + "_" + pattern_file + "_stitched.tif";
			SavePath = SaveDir + SaveName;

			// get image sequence files in directory
			ImgSeqFiles = getFileList(FolderPath);
			
			ImgSeqFiles = GetFilesWithExt(ImgSeqFiles, file_ending);
			ImgSeqFiles = GetFilesWithStr(ImgSeqFiles, pattern_file);
			
			// sort files
			ImgSeqFiles = Array.sort(ImgSeqFiles);
			FirstFile = ImgSeqFiles[0];
			
			idx = indexOf(FirstFile, first_index);
			prefix = FirstFile.substring(0, idx);
			suffix = FirstFile.substring(idx + first_index.length);
			FileName = prefix + index_prefix + index_sub + suffix;
			ConfigFile = "TileConfiguration_" + FolderName + "_" + pattern_file + ".txt";
			
			StitchImages(FolderPath, SavePath, FileName, ConfigFile);
		}
	}
	print("Finished stitching in all directories.");
}

function StitchImages(FolderPath, SavePath, FileName, ConfigFile)
{
	//Do stitching
	print("Start stitching...");
	if ( compute ){
		run("Grid/Collection stitching", "type=[" + type + "] " +
		"order=[" + order + "] " + 
		"grid_size_x=" + GridX + " " +
		"grid_size_y=" + GridY + " " + 
		"tile_overlap=" + overlap + " " + 
		"first_file_index_i=" + index_start + " " + 
		"directory=[" + FolderPath + "] " + 
		"file_names=[" + FileName + "] " + 
		"output_textfile_name=[" + ConfigFile + "] " + 
		"fusion_method=[Linear Blending] " + 
		"regression_threshold=0.30 " + 
		"max/avg_displacement_threshold=2.50 " + 
		"absolute_displacement_threshold=3.50 " + 
		"compute_overlap " +
		"computation_parameters=[Save memory (but be slower)] " + 
		"image_output=[Fuse and display]");
	}
	else {
		run("Grid/Collection stitching", "type=[" + type + "] " +
		"order=[" + order + "] " + 
		"grid_size_x=" + GridX + " " +
		"grid_size_y=" + GridY + " " + 
		"tile_overlap=" + overlap + " " + 
		"first_file_index_i=" + index_start + " " + 
		"directory=[" + FolderPath + "] " + 
		"file_names=[" + FileName + "] " + 
		"output_textfile_name=[" + ConfigFile + "] " + 
		"fusion_method=[Linear Blending] " + 
		"regression_threshold=0.30 " + 
		"max/avg_displacement_threshold=2.50 " + 
		"absolute_displacement_threshold=3.50 " + 
		"image_output=[Fuse and display]");
	}
	
	// change bit type as specified
	run(bit_type);
	print("Output bit type: " + bit_type);

	if (flip != "Do not flip"){
		run("Flip " + flip);
		print("Image flipped");
		RootName = GetRootName(SavePath,".");
		SavePath = RootName + "_flipped" + file_ending;
	}
	else{
		print("Image not flipped.");
	}

	if ( overwrite ){
		File.delete(SavePath);
		saveAs("Tiff", SavePath);
	}
	else {
		SavePath = NameExist(SavePath);
		saveAs("Tiff", SavePath);
	}

	close();

}

function repeatStr(string, repeats){
	newstring = "";
	for (i=0; i<repeats; i++){
		newstring = newstring + string;
	}
	return newstring;
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

function NameExist(Path) {
	// Checks whether the file exists and if yes it adds a running number.
	
	if( endsWith(Path, "\\") ){
		RootName = substring(Path, 0, lengthOf(Path)-1);
		Extension = "\\";
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

function GetDirectories(FileList, DirName)
{
	NewFileList = newArray;
	for (Cpt = 0; Cpt < FileList.length; Cpt++) {
		if (File.isDirectory(DirName + FileList[Cpt])) {
			currentdir = FileList[Cpt];
			NewFileList = Array.concat(NewFileList,currentdir);
		}
	}
	
	return NewFileList;
}

function BackslashOnly(FileList, DirName)
{
	NewFileList = newArray;
	for (Cpt = 0; Cpt < FileList.length; Cpt++) {
		currentdir = FileList[Cpt];
		currentdir = replace(currentdir, "/", "\\");
		NewFileList = Array.concat(NewFileList,currentdir);
	}
	
	return NewFileList;
}

