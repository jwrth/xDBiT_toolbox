// macro to stitch images from Configfiles

#@ File (label="Select directory with folders containing the image sequences: ", style="directory") DirPath
#@ File (label="Base path of folders with configuration files: ", style="directory", value="Default") ConfigBasePath
#@ String (label="Unique pattern in ConfigFiles: ", value="projection") unique_configpattern
#@ Boolean (label="Use registered Configurations file? ") registered
#@ String (label="Input file format: ", choices={".tif", ".png", ".jpg"}, style="listBox") file_ending
#@ String (label="Patterns in input files (Comma-separated, e.g. 'c1,c2,c3': ") patterns_file
#@ String (label="Pattern in target directories (Optional): ", value="") pattern_dir
#@ Integer (label="First file index i: ") index_start
//#@ String (label="Special ID prefix (E.g. 'm'): ", value="m") id_prefix
//#@ String (label="Special ID suffix (E.g. '_ORG'): ", value="_ORG") id_suffix
#@ String (label="Flip? ", choices={"Vertically", "Horizontally", "Do not flip"}, style="listBox") flip
#@ Boolean (label="Overwrite output file if it exists? ") overwrite

setBatchMode(true); // run all processes in background

macro  ApplyStitchingToMultipleFolders
{
	sep = File.separator;
	DirPath = DirPath + sep;
	DirList = getFileList(DirPath);
	DirList = GetDirectories(DirList, DirPath);
	DirList = BackslashOnly(DirList, DirPath); // substitue all forward slashes with backslashes

	BasePath = NavigateUp(DirPath, 1, sep);

	// check config file base path
	if (ConfigBasePath == "Default"){
		ConfigBasePath = BasePath + "Processed" + sep;
	}
	else {
		ConfigBasePath = ConfigBasePath + sep;
	}
	ConfigFolders = getFileList(ConfigBasePath);
	ConfigFolders = GetDirectories(ConfigFolders, ConfigBasePath);
	ConfigFolders = BackslashOnly(ConfigFolders, ConfigBasePath); // substitue all forward slashes with backslashes

	// extract patterns
	patterns_file_ar = split(patterns_file, ",");
	
	// check for pattern_dir
	if(pattern_dir.length>0){
		DirList = GetFilesWithStr(DirList, pattern_dir);
		ConfigFolders = GetFilesWithStr(ConfigFolders, pattern_dir);

	}
	
	for(Cpt=0; Cpt<DirList.length; Cpt++){
		// select one of the folders
		Folder = DirList[Cpt];
		FolderName = split(Folder, sep);
		FolderName = FolderName[0];
		FolderPath = DirPath + Folder;
		SaveDir = DirPath + "Stitched" + sep;

		ConfigFolder = ConfigFolders[Cpt];
		ConfigPath = ConfigBasePath + ConfigFolder;

		if(!File.exists(SaveDir)){
			File.makeDirectory(SaveDir);
		}

		for(i=0; i<patterns_file_ar.length; i++){
			pattern_file = patterns_file_ar[i];
			SaveName = FolderName + "_" + pattern_file + "_stitched.tif";
			SavePath = SaveDir + SaveName;

			// get image sequence files in directory
			AllFiles = getFileList(FolderPath);
			ImgSeqFiles = GetFilesWithExt(AllFiles, file_ending);
			ImgSeqFiles = GetFilesWithStr(ImgSeqFiles, pattern_file);
			
			// get list of all possible configuration files
			AllConfigFiles = getFileList(ConfigPath);
			InputConfigFiles = GetFilesWithStr(AllConfigFiles, "TileConfiguration");

			if ( registered ){
				InputConfigFiles = GetFilesWithExt(InputConfigFiles, ".registered.txt");
			}
			else {
				InputConfigFiles = GetFilesWithoutExt(InputConfigFiles, ".registered.txt");
			}
			
			InputConfigFile = GetFilesWithStr(InputConfigFiles, unique_configpattern);

			//InputConfigPath = ConfigBasePath + Folder + "TileConfiguration_" + FolderName + "_projection.registered.txt";
			if(InputConfigFile.length == 1){
				InputConfigFile = InputConfigFile[0];
				print(InputConfigFile);
				InputConfigPath = ConfigBasePath + Folder + InputConfigFile;
			}
			else{
				if(InputConfigFile.length < 1){
					exit("No configuration files found.");
				}
				else{
					exit("More than one possible configuration file found.");
				}
			}

			ConfigFile = "TileConfiguration_" + FolderName + "_" + pattern_file + ".registered.txt";
			OutputConfigPath = DirPath + Folder + ConfigFile;
			
			ChangeConfigFile(InputConfigPath, OutputConfigPath, ImgSeqFiles);

			StitchImagesByFile(FolderPath, SavePath, ConfigFile);
		}
	}
	print("Finished stitching in all directories.");
}

function ChangeConfigFile(InputConfigPath, OutputConfigPath, FileList){
	outfile=File.open(OutputConfigPath);
	print(InputConfigPath);
	filestring=File.openAsString(InputConfigPath);
	rows=split(filestring, "\n");

	// get number of entries
	n_entries = GetFilesWithStr(rows, file_ending);
	n_entries = n_entries.length;

	// sort file list
	FileList = Array.sort(FileList);

	// change file names in config files using list of files
	count=0;
	if (FileList.length == n_entries){
		for (i=0; i<rows.length; i++){
			if (count>3){
				file_id = i - 4;
				entries = split(rows[i],";"); // split line in parts
				filename = entries[0];
				
				// create output file name
				entries[0] = FileList[file_id]; // substitute old file name with new file name
				//entries[0] = new_name;
				newline = String.join(entries, ";"); // generate new line
				print(outfile, newline);
			}
			else{
				print(outfile, rows[i]);
			}
			count=count+1;
		}
	}
	else{
		exit("Error during modification of Tile configuration files: Length of file list [" + FileList.length + "] and number of entries [" + n_entries + "] not the same.");
	}
	File.close(outfile);
}

function StitchImagesByFile(FolderPath, SavePath, ConfigFile)
{
	//Do stitching
	//print(FolderPath);
	//print(SavePath);
	//print(ConfigFile);
	//FolderPath = replace(FolderPath, "\\", "/");
	//FolderPath = "N:/01 HPC/03 Team Meier/08_Projects/37_Spatial_Barcoding/37_30/images/37_30_10X_13x10_20%/B=0/S=0";
	//print(FolderPath);
	print("Start stitching...");
	run("Grid/Collection stitching", "type=[Positions from file] " + 
		"order=[Defined by TileConfiguration] " + 
		"directory=[" + FolderPath + "] " + 
		"layout_file=" + ConfigFile + " " + 
		"fusion_method=[Linear Blending] " + 
		"regression_threshold=0.30 " + 
		"max/avg_displacement_threshold=2.50 " + 
		"absolute_displacement_threshold=3.50 " + 
		"computation_parameters=[Save memory (but be slower)] " + 
		"image_output=[Fuse and display]");


	if (flip != "Do not flip"){
		run("Flip " + flip);
		print("Image flipped");
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

// Functions
function NavigateUp(directory, n, sep){
	dirsplit = split(directory, sep);
	outdir = Array.trim(dirsplit, dirsplit.length-n);
	outdir = String.join(outdir, sep) + sep;

	return outdir;
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

function GetFilesWithoutExt(FileList,Extension)
{
	NewFileList = newArray;
	for (Cpt = 0; Cpt < FileList.length; Cpt++) {
		if(endsWith(FileList[Cpt], Extension) == false){
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

function GetFilesWithoutStr(FileList,Pattern)
{
	NewFileList = newArray;
	for (Cpt = 0; Cpt < FileList.length; Cpt++) {
		Idx = lastIndexOf(FileList[Cpt], Pattern);
		if (Idx == -1) {
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

