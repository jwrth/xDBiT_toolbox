setBatchMode(true); // run all processes in background

macro  Generic_BatchProcess_Macro
{
	DirName = getDirectory("Choose Directory with files");
	FileList = getFileList(DirName);
	FileList = GetFilesWithExt(FileList,".tif");

	for(Cpt=0; Cpt<FileList.length; Cpt++){
		
		FileName = FileList[Cpt];
		RootName = GetRootName(FileName,".");
		
		SaveDir = "Flipped\\";
		SaveDir = DirName + SaveDir;
		if(!File.exists(SaveDir)){
			File.makeDirectory(SaveDir);
		}

		SaveName = RootName + "_flipped.tif";

		FilePath = DirName + FileName;
		SavePath = NameExist(SaveDir + SaveName);

		TreatSingleFile(FilePath,SavePath);
	}
}

function TreatSingleFile(FilePath,SavePath)
{

	CloseAllImagesExpect(0);
	open(FilePath);
	OrigID = getImageID();

	//Flip image
	print("Processing " + FilePath); 
	run("Flip Vertically");
	
	saveAs("Tiff", SavePath);
	print("Saved to " + SavePath);
	close();
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
