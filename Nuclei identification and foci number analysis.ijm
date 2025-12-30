
dir1 = getDirectory("Choose Source Directory ");
list = getFileList(dir1)
len = lengthOf(dir1);
path = substring(dir1, 0, len-1);
last_pos = indexOf(path, "\\");
second_last_pos = last_pos;
while (indexOf(path, "\\", second_last_pos+1) != -1) {
    second_last_pos = indexOf(path, "\\", second_last_pos+1);
}

// Load Bio-Formats library
run("Bio-Formats Macro Extensions");

dapiTh=8000;
GreenTh=4000;
RedTh=4000;
mfz=1;
mns=500;
viewD=5;
viewR=10;

		

Dialog.create("System control");
Dialog.addMessage("What you want to analyze : ");
Dialog.addChoice("Analyze:", newArray("Foci numbers","Signal intensity"));
Dialog.show();
AnalyzeChoice  = Dialog.getChoice();


if( AnalyzeChoice=="Foci numbers"){
	print("'----------------Run informaiton of Foci analysis");	
	getDateAndTime(year, month, week, day, hour, min, sec, msec);
	print("Date: "+day+"/"+month+1+"/"+year+"  Time: " +hour+":"+min+":"+sec);
		
	Dialog.create("System control");
	Dialog.addMessage("Analysis pararmeters : ");
	
	var min = 0;
	var max = 65535;
	Dialog.addSlider("Blue Threshold(x-16bit)", min, max, dapiTh);
	Dialog.addSlider("Red Threshold(x-16bit)", min, max, RedTh);
	Dialog.addNumber("Minimum foci size in pixels", mfz);
	Dialog.addNumber("Minimum nucleus size in pixels", mns);
	Dialog.addCheckbox("Slow analysis?", false);
	Dialog.addCheckbox("Save ROI?", false);
	Dialog.addCheckbox("Use ROI?", false);
	Dialog.addNumber(" Max 16bit for veiwer (Redx10^3):", viewR);
	Dialog.addNumber(" Max 16bit for veiwer (Dapix10^3):", viewD);
	//Dialog.addCheckbox("Auto DapiTH", false);
	Dialog.addChoice("Auto DapiTH:", newArray("None","MinimumTh","IsoTh"));
	Dialog.show();


	inNumber1 = Dialog.getNumber(); // Get numbers
	inNumber2 = Dialog.getNumber();
	inNumber4 = Dialog.getNumber();
	inNumber5 = Dialog.getNumber();
	inChoice  = " (red)";
	Slow = Dialog.getCheckbox();
	ROIs=Dialog.getCheckbox();
	ROIu=Dialog.getCheckbox();
	viewR=Dialog.getNumber()*1000;
	viewD=Dialog.getNumber()*1000;
	//Auto=Dialog.getCheckbox();
	ThAuto  = Dialog.getChoice();
	if( ROIs==1){
		Slow =1;
	}
	 
	print("Blue Threshold(x-16bit)", inNumber1);
	print("Red Threshold(x-16bit)", inNumber2);
	print("Minimum foci size in pixels", inNumber4);
	print("Minimum nucleus size in pixels", inNumber5);
	if( Slow==1){
		print("slow analysis"," Views(10^3): DAPI ",viewD/1000, "Red ",viewR/1000);
				}else{
					print("Fast analysis"," Views(10^3): DAPI ",viewD/1000, "Red ",viewR/1000);
				}
	if( ROIu==1){
		print("Run on saved ROI");
				}else{
					if(ThAuto!= "None" || ROIs ==1){
						print("Save ROI analysis using AutoTh DAPI" + ThAuto);
					}else{
					print("New ROI analysis");
					}
				}
	print(dir1);
}

if( AnalyzeChoice=="signal intensity"){
	print("'----------------Run informaiton of intensity analysis");
	getDateAndTime(year, month, week, day, hour, min, sec, msec);
	print("Date: "+day+"/"+month+1+"/"+year+"  Time: " +hour+":"+min+":"+sec);
	
	Dialog.create("System control");
	Dialog.addMessage("Analysis pararmeters : ");
	 
	var min = 0;
	var max = 65535;
	Dialog.addSlider("Blue Threshold(x-16bit)", min, max, 8000);
	Dialog.addNumber("Minimum nucleus size in pixels", 3000);
	Dialog.addCheckbox("Slow analysis?", false);
	Dialog.addCheckbox("Save ROI?", false);
	Dialog.addCheckbox("Use ROI?", false);
	Dialog.addNumber(" Max 16bit for veiwer (Redx10^3):", viewR);
	Dialog.addNumber(" Max 16bit for veiwer (Dapix10^3):", viewD);
	Dialog.addCheckbox("Print Autothresuold values?", false);
	//Dialog.addCheckbox("Auto DapiTH", false);
	Dialog.addChoice("Auto DapiTH:", newArray("None","MinimumTh","IsoTh"));
	Dialog.show();

	
	inNumber1 = Dialog.getNumber(); // Sliders are number too
	inChoice  = " (red)";
	inNumber5 = Dialog.getNumber();
	Slow = Dialog.getCheckbox();
	ROIs=Dialog.getCheckbox();
	ROIu=Dialog.getCheckbox();
	viewR=Dialog.getNumber()*1000;
	viewD=Dialog.getNumber()*1000;
	Mdark=Dialog.getCheckbox();
	//Auto=Dialog.getCheckbox();
	ThAuto  = Dialog.getChoice();
	if( ROIs==1){
		Slow =1;
	}
 
	print("Blue Threshold(x-16bit)", inNumber1);
	if( Slow==1){
		print("slow analysis"," Views(10^3): DAPI ",viewD/1000, "Red ",viewR/1000);
				}else{
					print("Fast analysis"," Views(10^3): DAPI ",viewD/1000, "Red ",viewR/1000);
				}
	if( ROIu==1){
		print("Run on saved ROI");

				}
	print(dir1);
}

if(isOpen("Summary") || isOpen("Results")){
	waitForUser("Warnning!: Summary or Results tables are not closed yet!");
}
if( AnalyzeChoice=="Foci numbers"){
	selectWindow("Log");
	saveAs("Text", substring(path, 0, second_last_pos)+"\\"+"Runs.txt");
}


for (i = 0; i < list.length; i++) {
	
	//working bit:
	call("ij.ImagePlus.setDefault16bitRange", 16);
			for (t=0; t<(3); t++) {
				setKeyDown("none");
				showProgress(t,3);
				wait(500);
				showStatus("Opening Imagej, (shift) stops the macro");
				selectWindow("Log");
			}
			if (isKeyDown("shift")){
					print("Analysis breaked!");
					break;	
			}

			open(dir1+list[i]);
			//run("Viewer", "open=["+dir1+n1+".tif]");
			run("Split Channels");
			selectWindow("C1-"+list[i]);
			run("Z Project...", "projection=[Max Intensity]");
			Image1=list[i];
			close("C1-"+list[i]);
			close("C3-"+list[i]);
			selectWindow("MAX_C1-"+list[i]);
			rename(Image1 + " (blue)");
			setMinAndMax(0, viewD);
			run("Duplicate...", "duplicate");
			rename("Adjust Th here Original blue");	
			setMinAndMax(0, viewD);
			run("Duplicate...", "duplicate");
			rename("Blue For View only");				
        	setMinAndMax(0, viewD);
        	run("Apply LUT");

			selectWindow("C2-"+list[i]);
			run("Z Project...", "projection=[Max Intensity]");	
			close("C2-"+list[i]);				
			selectWindow("MAX_C2-"+list[i]);
			rename(Image1 + " (red)");
			setMinAndMax(0, viewR);
			run("Duplicate...", "duplicate");
			rename("Adjust Th here Original red");
			setMinAndMax(0, viewR);
			run("Duplicate...", "duplicate");
			rename("Red For View only");				
        	setMinAndMax(0, viewR);
        	run("Apply LUT");
 	
			selectWindow(Image1 + " (blue)"); //For DAPI
			wait(1000);
			if( Slow==1){
				wait(1000);
				run("Threshold...");
				
			}
			if(ThAuto!="None"){
				if(ThAuto=="MinimumTh"){
				setAutoThreshold("Minimum dark");
				getThreshold(lower, upper);
				}else{
				setAutoThreshold("IsoData dark");
				getThreshold(lower, upper);
				}
			}else{
			setThreshold(inNumber1,65535); // set TH
			}
			selectWindow(Image1 + " (blue)"); //For DAPI
			wait(1000);
			if( Slow==1 && ROIu!=1 ){
				wait(1000);
				waitForUser("adjustTh");
				}
				
			// Set DAPI 
			selectWindow(Image1 + " (blue)");
			run("Convert to Mask");
			run("Fill Holes");
			run("Watershed");
			//waitForUser("Continue");
			run("Analyze Particles...", "size="+inNumber5+"-Infinity pixel circularity=0.00-1.00 show=Outlines exclude add");
			len = lengthOf(dir1);
			path = substring(dir1, 0, len-1);
			last_pos = indexOf(path, "\\");
			second_last_pos = last_pos;
			while (indexOf(path, "\\", second_last_pos+1) != -1) {
			    second_last_pos = indexOf(path, "\\", second_last_pos+1);
			}
				if( ROIs==1){
					if( Slow==1){
					selectWindow(Image1 + " (blue)");
					roiManager("Show All");
					selectWindow("Blue For View only");
					roiManager("Show All");
					selectWindow("Blue For View only");
					if ( roiManager("count") > 1){
						roiManager("Select", 0);
					}
					if(ThAuto=="MinimumTh"){
						Dialog.create("System control");
						Dialog.addChoice("Good ROI:", newArray("Yes","TryEdit"));
						Dialog.show();
						SaveRoI  = Dialog.getChoice();
						if (SaveRoI == "TryEdit"){
							selectWindow("ROI Manager");
							run("Close");
							selectWindow("Adjust Th here Original blue");	
							run("Duplicate...", "duplicate");
							rename("TryEditDapiTh");
							selectWindow("TryEditDapiTh");
							setMinAndMax(0, viewD);
							waitForUser("Delete high Dapi Regions & analyze again");
							selectWindow("TryEditDapiTh");
							setAutoThreshold("Minimum dark");
							selectWindow("TryEditDapiTh");
							run("Convert to Mask");
							run("Fill Holes");
							run("Watershed");
							//waitForUser("Continue");
							run("Analyze Particles...", "size="+inNumber5+"-Infinity pixel circularity=0.00-1.00 show=Outlines exclude add");
							selectWindow("Blue For View only");
							roiManager("Show All");
							selectWindow("Blue For View only");
							if ( roiManager("count") > 1){
								roiManager("Select", 0);						
								}
					}
					}
					waitForUser("Save ROI");	
					myfolder=substring(path, 0, second_last_pos)+"\\";
					// Create a new folder
					newFolderName = "ROI";
					File.makeDirectory(myfolder+"\\"+newFolderName);
						
					Saveat = substring(path, 0, second_last_pos)+"\\" +"ROI"+"\\"+substring(path,second_last_pos+1,len-1)+list[i]+".zip";	
					roiManager("Save", Saveat);
				}
				}
				if( ROIu==1){
					openat = substring(path, 0, second_last_pos)+"\\" +"ROI"+"\\"+substring(path,second_last_pos+1,len-1)+list[i]+".zip";
					run("ROI Manager...");
					selectWindow("ROI Manager");
					run("Close");
					roiManager("Open", openat);
					selectWindow(Image1 + " (blue)");
					roiManager("Show None");
					roiManager("Show All");
				}
			
				if( Slow==1 && ROIs == 0){
					selectWindow(Image1 + " (blue)");
					selectWindow("Blue For View only");
					roiManager("Show All");
					selectWindow("Adjust Th here Original blue");
					if(ThAuto!="None"){
						setThreshold(lower,65535);
						}else{
						setThreshold(inNumber1,65535); // set TH
					}
					roiManager("Show All");	
					waitForUser("ROI Continue");
				}
								
			if( AnalyzeChoice=="Foci numbers"){	
				selectWindow(Image1 + " (red)"); //For red
				wait(500);
				if( Slow==1){
					wait(1000);
					}
					
				setThreshold(inNumber2,65535); // set TH
				selectWindow(Image1 + " (red)"); //For red
				wait(500);
				if( Slow==1){
					wait(1000);
					}
				run("Convert to Mask");
				if( Slow==1 && ROIs == 0){
					if(inChoice ==" (red)"){
						selectWindow(Image1 + " (red)");
						selectWindow("Red For View only");
						roiManager("Show All");
						selectWindow("Adjust Th here Original red");
						setThreshold(inNumber2, 65535);
						roiManager("Show All");	
						waitForUser("Red Image Continue");
					}
				}
				
			}
			m = roiManager("count");
			// Count colocalized foci
			for (y=0; y<m ;y++) {
				selectWindow(Image1 + inChoice);
				run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
				roiManager("Select", y);
				run("Set Measurements...", "area integrated display redirect=None decimal=3");
				run("Measure");
				if( AnalyzeChoice=="Foci numbers"){
					roiManager("Select", y);
					run("Analyze Particles...", "size="+inNumber4+"-2000 pixel circularity=0.00-1.00 show=Nothing summarize");
					//waitForUser("Next nucleus");
				}	
			}
		
				
				
			if( AnalyzeChoice=="signal intensity"){
				selectWindow(Image1 + inChoice);
				roiManager("Show All");
				
				if( Mdark==1){
						selectWindow("Adjust Th here Original blue");
						setThreshold(inNumber1, 65535);
					if(ThAuto!="None"){
						if(ThAuto=="MinimumTh"){
						setAutoThreshold("Minimum dark");
						getThreshold(lower, upper);
						print(lower);
						}else{
						setAutoThreshold("IsoData dark");
						getThreshold(lower, upper);
						print(lower);
						}
					}
					}else{
					waitForUser("Take 4 Background...");
					}
			}		
			if( Slow==1 && ROIs == 0){
						waitForUser("Close?");
			}		
			run("Close All");
			selectWindow("ROI Manager");
			run("Close");
}
if( ROIs==1){
selectWindow("Log");
saveAs("Text", substring(path, 0, second_last_pos)+"\\" +"ROI"+"\\"+substring(path,second_last_pos+1,len-1)+"Log.txt");
}
waitForUser("The analysis compleated!");
