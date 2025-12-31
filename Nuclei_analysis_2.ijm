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
mns=500;
viewD=5000;
viewG=10;

Dialog.create("System control");
Dialog.addMessage("What you want to analyze : ");
Dialog.addChoice("Analyze:", newArray("Area and Circularity"));
Dialog.show();
AnalyzeChoice  = Dialog.getChoice();

if(AnalyzeChoice=="Area and Circularity"){
    print("'----------------Run information of Area and Circularity analysis");
    getDateAndTime(year, month, week, day, hour, min, sec, msec);
    print("Date: "+day+"/"+(month+1)+"/"+year+"  Time: " +hour+":"+min+":"+sec);

    Dialog.create("System control");
    Dialog.addMessage("Analysis parameters : ");

    var min = 0;
    var max = 65535;
    Dialog.addSlider("Blue Threshold(x-16bit)", min, max, dapiTh);
    Dialog.addSlider("Green Threshold(x-16bit)", min, max, GreenTh);
    Dialog.addNumber("Minimum nucleus size in pixels", mns);
    Dialog.addCheckbox("Slow analysis?", false);
    Dialog.addCheckbox("Save ROI?", false);
    Dialog.addCheckbox("Use ROI?", false);
    Dialog.addNumber("Max 16bit for viewer (Greenx10^3):", viewG);
    Dialog.addNumber("Max 16bit for viewer (Dapix10^3):", viewD);
    Dialog.addCheckbox("Print Autothreshold values?", false);
    Dialog.addChoice("Auto DapiTH:", newArray("None","MinimumTh","IsoTh"));
    Dialog.show();

    inNumber1 = Dialog.getNumber(); // Get numbers
    inNumber2 = Dialog.getNumber();
 //   inNumber4 = Dialog.getNumber();
    inNumber5 = Dialog.getNumber();
    Slow = Dialog.getCheckbox();
    ROIs = Dialog.getCheckbox();
    ROIu = Dialog.getCheckbox();
    viewG = Dialog.getNumber() * 1000;
  //  viewD = Dialog.getNumber() * 1000;
    ThAuto = Dialog.getChoice();

    print("Blue Threshold(x-16bit)", inNumber1);
    print("Green Threshold(x-16bit)", inNumber2);
 //   print("Minimum nucleus size in pixels", inNumber4);
    print("Minimum nucleus size in pixels", inNumber5);
    if( Slow==1){
        print("slow analysis"," Views(10^3): DAPI ",viewD/1000, "Green ",viewG/1000);
    }else{
        print("Fast analysis"," Views(10^3): DAPI ",viewD/1000, "Green ",viewG/1000);
    }
    if( ROIu==1){
        print("Run on saved ROI");
    }else{
        if(ThAuto != "None" || ROIs ==1){
            print("Save ROI analysis using AutoTh DAPI" + ThAuto);
        }else{
            print("New ROI analysis");
        }
    }
    print(dir1);
}

// Loop to process each file in the directory
for (i = 0; i < list.length; i++) {
    // Working bit:
    call("ij.ImagePlus.setDefault16bitRange", 16);
    for (t = 0; t < (3); t++) {
        setKeyDown("none");
        showProgress(t, 3);
        wait(500);
        showStatus("Opening Imagej, (shift) stops the macro");
        selectWindow("Log");
    }
    if (isKeyDown("shift")) {
        print("Analysis broken!");
        break;    
    }

    open(dir1 + list[i]); // Open the .oib image
    // Split channels for .oib (Bio-Formats will load this correctly)
    run("Split Channels");
    selectWindow("C1-" + list[i]);
//    run("Z Project...", "projection=[Max Intensity]");

    Image1 = list[i];
    close("C2-" + list[i]);
//    close("C3-" + list[i]);
//    selectWindow("MAX_C1-" + list[i]);
    rename(Image1 + " (blue)");
    setMinAndMax(0, viewD);
    run("Duplicate...", "duplicate");
    rename("Adjust Th here Original blue");    
    setMinAndMax(0, viewD);
    run("Duplicate...", "duplicate");
    rename("Blue For View only");                
    setMinAndMax(0, viewD);
    run("Apply LUT");

//    selectWindow("C2-" + list[i]);
//    run("Z Project...", "projection=[Max Intensity]");    
//    close("C2-" + list[i]);                
//    selectWindow("MAX_C2-" + list[i]);
//    rename(Image1 + " (green)");
//    setMinAndMax(0, viewG);
//    run("Duplicate...", "duplicate");
 //  rename("Adjust Th here Original green");
//    setMinAndMax(0, viewG);
//    run("Duplicate...", "duplicate");
//    rename("Green For View only");                
//    setMinAndMax(0, viewG);
//    run("Apply LUT");

    // Processing the DAPI (blue channel)
    selectWindow(Image1 + " (blue)"); // For DAPI
    wait(1000);
    if(Slow == 1) {
        wait(1000);
        run("Threshold...");
    }
    if (ThAuto != "None") {
        if (ThAuto == "MinimumTh") {
            setAutoThreshold("Minimum dark");
            getThreshold(lower, upper);
        } else {
            setAutoThreshold("IsoData dark");
            getThreshold(lower, upper);
        }
    } else {
        setThreshold(inNumber1, 65535); // Set Threshold
    }

    selectWindow(Image1 + " (blue)"); // For DAPI
    wait(1000);
    if(Slow == 1 && ROIu != 1) {
        wait(1000);
        waitForUser("Adjust Threshold");
    }

    // Set DAPI Mask
    selectWindow(Image1 + " (blue)");
    run("Convert to Mask");
    run("Fill Holes");
    run("Watershed");
    run("Analyze Particles...", "size=" + inNumber5 + "-Infinity pixel circularity=0.00-1.00 show=Outlines exclude add");

    // Process ROI if selected
    if (ROIs == 1) {
        if (Slow == 1) {
            selectWindow(Image1 + " (blue)");
            roiManager("Show All");
            selectWindow("Blue For View only");
            roiManager("Show All");
            selectWindow("Blue For View only");
            if (roiManager("count") > 1) {
                roiManager("Select", 0);
            }
            if (ThAuto == "MinimumTh") {
                Dialog.create("System control");
                Dialog.addChoice("Good ROI:", newArray("Yes", "TryEdit"));
                Dialog.show();
                SaveRoI = Dialog.getChoice();
                if (SaveRoI == "TryEdit") {
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
                    run("Analyze Particles...", "size=" + inNumber5 + "-Infinity pixel circularity=0.00-1.00 show=Outlines exclude add");
                    selectWindow("Blue For View only");
                    roiManager("Show All");
                    selectWindow("Blue For View only");
                    if (roiManager("count") > 1) {
                        roiManager("Select", 0);                        
                    }
                }
            }
        }
        waitForUser("Save ROI");
        myfolder = substring(path, 0, second_last_pos) + "\\";
        // Create a new folder
        newFolderName = "ROI";
        File.makeDirectory(myfolder + "\\" + newFolderName);

        Saveat = substring(path, 0, second_last_pos) + "\\" + "ROI" + "\\" + substring(path, second_last_pos + 1, len - 1) + list[i] + ".zip";    
        roiManager("Save", Saveat);
    }
    if (ROIu == 1) {
        openat = substring(path, 0, second_last_pos) + "\\" + "ROI" + "\\" + substring(path, second_last_pos + 1, len - 1) + list[i] + ".zip";
        run("ROI Manager...");
