# Script collection for ams AG Process Relibaility team.
# Allow the creation of exportfiles based on
# data files provided by equipements.
# September 2015
# Emmanuel Chery
# Version 0.4




CreateExportFiles.Ace <- function(DegFileName,TCRFileName)
# Function called to create the exportfiles
# from a Qualitau Ace experiment.
# Use a Degradation and a TCR file to produce the exportfile
{
    # Device and Width parameters
    ListDevice <- read.delim("//fsup04/fntquap/Common/Qual/Process_Reliability/Process/amsReliability_R_Package/ListDeviceName.txt")

    DegFile <- read.delim(DegFileName,sep="")
    names(DegFile) <- c("Device","Failed","Lifetime[s]","StressTemp[K]","Positive Current[A]","Negative Current[A]","Duty Cycle","Pulse Width","Stress Type")
    TCRFile <- read.delim(TCRFileName,sep="",skip=1)
    names(TCRFile) <- c("Device","Rref[Ohm]","TCR[%/°C]","Rref[Ohm](High I)","TCR[%/°C](High I)","Rref[Ohm](Low I)","TCR[%/°C](Low I)","R[Ohm] (stress)","Temperature[°C](High I)","Temperature[°C](Low I)","Temperature[°C](Opt)","R[Ohm](Init Temp Low I)","R[Ohm](stress Temp Low I)","R[Ohm](stress Temp High I)")

    # Merge the files with additional data:
    # Split
    Split<-c(1,2)
    # Istress extracted from filename and mA is removed
    Istress <- strsplit(DegFileName,split="_")[[1]][3]
    Istress <- as.numeric(substr(Istress, 1, nchar(Istress)-2))
    # Temp extracted from filename and C is removed
    Temp <- strsplit(DegFileName,split="_")[[1]][4]
    Temp <- as.numeric(substr(Temp, 1, nchar(Temp)-1))
    # DeviceID, Length and width
    L <- c(200,400)
    DeviceID <- strsplit(DegFileName,split="_")[[1]][2]
    W <- ListDevice$Width[ListDevice$Device==DeviceID]

    # DataFrame creation
    NewFile <- data.frame(DegFile[,1:3],Split,Istress,L,W,Temp,DeviceID,TCRFile[,2:14])
    names(NewFile) <- c("Device","Failed","Lifetime[s]","Split","Istress","L","W","Temp","DeviceID","Rref[Ohm]","TCR[%/°C]","Rref[Ohm](High I)","TCR[%/°C](High I)","Rref[Ohm](Low I)","TCR[%/°C](Low I)","R[Ohm] (stress)","Temperature[°C](High I)","Temperature[°C](Low I)","Temperature[°C](Opt)","R[Ohm](Init Temp Low I)","R[Ohm](stress Temp Low I)","R[Ohm](stress Temp High I)")

    # Saving in a file
    FileName <- paste(substr(DegFileName, 1, nchar(DegFileName)-7),"exportfile.txt",sep="")
    write.table(NewFile,file=FileName,sep="\t",row.names=FALSE,quote=FALSE)
}


CreateExportFiles.Mira <- function(DegFileName,TCRFileName)
# Function called to create the exportfiles
# from a Qualitau Mira experiment.
# Use a Degradation and a TCR file to produce the exportfile
{
    # Device and Width parameters
    ListDevice <- read.delim("//fsup04/fntquap/Common/Qual/Process_Reliability/Process/amsReliability_R_Package/ListDeviceName.txt")

    DegFile <- read.delim(DegFileName,sep="")
    names(DegFile) <- c("#RESISTANCE#pkgNum", "ValueAtTimeZero","MinAchieved","MinUsedForFailure","FailureIteration","LastValue","NowValue","Split","FailureTime","FailedDevice","DUT","Extr1","Extr2","Extr3","Extr4","Ref")
    TCRFile <- read.delim(TCRFileName,sep="",skip=1)
    names(TCRFile) <- c("Pkg","Device","Rref(in Ohms)","TCR(in %/°C)","R(stress)","T(stress)","T(stress)-OvenT","Theta(°C/W)","R,50C,LowI","R,240C,LowI","R,240C,HiI")

    # Merge the files with additional data:
    # Split
    #Split<-c(1,2)
    # Istress extracted from filename and mA is removed
    Istress <- strsplit(DegFileName,split="_")[[1]][3]
    Istress <- as.numeric(substr(Istress, 1, nchar(Istress)-2))
    # Temp extracted from filename and C is removed
    Temp <- strsplit(DegFileName,split="_")[[1]][4]
    Temp <- as.numeric(substr(Temp, 1, nchar(Temp)-1))
    # DeviceID, Length and width
    L <- c(200,400)
    DeviceID <- strsplit(DegFileName,split="_")[[1]][2]
    W <- ListDevice$Width[ListDevice$Device==DeviceID]

    # Mira is made compatible with ACE for the 8 first colum. TTF are in hours. COnversion is made to seconds with 3600 factor. 
    NewFile <- data.frame(DegFile[DegFile$FailureIteration>0,1],DegFile[DegFile$FailureIteration>0,10],3600*DegFile[DegFile$FailureIteration>0,9],DegFile[DegFile$FailureIteration>0,8],Istress,L,W,Temp,DeviceID,DegFile[DegFile$FailureIteration>0,2],TCRFile[TCRFile$Device==4,3:8])
    #names(NewFile) <- c("#RESISTANCE#pkgNum","FailedDevice","FailureTime","Split","Istress","L","W","Temp","DeviceID","ValueAtTimeZero","Rref(in Ohms)","TCR(in %/°C)","R(stress)","T(stress)","T(stress)-OvenT","Theta(°C/W)")
    names(NewFile) <- c("Device","Failed","Lifetime[s]","Split","Istress","L","W","Temp","DeviceID","ValueAtTimeZero","Rref(in Ohms)","TCR(in %/°C)","R(stress)","T(stress)","T(stress)-OvenT","Theta(°C/W)")

    # Saving in a file
    #NewFile[,11:16] <-format(NewFile[,11:16], scientific = FALSE)
    FileName <- paste(substr(DegFileName, 1, nchar(DegFileName)-7),"exportfile.txt",sep="")
    write.table(NewFile,file=FileName,sep="\t",row.names=FALSE,quote=FALSE)

}


CreateExportFiles.EM <- function()
# Main function called to create the exportfiles
# from a set of Electromigration experiments.
# Open all deg and TCR files from the workfolder
{
    # Device and Width parameters
    ListDevice <- read.delim("//fsup04/fntquap/Common/Qual/Process_Reliability/Process/amsReliability_R_Package/ListDeviceName.txt")
    # File to be read
    ListDegFiles <- list.files(pattern="*deg.txt")
    ListTCRFiles <- list.files(pattern="*TCR.txt")

    # File number verification
    if (length(ListDegFiles) == 0 || length(ListTCRFiles) == 0 || length(ListDegFiles) != length(ListTCRFiles))  {
        print("Number of files is wrong. Please check!")
    } else {
        # We proceed

        for (i in 1:length(ListDegFiles)){

            # Check if the common part of the file name is the same. If not, nothing is produced.
            if (substr(ListDegFiles[i], 1, nchar(ListDegFiles[i])-7) == substr(ListTCRFiles[i], 1, nchar(ListDegFiles[i])-7)) {
              # Check if the structure is present in the list. If not user is warned
                DeviceID <- strsplit(ListDegFiles[i],split="_")[[1]][2]
                if (length(ListDevice$Width[ListDevice$Device==DeviceID])>0){

                    # Read the first line (headers) & distinguish between Mira and ACE
                    Headers <- scan(ListDegFiles[i], what="character", nlines = 1)
                    if (Headers[1] == "#RESISTANCE#pkgNum") {
                        # This is a Mira file
                        CreateExportFiles.Mira(ListDegFiles[i],ListTCRFiles[i])
                    } else {
                        # This is an ACE file
                        CreateExportFiles.Ace(ListDegFiles[i],ListTCRFiles[i])
                    }
                # Length is 0 thus the structure is not in the list
                } else {
                        print(paste("Structure",DeviceID, "is not present in the list. Please fill the list!"))
                }
            } else {
                print(paste("File names are not matching. Please verify!"))
            }
        }
    }
}
