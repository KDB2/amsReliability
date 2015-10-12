# Script collection for ams AG Process Reliability team.
# Allow the creation of exportfiles based on 
# data files provided by equipements.
# September 2015
# Emmanuel Chery
# Version 0.3 





CreateExportFiles.Ace <- function()
# Main function called to create the exportfiles.
# Open all deg and TCR files from the workfolder
{
    # Device and Width parameters
    ListDevice <- read.delim("//fsup04/fntquap/Common/Qual/Process_Reliability/Process/0.18_um_Technology/0.18 FabB/BEOL/Elmig/ListDeviceName.txt")
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
                DegFile <- read.delim(ListDegFiles[i],sep="")
                names(DegFile) <- c("Device","Failed","Lifetime[s]","StressTemp[K]","Positive Current[A]","Negative Current[A]","Duty Cycle","Pulse Width","Stress Type")
                TCRFile <- read.delim(ListTCRFiles[i],sep="",skip=1)
                names(TCRFile) <- c("Device","Rref[Ohm]","TCR[%/°C]","Rref[Ohm](High I)","TCR[%/°C](High I)","Rref[Ohm](Low I)","TCR[%/°C](Low I)","R[Ohm] (stress)","Temperature[°C](High I)","Temperature[°C](Low I)","Temperature[°C](Opt)","R[Ohm](Init Temp Low I)","R[Ohm](stress Temp Low I)","R[Ohm](stress Temp High I)")

                # Merge the files with new columns Split, Istress...
                Split<-c(1,2)
                # Istress extracted from filename and mA is removed
                Istress <- strsplit(ListDegFiles[i],split="_")[[1]][3]
                Istress <- as.numeric(substr(Istress, 1, nchar(Istress)-2))
                # Temp extracted from filename and C is removed
                Temp <- strsplit(ListDegFiles[i],split="_")[[1]][4]
                Temp <- as.numeric(substr(Temp, 1, nchar(Temp)-1))

                L <- c(200,400)
                DeviceID <- strsplit(ListDegFiles[i],split="_")[[1]][2]
                W <- ListDevice$Width[ListDevice$Device==DeviceID]

                NewFile <- data.frame(DegFile[,1:3],Split,Istress,L,W,Temp,DeviceID,TCRFile[,2:14])
                names(NewFile) <- c("Device","Failed","Lifetime[s]","Split","Istress","L","W","Temp","DeviceID","Rref[Ohm]","TCR[%/°C]","Rref[Ohm](High I)","TCR[%/°C](High I)","Rref[Ohm](Low I)","TCR[%/°C](Low I)","R[Ohm] (stress)","Temperature[°C](High I)","Temperature[°C](Low I)","Temperature[°C](Opt)","R[Ohm](Init Temp Low I)","R[Ohm](stress Temp Low I)","R[Ohm](stress Temp High I)")

                # Saving in a file
                FileName <- paste(substr(ListDegFiles[i], 1, nchar(ListDegFiles[i])-7),"exportfile.txt",sep="")
                write.table(NewFile,file=FileName,sep="\t",row.names=FALSE,quote=FALSE)
            }
        }
    }
}
