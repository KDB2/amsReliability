################################################################################
###                                                                          ###
###    INFORMATIONS                                                          ###
###    ---------------------------------                                     ###
###                                                                          ###
###       PACKAGE NAME        amsReliability                                 ###
###       MODULE NAME         exportFiles.r                                  ###
###       VERSION             0.10                                           ###
###                                                                          ###
###       AUTHOR              Emmanuel Chery                                 ###
###       MAIL                emmanuel.chery@ams.com                         ###
###       DATE                2016/02/24                                     ###
###       PLATFORM            Windows 7 & Gnu/Linux 3.16                     ###
###       R VERSION           R 3.1.1                                        ###
###       REQUIRED PACKAGES   ggplot2, grid, MASS, nlstools, scales          ###
###       LICENSE             GNU GENERAL PUBLIC LICENSE                     ###
###                           Version 3, 29 June 2007                        ###
###                                                                          ###
###                                                                          ###
###    DESCRIPTION                                                           ###
###    ---------------------------------                                     ###
###                                                                          ###
###       This package is a collection of scripts dedicated to help          ###
###    the process reliability team of ams AG. It includes tools to          ###
###    quickly visualize data and extract model parameters in order          ###
###    to predict device lifetimes.                                          ###
###                                                                          ###
###       This module includes functions used to automatize export           ###
###    files creation.                                                       ###
###                                                                          ###
###                                                                          ###
###    FUNCTIONS                                                             ###
###    ---------------------------------                                     ###
###                                                                          ###
###       CreateExportFiles.Ace     Subfunction for Ace equipement           ###
###       CreateExportFiles.Mira    Subfunction for Mira equipement          ###
###       CreateExportFiles         Main function for electromigration       ###
###                                                                          ###
################################################################################


CreateExportFiles.Ace <- function(DegFileName,TCRFileName)
# Function called to create the exportfiles
# from a Qualitau Ace experiment.
# Use a Degradation and a TCR file to produce the exportfile
{
    # Device and Width parameters
    ListDevice <- read.delim("//fsup04/fntquap/Common/Qual/Process_Reliability/Process/amsReliability_R_Package/ListDeviceName.txt")
    # List of exportfile already present in the folder
    ListExpFiles <- list.files(pattern="*exportfile.txt")

    DegFile <- read.delim(DegFileName,sep="")
    names(DegFile) = c("Device","Failed","Lifetime[s]","StressTemp[K]","Positive Current[A]","Negative Current[A]","Duty Cycle","Pulse Width","Stress Type")
    TCRFile <- read.delim(TCRFileName,sep="",skip=1)
    names(TCRFile) = c("Device","Rref[Ohm]","TCR[%/°C]","Rref[Ohm](High I)","TCR[%/°C](High I)","Rref[Ohm](Low I)","TCR[%/°C](Low I)","R[Ohm] (stress)","Temperature[°C](High I)","Temperature[°C](Low I)","Temperature[°C](Opt)","R[Ohm](Init Temp Low I)","R[Ohm](stress Temp Low I)","R[Ohm](stress Temp High I)", "Lot Number", "Wafer Number")

    # This section allows to work with degradation file where several conditions have been stored.
    # Creation of the condition vector. A condition is a couple current/Temperature
    Conditions <- factor(paste(DegFile[,5]*1E3,"mA/",DegFile[,4]-273.16,"°C",sep=""))
    # Add Conditions to both tables. Will be used to sort them.
    DegFile[,"Conditions"] <- Conditions
    TCRFile[,"Conditions"] <- Conditions

    # Add Split column
    DegFile[DegFile$Device<33,"Split"] <- 1
    DegFile[DegFile$Device>32,"Split"] <- 2
    # Add Length column
    DegFile[DegFile$Device<33,"L"] <- 200
    DegFile[DegFile$Device>32,"L"] <- 400

    # Creation of one export file per condition
    for (i in 1:length(levels(Conditions))){
        # Merge the files with additional data:
        # Split Istress Temp
        Split <- DegFile$Split[DegFile$Conditions == levels(Conditions)[i]]
        Istress <- DegFile[DegFile$Conditions == levels(Conditions)[i],5]*1E3 # mA
        Temp <- DegFile[DegFile$Conditions == levels(Conditions)[i],4]-273.16 # °C
        # DeviceID, Length and width
        L <- DegFile$L[DegFile$Conditions == levels(Conditions)[i]]
        DeviceID <- strsplit(DegFileName,split="_")[[1]][2]
        W <- ListDevice$Width[ListDevice$Device==DeviceID]

        # DataFrame creation
        NewFile <- data.frame(DegFile[DegFile$Conditions == levels(Conditions)[i],1:3],Split,Istress,L,W,Temp,DeviceID,TCRFile[TCRFile$Conditions == levels(Conditions)[i],2:16])
        names(NewFile) <- c("Device","Failed","Lifetime[s]","Split","Istress","L","W","Temp","DeviceID","Rref[Ohm]","TCR[%/°C]","Rref[Ohm](High I)","TCR[%/°C](High I)","Rref[Ohm](Low I)","TCR[%/°C](Low I)","R[Ohm] (stress)","Temperature[°C](High I)","Temperature[°C](Low I)","Temperature[°C](Opt)","R[Ohm](Init Temp Low I)","R[Ohm](stress Temp Low I)","R[Ohm](stress Temp High I)", "Lot Number", "Wafer Number")

        # Saving in a file
        FileName <- paste(strsplit(DegFileName,split="_")[[1]][1],"_",DeviceID,"_",Istress[1],"mA_",Temp[1],"C_exportfile.txt",sep="")
        # Let's check if the file exists already
        ShortFileName <- paste(strsplit(DegFileName,split="_")[[1]][1],"_",DeviceID,"_",Istress[1],"mA_",Temp[1],"C",sep="")

        if (length(grep(ShortFileName,ListExpFiles)) > 0) {

            UserChoice <- readline(prompt=paste("File ",FileName, " is already present. Keep[K], Replace[R], Merge[M]? (Default=K)",sep=""))

              if (UserChoice=="M" || UserChoice=="m"){
                  # Store the old File and stack it with NewFile if they are both from ACE
                  OldFile <- read.delim(FileName)

                  if ( length(names(OldFile)) == length(names(NewFile)) ) {
                      names(OldFile) <- names(NewFile)
                      NewFile <- rbind(OldFile,NewFile)
                      write.table(NewFile,file=FileName,sep="\t",row.names=FALSE,quote=FALSE)
                      print(paste("Data have been added to",FileName,sep=" "))
                  } else {
                      print("You are trying to merge data from ACE and MIRA. Impossible operation. Old file has been kept.")
                  }

              } else if (UserChoice=="R" || UserChoice=="r"){
                  write.table(NewFile,file=FileName,sep="\t",row.names=FALSE,quote=FALSE)
                  print(paste(FileName, "created.",sep=" "))
              } else {
                  print(paste("File ",FileName," was kept.",sep=""))
              }

        } else {
            write.table(NewFile,file=FileName,sep="\t",row.names=FALSE,quote=FALSE)
            print(paste(FileName, "created.",sep=" "))
        }
    }
}


CreateExportFiles.Mira <- function(DegFileName,TCRFileName)
# Function called to create the exportfiles
# from a Qualitau Mira experiment.
# Use a Degradation and a TCR file to produce the exportfile
{
    # Device and Width parameters
    ListDevice <- read.delim("//fsup04/fntquap/Common/Qual/Process_Reliability/Process/amsReliability_R_Package/ListDeviceName.txt")
    # List of exportfile already present in the folder
    ListExpFiles <- list.files(pattern="*exportfile.txt")

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
    DegFile[DegFile[,1]<31,"L"] <- 200
    DegFile[DegFile[,1]>30,"L"] <- 400
    DeviceID <- strsplit(DegFileName,split="_")[[1]][2]
    W <- ListDevice$Width[ListDevice$Device==DeviceID]
    # if MinUsedForFailure is -1 device is not fail -> 0
    DegFile[DegFile$MinUsedForFailure==-1,10] <- 0

    # Mira is made compatible with ACE for the 8 first colum. TTF are in hours. COnversion is made to seconds with 3600 factor.
    NewFile <- data.frame(DegFile[DegFile$ValueAtTimeZero>0,1],DegFile[DegFile$ValueAtTimeZero>0,10],3600*DegFile[DegFile$ValueAtTimeZero>0,9],DegFile[DegFile$ValueAtTimeZero>0,8],Istress,DegFile[DegFile$ValueAtTimeZero>0,17],W,Temp,DeviceID,DegFile[DegFile$ValueAtTimeZero>0,2],TCRFile[TCRFile$Device==4,3:8])
    #names(NewFile) <- c("#RESISTANCE#pkgNum","FailedDevice","FailureTime","Split","Istress","L","W","Temp","DeviceID","ValueAtTimeZero","Rref(in Ohms)","TCR(in %/°C)","R(stress)","T(stress)","T(stress)-OvenT","Theta(°C/W)")
    names(NewFile) <- c("Device","Failed","Lifetime[s]","Split","Istress","L","W","Temp","DeviceID","ValueAtTimeZero","Rref(in Ohms)","TCR(in %/°C)","R(stress)","T(stress)","T(stress)-OvenT","Theta(°C/W)")

    # Saving in a file
    FileName <- paste(substr(DegFileName, 1, nchar(DegFileName)-7),"exportfile.txt",sep="")
    # Let's check if the file exists already
    ShortFileName <- paste(strsplit(DegFileName,split="_")[[1]][1],"_",DeviceID,"_",Istress[1],"mA_",Temp[1],"C",sep="")

    if (length(grep(ShortFileName,ListExpFiles)) > 0) {

        UserChoice <- readline(prompt=paste("File ",FileName, " is already present. Keep[K], Replace[R], Merge[M]? (Default=K)",sep=""))

          if (UserChoice=="M" || UserChoice=="m"){
              # Store the old File and stack it with NewFile if they are both from MIRA
              OldFile <- read.delim(FileName)

              if ( length(names(OldFile)) == length(names(NewFile)) ) {
                  names(OldFile) <- names(NewFile)
                  NewFile <- rbind(OldFile,NewFile)
                  write.table(NewFile,file=FileName,sep="\t",row.names=FALSE,quote=FALSE)
                  print(paste("Data have been added to",FileName,sep=" "))
              } else {
                  print("You are trying to merge data from ACE and MIRA. Impossible operation. Old file has been kept.")
              }

          } else if (UserChoice=="R" || UserChoice=="r"){
              write.table(NewFile,file=FileName,sep="\t",row.names=FALSE,quote=FALSE)
              print(paste(FileName, "created.",sep=" "))
          } else {
              print(paste("File ",FileName," was kept.",sep=""))
          }

    } else {
        write.table(NewFile,file=FileName,sep="\t",row.names=FALSE,quote=FALSE)
        print(paste(FileName, "created.",sep=" "))
    }
}


#' Exportfiles creation
#'
#' Automatically creates the exportfiles for a set of electromigration experiments
#' by using degradation and TCR files provided by Qualitau ACE and MIRA equipments.
#' All files should be placed in the folder where the function is run.
#' The function matches file names to create the exportfiles but is able to detect
#' if several conditions are stored in one file. One exportfile is created for each condition.
#' If an exportfile is already present, user is asked if the old file is to be kept, replaced or if data have to be merged.
#' Times to failure are converted in seconds if they are provided in hours.
#'
#' @param None
#'
#' @return None
#'
#' @examples
#' CreateExportFiles()
#' @author Emmanuel Chery, \email{emmanuel.chery@@ams.com}
#' @export
CreateExportFiles <- function()
{
    # Device and Width parameters
    ListDevice <- try(read.delim("//fsup04/fntquap/Common/Qual/Process_Reliability/Process/amsReliability_R_Package/ListDeviceName.txt"),silent=TRUE)
    # File to be read
    # ListDegFiles <- list.files(pattern="*deg.txt")
    # ListTCRFiles <- list.files(pattern="*TCR.txt")
    listFilesSelected <- ExportFilesCleverSelection()
    ListDegFiles <- listFilesSelected[[1]]
    ListTCRFiles <- listFilesSelected[[2]]

    # Check if ListDeviceName.txt was available
    if (class(ListDevice)=="try-error"){
      print("File //fsup04/fntquap/Common/Qual/Process_Reliability/Process/amsReliability_R_Package/ListDeviceName.txt not found.")
    } else{

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
}


ExportFilesCleverSelection <- function()
# Handle in a clever way the selection of ExportFiles
# Will try to match files so that user can only select the deg or the TCR file.
# return a vector with both lists of files: list(listDeg, listTCR)
{
    filters <- matrix(c("All files", "*", "Text", ".txt", "TCR Files", "*TCR.txt", "deg Files", "*deg.txt"),4, 2, byrow = TRUE)
    fileSelected <- SelectFilesAdvanced(filters)

    # Remove the deg or TCR part and keep only one occurence of each name
    fileSelected <- unique(sapply(fileSelected,function(x){substr(x,1,nchar(x)-7)}))
    # add TCR and deg
    listTCR <- paste(fileSelected,"TCR.txt",sep="")
    listDeg <- paste(fileSelected,"deg.txt",sep="")
    # check if existant in the global list
    globalList <- list.files(pattern="*.txt")

    rangTCR <- lapply(listTCR,function(x){grep(x, globalList)}) # list of positions if present. 0 otherwise.
    rangTCR <- sapply(rangTCR, function(x) length(x) > 0) # Vector of booleans
    rangDeg <- lapply(listDeg,function(x){grep(x, globalList)})
    rangDeg <- sapply(rangDeg, function(x) length(x) > 0)

    # Inform user about missing files:
    if (sum(rangTCR == 0) > 0){
        print(paste("File", listTCR[!rangTCR], "is missing!", sep=" "))
    }

    if (sum(rangDeg == 0) > 0){
        print(paste("File", listDeg[!rangDeg], "is missing!", sep=" "))
    }

    # rangDeg & rangTCR return only the files where both the TCR and the deg files exists.
    listTCR <- listTCR[rangDeg & rangTCR]
    listDeg <- listDeg[rangDeg & rangTCR]
    return(list(listDeg, listTCR))
}
