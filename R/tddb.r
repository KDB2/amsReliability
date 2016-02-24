################################################################################
###                                                                          ###
###    INFORMATIONS                                                          ###
###    ---------------------------------                                     ###
###                                                                          ###
###       PACKAGE NAME        amsReliability                                 ###
###       MODULE NAME         tddb.r                                         ###
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
###       This module is dedicated to tddb experiments.                      ###
###    Extraction of model parameters is performed.                          ###
###                                                                          ###
###                                                                          ###
###    FUNCTIONS                                                             ###
###    ---------------------------------                                     ###
###                                                                          ###
###       OxideLifetimeModelization Extraction of oxide lifetime parameters  ###
###       OxideTDDB                 Main function for TDDB data analysis     ###
###       ReadDataTDDB              Read Exportfile and create data table    ###
###                                                                          ###
################################################################################


ReadDataTDDB <- function(ListFiles)
#Read the exportfiles listed in ListFile and store them in a dataframe.
# First read all the files and then calculate the probability scale
# for each condition. This allows to work with conditions splitted in different files.
# Data are cleaned to remove bad units
{

    # ResTable initialisation
    ResTable <- data.frame()

    for (FileName in ListFiles){

        # Find the line where the data are stored. They start after [DATA].
        StartLine <- grep('DATA]', readLines(FileName) )
        # Read the file. Data start two lines after [DATA]
        File <- read.delim(FileName, skip=StartLine+2, header=FALSE)

        # Store important data
        TTF <- File[,14]
        Temperature <- File[,9]
        Area <- File[,17]
        Stress <- File[,12]
        Status <- File[,15]
        # Make status similar to Electromigration: -1 wrong device, not to be considered in statistic; 1 failed; 0 not failed.
        # -1 bad sample set manually we leave them as it is.
        # 1 is failed before exp, we place them as -1
        Status[Status == 1] <- -1
        # 100 is still ongoing, it becomes a 0.
        Status[Status == 100] <- 0
        # other guys are considered as failed. They become 1.
        Status[Status != 0 & Status != -1] <- 1
        # Device dead during ramp are not considered in statistique (-1). They have a negative TTF
        Status[TTF < 0] <- -1
        # Device with Status -1 have wrong voltage
        Stress[Status==-1] <- as.numeric(substr(as.character(File[Status==-1,2]), 2, nchar(as.character(File[Status==-1,2]))))
        # Device with Status 0 have wrong voltage and TTF = 0. We will force good voltage and 1E30 for TTF.
        TTF[Status==0] <- 1E30
        Stress[Status==0] <- as.numeric(substr(as.character(File[Status==0,2]), 2, nchar(as.character(File[Status==0,2]))))
        # Condition stickers
        Conditions <- paste(Stress,"V/",Temperature,"°C",sep="")
        # Creation of a dataframe to store the data
        TempDataFrame <- data.frame(TTF,Status,Conditions,Stress,Temperature,Area)
        # Force the column names
        names(TempDataFrame) <- c("TTF", "Status", "Conditions", "Stress", "Temperature","Area")
        # Store the data in the final table
        ResTable <- rbind(ResTable,TempDataFrame)
    }

    # Now that all the files have been opened, let's check if different sizes were used.
    # If yes, we will change the Condition sticker to include the area.
    if (length(levels(factor(ResTable$Area))) > 1 ){
        ResTable$Conditions <- paste(ResTable$Conditions,"/",ResTable$Area /1000, "k" , sep="")
    }

    # Cleaning to remove units where status is not 1 or 0.
    ResTable <- Clean(ResTable)

    # List the conditions present in ResTable
    CondList <- levels(factor(ResTable$Conditions))

    # Probability is missing. Let's add it.
    # Final dataframe is ExpDataTable
    # Initialisation
    ExpDataTable <- data.frame()

    # For each condition found, we calculate the probability of failure. Data are stacked in ExpDataFrame. Weibull scale is used.
    for (condition in CondList){

        TempDataTable <- CreateDataFrame(ResTable$TTF[ResTable$Conditions==condition], ResTable$Status[ResTable$Conditions==condition],
          ResTable$Condition[ResTable$Conditions==condition], ResTable$Stress[ResTable$Conditions==condition], ResTable$Temperature[ResTable$Conditions==condition], Scale="Weibull",ResTable$Area[ResTable$Conditions==condition])

        ExpDataTable <- rbind(ExpDataTable,TempDataTable)
    }

    # We force the new names here as a security check. Area replace Dimension (Dimension is forced by CreateDataFrame)
    names(ExpDataTable) <- c("TTF", "Status", "Probability", "Conditions", "Stress", "Temperature","Area")
    # Order the condition in numerical/alphabetical order
    ExpDataTable <- ExpDataTable[OrderConditions(ExpDataTable),]
    # Do the same for conditions levels
    ExpDataTable$Conditions <- factor(ExpDataTable$Conditions, SortConditions(levels(ExpDataTable$Conditions)))


    return(ExpDataTable)

}


OxideLifetimeModelization <- function(DataTable, DeviceID)
# Modelize the data using a TDDB lifetime model
# Extract the parameters: t0, A and Ea
# as well as the Weibull slope
# TTF =
#
# Data(TTF,Status,Probability,Conditions,Stress,Temperature, Dimension)
# Data have to be cleaned upfront. Only valid data (status==1) should be given.
{
    # Modelization
    Model <- ModelFit(DataTable, Law="TDDB")

    # List the number of different area used during tests.
    ListArea <- levels(as.factor(DataTable$Area))

    # Table initialization
    ModelDataTable <- data.frame()

    # Model calculation for each area and for each condition.
    for (area in ListArea){
        ListConditions <- levels(DataTable$Conditions[DataTable$Area == area])
        ModelDataTable <- rbind(ModelDataTable, CreateModelDataTable(Model, ListConditions, as.numeric(area), Law="TDDB", Scale="Weibull"))
    }

    # if several length are used, stickers are adapted
    if (length(ListArea) != 1){
        ModelDataTable$Conditions <- paste(ModelDataTable$Conditions,"/", ModelDataTable$Area /1000, "k" , sep="")
    }

    # Display information about fit: parameters, goodness of fit...
    FitResultsDisplay(Model, DataTable, DeviceID)

    return(ModelDataTable)
}


#' Oxide breakdown data analysis
#'
#' Extract oxide lifetime parameters from a set of Time Dependant Dielectric
#' Breakdown (TDDB) experiments.
#' The experimental data as well as the resulting model are displayed and
#' can be saved. Extracted parameters are saved in a fit.txt file.
#'
#' @param ErrorBand displays the confidence intervals if set to TRUE.
#' @param ConfidenceValue percentage used in the confidence interval calculation
#' @param Save saves the chart as .png if set to TRUE.
#'
#' @return None
#'
#' @examples
#' OxideTDDB()
#' OxideTDDB(ErrorBand=FALSE)
#' @author Emmanuel Chery, \email{emmanuel.chery@@ams.com}
#' @import ggplot2 MASS scales nlstools tcltk
#' @export
OxideTDDB <- function(ErrorBand=FALSE, ConfidenceValue=0.95, Save=TRUE)
{
    # Disable warning for this function.
    oldw <- getOption("warn")
    options(warn = -1)
    #rm(list=ls())
    # ListFiles <- list.files(pattern="k_T.*txt$")

    # Filters for file selection
    Filters <- matrix(c("All files", "*", "Text", ".txt"),2, 2, byrow = TRUE)
    ListFiles <- SelectFilesAdvanced(Filters)
    #DeviceID <- strsplit(ListFiles[1],split="_")[[1]][2]
    # case 1, there are one or several files available
    if (length(ListFiles) != 0){
          # List of DeviceID available in the selected exportfiles
          DeviceIDList <- levels(sapply(ListFiles,function(x){factor(strsplit(x,split="_")[[1]][1])}))

          for (DeviceID in DeviceIDList){
              SubListFiles <- ListFiles[grep(DeviceID,ListFiles)]
              # Import the file(s) and create the 3 dataframes + display data
              DataTable <- try(ReadDataTDDB(SubListFiles), silent=TRUE)

              if (class(DataTable) != "try-error"){
                  # Reading the file was ok.

                  # Modelization, errorBands calculation and Graph is made with a clean table where only failed samples are kept.
                  # DataTable is kept in order to be saved in fit.txt
                  CleanExpDataTable <- KeepOnlyFailed(DataTable)

                  # Attempt to modelize. If succes, we plot the chart, otherwise we only plot the data.
                  ModelDataTable <- try(OxideLifetimeModelization(CleanExpDataTable, DeviceID),silent=TRUE)
                  # Check if the modelization is a succes
                  if (class(ModelDataTable) != "try-error"){
                      if (ErrorBand) {
                          ErrorDataTable <- ErrorEstimation(CleanExpDataTable, ModelDataTable, ConfidenceValue, Scale="Weibull")
                      } else {
                          ErrorDataTable <- NULL
                      }
                      CreateGraph(CleanExpDataTable, ModelDataTable, ErrorDataTable, aesVec = c("TTF", "Probability", "Conditions"), title = DeviceID,
                          axisTitles = c("Time to Failure (s)","Probability (%)"), scale.x = "Log", scale.y = "Weibull", save = Save)
                      # ExpData are added to the fit.txt file created during modelization
                      SaveData2File(DataTable, "fit.txt")
                  } else { # if modelization is not a success, we display the data and return parameters of the distribution in the console (scale and loc) in case user need them.
                      ModelDataTable <- FitDistribution(CleanExpDataTable,Scale="Weibull")
                      CreateGraph(CleanExpDataTable, ModelDataTable, aesVec = c("TTF", "Probability", "Conditions"), title = DeviceID,
                          axisTitles = c("Time to Failure (s)","Probability (%)"), scale.x = "Log", scale.y = "Weibull", save = FALSE)
                  }
              } else { # reading files returned an error
                  print("Error detected in the file(s) you selected. Please check your selection.")
              }
          }

    } else { # case 2, there are no files available
          print("You need to create the export files first!")
    }
    # return(DataTable)
    # Warning are set on again.
    options(warn = oldw)
}
