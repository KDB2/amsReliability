################################################################################
###                                                                          ###
###    INFORMATIONS                                                          ###
###    ---------------------------------                                     ###
###                                                                          ###
###       PACKAGE NAME        amsReliability                                 ###
###       MODULE NAME         tddb.r                                         ###
###       VERSION             0.8                                            ###
###                                                                          ###
###       AUTHOR              Emmanuel Chery                                 ###
###       MAIL                emmanuel.chery@ams.com                         ###
###       DATE                2015/11/10                                     ###
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
# Read the exportfiles listed in ListFile and store them in a dataframe.
# First read all the files and then calculate the probability scale
# for each condition. This allows to work with conditions splitted in different files.
# Data are cleaned to remove bad units
{

    # ResTable initialisation
    ResTable <- data.frame()

    for (i in seq_along(ListFiles)){

        # Find the line where the data are stored. They start after [DATA].
        StartLine <- grep('DATA]', readLines(ListFiles[i]) )
        # Read the file. Data start two lines after [DATA]
        File <- read.delim(ListFiles[i], skip=StartLine+2, header=FALSE)

        # Store important data
        TTF <- File[,14]
        Temperature <- File[,9]
        Dimension <- File[,17]
        Stress <- File[,12]
        Status <- File[,15]
        # Make status similar to Electromigration: -1 wrong device, not to be considered in statistic; 1 failed; 0 not failed.
        Status[Status != 100 & Status != -1] <- 1
        Status[Status == 100] <- 0
        # Device dead during ramp are not considered in statistique (-1). They have a negative TTF
        Status[TTF < 0] <- -1
        # Device with Status 0 have wrong voltage and TTF = 0. We will force good voltage and 1E30 for TTF.
        TTF[Status==0] <- 1E30
        Stress[Status==0] <- as.numeric(substr(as.character(File[Status==0,2]), 2, nchar(as.character(File[Status==0,2]))))
        # Condition stickers
        Conditions <- paste(Stress,"V/",Temperature,"C",sep="")
        # Creation of a dataframe to store the data
        TempDataFrame <- data.frame(TTF,Status,Conditions,Stress,Temperature,Dimension)
        # Force the column names
        names(TempDataFrame) <- c("TTF", "Status", "Conditions", "Stress", "Temperature","Dimension")
        # Store the data in the final table
        ResTable <- rbind(ResTable,TempDataFrame)
    }

    # Now that all the files have been opened, let's check if different sizes were used.
    # If yes, we will change the Condition sticker to include the area.
    if (length(levels(factor(ResTable$Dimension))) > 1 ){
        ResTable$Conditions <- paste(ResTable$Conditions,"/",ResTable$Dimension /1000, "k" , sep="")
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
    for (i in seq_along(CondList)){

        TempDataTable <- CreateDataFrame(ResTable$TTF[ResTable$Conditions==CondList[i]], ResTable$Status[ResTable$Conditions==CondList[i]],
          ResTable$Condition[ResTable$Conditions==CondList[i]], ResTable$Stress[ResTable$Conditions==CondList[i]], ResTable$Temperature[ResTable$Conditions==CondList[i]], Scale="Weibull",ResTable$Dimension[ResTable$Conditions==CondList[i]])

        ExpDataTable <- rbind(ExpDataTable,TempDataTable)
    }

    # We force the new names here as a security check.
    names(ExpDataTable) <- c("TTF", "Status", "Probability", "Conditions", "Stress", "Temperature","Dimension")
    return(ExpDataTable)

}


OxideLifetimeModelization <- function(DataTable,DeviceID)
# Modelize the data using a TDDB lifetime model
# Extract the parameters: t0, A and Ea
# as well as the Weibull slope
# TTF =
#
# Data(TTF,Status,Probability,Conditions,Stress,Temperature, Dimension)
{

    Area <- DataTable$Dimension * 1E-12 # m^2

    # Physical constants
    k <- 1.38E-23 # Boltzmann
    e <- 1.6E-19 # electron charge

    # Remove the units where status is 0
    CleanDataTable <- DataTable[DataTable$Status==1,]

    # Black model / Log scale: use of log10 to avoid giving too much importance to data with a high TTF
    Model <- nls(log10(TTF) ~ log10(exp(t0)*exp(-g*Stress)*exp((Ea*e)/(k*(Temperature+273.15)))*Area^(-1/beta)*exp(Probability/beta)), CleanDataTable, start=list(t0=30,g=1,Ea=0.2,beta=1),control= list(maxiter = 50, tol = 1e-6))#, minFactor = 1E-5, printEval = FALSE, warnOnly = FALSE))#,trace = T)
    # Parameters Extraction
    t0 <- coef(Model)[1]
    g <- coef(Model)[2]
    Ea <-coef(Model)[3]
    beta <- coef(Model)[4]
    # Residual Sum of Squares
    RSS <- sum(resid(Model)^2)
    # Total Sum of Squares: TSS <- sum((TTF - mean(TTF))^2))
    TSS <- sum(sapply(split(CleanDataTable[,1],CleanDataTable$Conditions),function(x) sum((x-mean(x))^2)))
    Rsq <- 1-RSS/TSS # R-squared measure

    # Using the parameters and the conditions, theoretical distributions are created
    ListConditions <- levels(CleanDataTable$Conditions)

    # Initialisation
    ModelDataTable <- data.frame()
    # y axis points are calculated. (limits 0.01% -- 99.99%) Necessary to have nice confidence bands.
    Proba <-  seq(log(-log(1-0.001)),log(-log(1-0.9999)),0.05) # Weibit

    for (i in seq_along(ListConditions)){
        # Experimental conditions:
        Condition <- ListConditions[i]
        V <- CleanDataTable$Stress[CleanDataTable$Conditions==Condition][1]
        Temp <- CleanDataTable$Temperature[CleanDataTable$Conditions==Condition][1]  # °C
        Area <- CleanDataTable$Dimension[CleanDataTable$Conditions==Condition][1] * 1E-12 # m^2

        # TTF calculation
        TTF <- exp(t0)*exp(-g*V)*exp((Ea*e)/(k*(Temp+273.15)))*Area^(-1/beta)*exp(Proba/beta)

        # Dataframe creation
        ModelDataTable <- rbind(ModelDataTable, data.frame('TTF'=TTF,'Status'=1,'Probability'=Proba,'Conditions'=Condition,'Stress'=V,'Temperature'=Temp))
    }

    # Drawing of the residual plots
    plot(nlsResiduals(Model))
    # Display of fit results
    print(DeviceID)
    print(summary(Model))
    print(paste("Residual squared sum: ",RSS,sep=""))
    #print(coef(Model))
    #print(sd(resid(Model)))

    # Save in a file
    capture.output(summary(Model),file="fit.txt")
    cat("Residual Squared sum:\t",file="fit.txt",append=TRUE)
    cat(RSS,file="fit.txt",append=TRUE)
    cat("\n \n",file="fit.txt",append=TRUE)
    cat("Experimental Data:",file="fit.txt",append=TRUE)
    cat("\n",file="fit.txt",append=TRUE)
    capture.output(DataTable,file="fit.txt",append=TRUE)
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
#' @import ggplot2 MASS scales grid nlstools
#' @export
OxideTDDB <- function(ErrorBand=FALSE, ConfidenceValue=0.95, Save=TRUE)
{
    # Disable warning for this function.
    oldw <- getOption("warn")
    options(warn = -1)
    #rm(list=ls())
    ListFiles <- list.files(pattern="k_T.*txt$")
    #DeviceID <- strsplit(ListFiles[1],split="_")[[1]][2]
    # case 1, there are one or several files available
    if (length(ListFiles) != 0){
          # List of DeviceID available in the selected exportfiles
          DeviceID <- levels(sapply(ListFiles,function(x){factor(strsplit(x,split="_")[[1]][1])}))

          for (i in seq_along(DeviceID)){
              SubListFiles <- ListFiles[grep(DeviceID[i],ListFiles)]
              # Import the file(s) and create the 3 dataframes + display data
              DataTable <- ReadDataTDDB(SubListFiles)
              # Attempt to modelize. If succes, we plot the chart, otherwise we only plot the data.
              ModelDataTable <- try(OxideLifetimeModelization(DataTable, DeviceID[i]),silent=TRUE)
              # Check if the modelization is a succes
              if (class(ModelDataTable) != "try-error"){
                    ErrorDataTable <- ErrorEstimation(DataTable, ModelDataTable, ConfidenceValue)
                    CreateGraph(DataTable,ModelDataTable,ErrorDataTable,DeviceID[i],Scale="Weibull",ErrorBand,Save)
              } else { # if modelization is not a success, we display the data and return parameters of the distribution in the console (scale and loc) in case user need them.
                    ModelDataTable <- FitDistribution(DataTable,Scale="Weibull")
                    CreateGraph(DataTable,ModelDataTable,DataTable,DeviceID[i],Scale="Weibull",ErrorBand=FALSE,Save=FALSE)
              }
          }

    } else { # case 2, there are no files available
          print("You need to create the export files first!")
    }
    return(DataTable)
    # Warning are set on again.
    options(warn = oldw)
}
