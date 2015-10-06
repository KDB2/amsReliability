# Script collection for the process reliability team of ams.
# First part is allowing the extraction of Black's parameters
# using a serie of electromigration experiments.
# September 2015
# Emmanuel Chery
# Version 0.8 RC

# Required Packages
library('ggplot2')
library('MASS')
library('scales')
library('grid')



Ranking <- function(TTF)
# Fraction estimator calculation
# rk(i)=(i-0.3)/(n+0.4)
# TTF is a vector.
{
    # ties.method="random" handles identical TTFs and provide a unique ID
    rk <- (rank(TTF, ties.method="random")-0.3)/(length(TTF)+0.4)
}


CalculProbability <- function(Probability, Scale="Lognormale")
# Given a vector Probability of probabilities, the function calculates
# the correspondence in standard deviations for the lognormale case.
# Calculation of the Weibit is made for the Weibull case.
{
  if (Scale=="Weibull") {
      Proba <- log(-log(1-Probability)) # Weibull
  } else {
      Proba <- qnorm(Probability) # Lognormal
  }
  return(Proba)
}


Clean <- function(DataTable)
# Take a datatable provided by CreateDataFrame and clean it.
# Cleaning of the data. Only lines with a status 1 or 0 are kept.
# Finally TTF are sorted from the smallest to the largest.
# Qualitau column name is 'Failed'
{
    CleanedTable <- DataTable[DataTable$Status==1 | DataTable$Status==0,]
    CleanedTable <- CleanedTable[order(CleanedTable$"TTF"),] # Sort TTF
    return(CleanedTable)
}


CreateDataFrame <- function(TTF, Status, Condition, Stress, Temperature, Scale="Lognormale")
# Creation of the dataframe assembling the TTF, the status of the samples,
# the probability, the condition (stickers for charts),
# the stress condition and the temperature used durng the stress.
# The probability is calculated according to lognormale or Weibull distribution.
# Data are given clean.

# Data(TTF,Status,Probability,Conditions,Stress,Temperature)
{
    rk <- Ranking(TTF) # Fraction estimator calculation
    if (Scale=="Weibull") {
        Proba <- CalculProbability(rk,Scale="Weibull") # Probability calculation Weibull
    } else {
        Proba <- CalculProbability(rk,Scale="Lognormale") # Probability calculation Lognormale
    }
    # Generation of the final data frame
    DataTable <- data.frame('TTF'=TTF,'Status'=Status,'Probability'=Proba,'Conditions'=Condition, 'Stress'=Stress, 'Temperature'=Temperature)
    return(DataTable)
}


ReadDataAce <- function(FileName, Scale="Lognormale")
# Read the file exportfile and store it in a dataframe
# Data are cleaned to remove bad units
{
    # Read the file and store it
    ResTable <- read.delim(FileName)
    # Creation of the new dataframe
    TTF <- ResTable["Lifetime.s."]
    Status <- ResTable["Failed"]
    Stress <- ResTable["Istress"]
    Temperature <- ResTable["Temp"]
    Condition <- paste(ResTable[,"Istress"],"mA/",ResTable[,"Temp"],"°C",sep="") #paste(ResTable[,5],"mA/",ResTable[,8],"C",sep="")
    ResTable <- data.frame(TTF,Status,Condition,Stress,Temperature)
    # Force the column names
    names(ResTable) <- c("TTF", "Status", "Conditions", "Stress", "Temperature")
    # Cleaning
    ResTable <- Clean(ResTable)

    # Probability is missing. Let's add it.
    if (Scale=="Weibull") {
        ExpDataTable <- CreateDataFrame(ResTable$TTF, ResTable$Status, ResTable$Condition, ResTable$Stress, ResTable$Temperature, Scale="Weibull")
    } else {
        ExpDataTable <- CreateDataFrame(ResTable$TTF, ResTable$Status, ResTable$Condition, ResTable$Stress, ResTable$Temperature, Scale="Lognormale")
    }
    # We force the new names here as a security check.
    names(ExpDataTable) <- c("TTF", "Status", "Probability", "Conditions", "Stress", "Temperature")
    return(ExpDataTable)
}


Modelization <- function(DataTable, Type="Lognormale")
# Using experimental data the theoretical distribution is generated
# Default is Lonormale scale but Weibull is available as an option.
{
    # Condition, Stress and Temperature stickers
    ModelCondition <- DataTable[1,"Conditions"]
    ModelStress <- DataTable[1,"Stress"]
    ModelTemperature <- DataTable[1,"Temperature"]

    # x axis limits are calculated
    lim <- range(DataTable$TTF)
    lim.high <- 10^(ceiling(log(lim[2],10)))
    lim.low <- 10^(floor(log(lim[1],10)))
    # Generation of a vector for the calculation of the model. 200pts/decades
    x <- 10^seq(log(lim.low,10),log(lim.high,10),0.005)
    # Model calculation with the experimental TTF
    if (Type=="Weibull") { # Weibull
          #x <- 10^seq(log(lim.low,10),log(lim.high,10),0.00005)
          fit <- fitdistr(DataTable$TTF[DataTable$Status==1],"weibull")
          Shape <- fit$estimate[1]  # Beta
          Scale <- fit$estimate[2]  # Characteristic time (t_63%)
          y <- CalculProbability(pweibull(x, Shape, Scale),"Weibull")
    } else { # Lognormale
          fit <- fitdistr(DataTable$TTF[DataTable$Status==1],"lognormal")
          Scale <- fit$estimate[1]  # meanlog
          Shape <- fit$estimate[2]  # sdlog
          y <- CalculProbability(plnorm(x, Scale, Shape),"Lognormale")
    }

    ModelDataTable <- data.frame('TTF'=x,'Status'=1,'Probability'=y,'Conditions'=ModelCondition,'Stress'=ModelStress,'Temperature'=ModelTemperature)
    return(ModelDataTable)
}


ErrorEstimation <- function(ExpDataTable, ModelDataTable, ConfidenceValue=0.95)
# Genration of confidence intervals
{
    NbData <- length(ExpDataTable$TTF)
    if (NbData > 30) {
        mZP_Value <- qnorm((1 - ConfidenceValue) / 2) # Normal case. Valid if sample size > 30.
    } else {
        mZP_Value <- qt((1 - ConfidenceValue) / 2, df=(NbData -1) ) # t-test statistic for low sample size
    }
    CDF <- pnorm(ModelDataTable$Probability) # TO BE CHECKED
    sef <- sqrt(CDF * (1 - CDF)/NbData) # TO BE CHECKED
    LowerLimit <- qnorm(CDF - sef * mZP_Value)
    HigherLimit <- qnorm(CDF + sef * mZP_Value)

    ConfidenceDataTable <- data.frame('TTF'=ModelDataTable$TTF,'LowerLimit'=LowerLimit,'HigherLimit'=HigherLimit,'Conditions'=ModelDataTable$Conditions)
    return(ConfidenceDataTable)
}


StackData <- function(DataTable1, DataTable2)
# Merge 2 DataTable
{
    NewDataTable <- merge(DataTable1, DataTable2, all=TRUE)
    NewDataTable <- NewDataTable[order(NewDataTable$"Conditions"),]
    return(NewDataTable)
}


CreateGraph <- function(ExpDataTable, ModelDataTable, ConfidenceDataTable, Title="", Scale="Lognormale", ErrorBands=TRUE)
# Use the table prepared with CreateDataFrame and create the probability plot.
# Default is Lonormale scale but Weibull is available as an option.
{
    # x scale limits calculation based on the data.
    lim <- range(ExpDataTable$TTF[ExpDataTable$Status==1]) # Min of the values is stored in [1] and max in  [2]
    lim.high <- 10^(ceiling(log(lim[2],10)))
    lim.low <- 10^(floor(log(lim[1],10)))
    # Now that we have the limits, we create the graph labels for x axis.
    GraphLabels <- 10^(seq(floor(log(lim[1],10)),ceiling(log(lim[2],10))))

    # Label for y axis
    # Dynamique labels as a function of the minimal probability observed.
    # Minimal proba is 0.01 %

    #  Weibull
    if (Scale == "Weibull") {
        if (ExpDataTable[1,"Probability"]<= CalculProbability(0.1/100,Scale)){ # Case 1: lower than 0.1%
            ListeProba <- c(0.01,0.1,1,2,3,5,10,20,30,40,50,63,70,80,90,95,99,99.9,99.99)
        }
        if (ExpDataTable[1,"Probability"]<= CalculProbability(1/100,Scale) && ExpDataTable[1,"Probability"]>= CalculProbability(0.1/100,Scale)){ # Case 2: lower than 1% but higher than 0.1%
            ListeProba <- c(0.1,1,2,3,5,10,20,30,40,50,63,70,80,90,95,99,99.9)
        }
        if (ExpDataTable[1,"Probability"] >= CalculProbability(1/100,Scale)) { # Case 3: higher than 1%
            ListeProba <- c(1,2,3,5,10,20,30,40,50,63,70,80,90,95,99)
        }
    ProbaNorm <- CalculProbability(ListeProba/100,Scale)

    } else { # Lognormale
        if (ExpDataTable[1,"Probability"]<= qnorm(0.1/100)){ # Case 1: lower than 0.1%
            ListeProba <- c(0.01,0.1,1,5,10,20,30,40,50,60,70,80,90,95,99,99.9,99.99)
        }
        if (ExpDataTable[1,"Probability"]<= qnorm(1/100) && ExpDataTable[1,"Probability"]>= qnorm(0.1/100)){ # Case 2: lower than 1% but higher than 0.1%
            ListeProba <- c(0.1,1,5,10,20,30,40,50,60,70,80,90,95,99,99.9)
        }
        if (ExpDataTable[1,"Probability"] >= qnorm(1/100)) { # Case 3: higher than 1%
            ListeProba <- c(1,5,10,20,30,40,50,60,70,80,90,95,99)
        }
    ProbaNorm <- qnorm(ListeProba/100)
    }

    # We are only going to plot samples where status is '1' (experiment is finished).
    # Table is sorted & conditions stay togeteher.
    CleanExpTable <- ExpDataTable[ExpDataTable$Status==1,]
    CleanExpTable <- CleanExpTable[order(CleanExpTable$"Conditions"),]

    # Graph creation with CleanTable
    Graph <- ggplot(data=CleanExpTable, aes(x=TTF, y=Probability, colour=Conditions, shape=Conditions))
    # box around chart + background
    Graph <- Graph + theme_linedraw() + theme(panel.background = element_rect(fill="gray90", color="black"))
    # Grid definitions
    Graph <- Graph + theme(panel.grid.major = element_line(colour="white", size=0.25, linetype=1))
    Graph <- Graph + theme(panel.grid.minor = element_line(linetype=0, colour="white", size = 0.25))
    # Definition of scales
    Graph <- Graph + scale_x_log10(limits = c(lim.low,lim.high),breaks = GraphLabels,labels = trans_format("log10", math_format(10^.x)))
    Graph <- Graph + scale_y_continuous(limits=range(ProbaNorm), breaks=ProbaNorm, labels=ListeProba )
    # Controled symbol list -- Max is 20 conditions on the chart.
    Graph <- Graph + scale_shape_manual(values=c(19,15,17,16,19,15,17,16,19,15,17,16,19,15,17,16,19,15,17,16))
    Graph <- Graph + scale_colour_manual(values = c("#d53e4f","#3288bd","#66a61e","#f46d43","#e6ab02","#8073ac","#a6761d","#666666","#bc80bd","#d53e4f","#3288bd","#66a61e","#f46d43","#e6ab02","#8073ac","#a6761d","#666666","#bc80bd","#d53e4f","#3288bd")) # "#5e4fa2" ,"#66c2a5", "#fec44f",
    Graph <- Graph + geom_point(size=4)+annotation_logticks(sides='tb')
    # Add the theoretical model
    Graph <- Graph + geom_line(data=ModelDataTable, aes(color=Conditions), size=0.8)
    # Add the confidence intervals
    if (ErrorBands==TRUE) {
        Graph <- Graph + geom_line(data=ConfidenceDataTable, aes(x=TTF, y=LowerLimit, color=Conditions), linetype="dashed", size=0.8)
        Graph <- Graph + geom_line(data=ConfidenceDataTable, aes(x=TTF, y=HigherLimit, color=Conditions), linetype="dashed",size=0.8)
    }
    # Font size & x/y titles...
    Graph <- Graph + xlab("Time to Failure (s)") + ylab("Probability (%)")
    Graph <- Graph + theme(axis.title.x = element_text(face="bold", size=16))
    Graph <- Graph + theme(axis.title.y = element_text(face="bold", size=16))
    # legend size
    Graph <- Graph + theme(legend.title = element_text(size=14, face="bold"))
    Graph <- Graph + theme(legend.text = element_text(size = 12))
    # Box around legend
    Graph <- Graph + theme(legend.background = element_rect())
    Graph <- Graph + theme(legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"))
    #Box around the conditions in legend
    Graph <- Graph + theme(legend.key = element_rect(fill="gray90", colour = "black", linetype=0))
    # Label/ticks size
    Graph <- Graph + theme(axis.text.x = element_text(face="bold", size=16))
    Graph <- Graph + theme(axis.text.y = element_text(size=16))
    Graph <- Graph + theme(axis.ticks.length = unit(-.25, "cm"), axis.ticks.margin=unit(0.4, "cm"))
    # Add a title
    Graph <- Graph + ggtitle(Title)
    Graph <- Graph + theme(plot.title = element_text(face="bold", size=18))

    print(Graph)

    # Save as png
    if (Title != ""){
        ggsave(paste(Title,"png",sep="."))
    } else {
        ggsave("Chart.png")
    }
}


BlackAnalysis <- function(Scale="Lognormale",ErrorBand=TRUE)
# Main function calling the other. The one to use to open all the files.
# Open all the exportfiles from the workfolder
{
    #rm(list=ls())
    ListFiles <- list.files(pattern="*exportfile.txt")
    StructureName <- strsplit(ListFiles[1],split="_")[[1]][2]
    # case 1, there are one or several files available
    if (length(ListFiles) != 0){
          # Import the first file to create the 3 dataframes
          DataTable <- ReadDataAce(ListFiles[1],Scale)
          ModelDataTable <- Modelization(DataTable,Scale)
          ErrorDataTable <- ErrorEstimation(DataTable,ModelDataTable)

          # Let's now check if other files are available
          if (length(ListFiles) > 1){
                # loop to open all the files and stack them in the dataframe
                for (i in 2:length(ListFiles)){
                    NewDataTable <- ReadDataAce(ListFiles[i],Scale)
                    NewModelDataTable <- Modelization(NewDataTable,Scale)
                    NewErrorDataTable <- ErrorEstimation(NewDataTable,NewModelDataTable)

                    # Merging the tables
                    DataTable <- StackData(DataTable,NewDataTable)
                    ModelDataTable <- StackData(ModelDataTable,NewModelDataTable)
                    ErrorDataTable <- StackData(ErrorDataTable,NewErrorDataTable)
                }
          }
    } else { # case 2, there are no files available
          print("You need to create the export files first!")
    }
    CreateGraph(DataTable,ModelDataTable,ErrorDataTable,StructureName,Scale,ErrorBand)
    #return(DataTable)
}


CreateExportFiles <- function()
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
