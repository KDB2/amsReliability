# Script collection for ams AG process reliability team.
# Allow to display and process electromigration data.
# Extraction of Black's parameters is performed.
# September 2015
# Emmanuel Chery
# Version 0.6


Ranking <- function(TTF)
# Fraction estimator calculation
# rk(i)=(i-0.3)/(n+0.4)
# TTF is a vector.
{
    # ties.method="random" handles identical TTFs and provide a unique ID
    rk <- (rank(TTF, ties.method="random")-0.3)/(length(TTF)+0.4)
}


CalculProbability <- function(Probability, Scale="Lognormal")
# Given a vector Probability of probabilities, the function calculates
# the correspondence in standard deviations for the Lognormal case.
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


CreateDataFrame <- function(TTF, Status, Condition, Stress, Temperature, Scale="Lognormal")
# Creation of the dataframe assembling the TTF, the status of the samples,
# the probability, the condition (stickers for charts),
# the stress condition and the temperature used durng the stress.
# The probability is calculated according to Lognormal or Weibull distribution.
# Data are given clean.

# Data(TTF,Status,Probability,Conditions,Stress,Temperature)
{
    rk <- Ranking(TTF) # Fraction estimator calculation
    if (Scale=="Weibull") {
        Proba <- CalculProbability(rk,Scale="Weibull") # Probability calculation Weibull
    } else {
        Proba <- CalculProbability(rk,Scale="Lognormal") # Probability calculation Lognormal
    }
    # Generation of the final data frame
    DataTable <- data.frame('TTF'=TTF,'Status'=Status,'Probability'=Proba,'Conditions'=Condition, 'Stress'=Stress, 'Temperature'=Temperature)
    return(DataTable)
}


ReadDataAce <- function(ListFileName, Scale="Lognormal")
# Read the exportfiles listed in ListFileName and store them in a dataframe.
# First read all the files and then calculate the probability scale
# for each condition. This allows to work with conditions splitted in different files.
# Data are cleaned to remove bad units
# Exportfile from Ace and Mira have different headers,
# therefore column numbers are used
{
    # ResTable initialisation
    ResTable <- data.frame()

    for (i in seq_along(ListFileName)){

        # Read the file and store it in a temporary dataframe
        TempTable <- read.delim(ListFileName[i])
        # Creation of the new dataframe
        TTF <- TempTable[,3]
        Status <- TempTable[,2]
        Stress <- TempTable[,5]
        Temperature <- TempTable[,8]
        Condition <- paste(TempTable[,5],"mA/",TempTable[,8],"C",sep="") #paste(TempTable[,"Istress"],"mA/",TempTable[,"Temp"],"°C",sep="")
        # Creation of a dataframe to store the data
        TempDataFrame <- data.frame(TTF,Status,Condition,Stress,Temperature)
        # Force the column names
        names(TempDataFrame) <- c("TTF", "Status", "Conditions", "Stress", "Temperature")
        # Store the data in the final table
        ResTable <- rbind(ResTable,TempDataFrame)
    }
    # security check, we force again the name on ResTable
    names(ResTable) <- c("TTF", "Status", "Conditions", "Stress", "Temperature")

    # Cleaning to remove units where status is not 1 or 0.
    ResTable <- Clean(ResTable)

    # List the conditions present in ResTable
    CondList <- levels(factor(ResTable$Conditions))

    # Probability is missing. Let's add it.
    # Final dataframe is ExpDataTable
    # Initialisation
    ExpDataTable <- data.frame()

    for (i in seq_along(CondList)){
        if (Scale=="Weibull") {
          TempDataTable <- CreateDataFrame(ResTable$TTF[ResTable$Conditions==CondList[i]], ResTable$Status[ResTable$Conditions==CondList[i]],
            ResTable$Condition[ResTable$Conditions==CondList[i]], ResTable$Stress[ResTable$Conditions==CondList[i]], ResTable$Temperature[ResTable$Conditions==CondList[i]], Scale="Weibull")
        } else {
          TempDataTable <- CreateDataFrame(ResTable$TTF[ResTable$Conditions==CondList[i]], ResTable$Status[ResTable$Conditions==CondList[i]],
            ResTable$Condition[ResTable$Conditions==CondList[i]], ResTable$Stress[ResTable$Conditions==CondList[i]], ResTable$Temperature[ResTable$Conditions==CondList[i]], Scale="Lognormal")
        }
        ExpDataTable <- rbind(ExpDataTable,TempDataTable)
    }

    # We force the new names here as a security check.
    names(ExpDataTable) <- c("TTF", "Status", "Probability", "Conditions", "Stress", "Temperature")
    return(ExpDataTable)
}


BlackModelization <- function(DataTable, DeviceID)
# Modelize the data using Black equation
# Extract the parameters: A, n and Ea
# as well as the lognormal slope
# TTF = A j^(-n) exp(Ea/kT + Scale * Proba)
# Proba in standard deviations
# Data(TTF,Status,Probability,Conditions,Stress,Temperature)
{
    # Read the list of device to retrieve the section parameters.
    ListDevice <- read.delim("//fsup04/fntquap/Common/Qual/Process_Reliability/Process/amsReliability_R_Package/ListDeviceName.txt")
    W <- ListDevice$Width[ListDevice$Device==DeviceID] # micrometers
    H <- ListDevice$Height[ListDevice$Device==DeviceID] # micrometers
    S <- W*H*1E-12 # m^2

    # if S is a positive number different from 0, we can proceed:
    if (is.na(S) || S<=0 ) {
        print(paste("Structure",DeviceID, "is not present in the list. Please fill the list!"))
        return(DataTable)
    } else { # we proceed

      # Physical constants
      k <- 1.38E-23 # Boltzmann
      e <- 1.6E-19 # electron charge

      # Remove the units where status is 0
      CleanDataTable <- DataTable[DataTable$Status==1,]

      # Black model / Log scale: use of log10 to avoid giving too much importance to data with a high TTF
      Model <- nls(log10(TTF) ~ log10(exp(A)*(Stress*1E-3/S)^(-n)*exp((Ea*e)/(k*(Temperature+273.15))+Scale*Probability)), CleanDataTable, start=list(A=30,n=1,Ea=0.7,Scale=0.3),control= list(maxiter = 50, tol = 1e-7))#, minFactor = 1E-5, printEval = FALSE, warnOnly = FALSE))#,trace = T)
      #Model <- nls(TTF ~ exp(A)*(Stress*1E-3/S)^(-n)*exp((Ea*e)/(k*(Temperature+273.15))+Scale*Probability), DataTable, start=list(A=30,n=1,Ea=0.7,Scale=0.3))
      # Parameters Extraction
      A <- coef(Model)[1]
      n <- coef(Model)[2]
      Ea <-coef(Model)[3]
      Scale <- coef(Model)[4]
      # Residual Sum of Squares
      RSS <- sum(resid(Model)^2)
      # Total Sum of Squares: TSS <- sum((TTF - mean(TTF))^2))
      TSS <- sum(sapply(split(CleanDataTable[,1],CleanDataTable$Conditions),function(x) sum((x-mean(x))^2)))
      Rsq <- 1-RSS/TSS # R-squared measure
      #print(paste("Size on 150 rows:", format(object.size(Model), unit="Mb")))

      # Using the parameters and the conditions, theoretical distributions are created
      ListConditions <- levels(CleanDataTable$Conditions)

      # Initialisation
      ModelDataTable <- data.frame()
      # y axis points are calculated. (limits 0.01% -- 99.99%) Necessary to have nice confidence bands.
      Proba <- seq(qnorm(0.0001),qnorm(0.9999),0.05)

      for (i in seq_along(ListConditions)){
          # Experimental conditions:
          Condition <- ListConditions[i]
          I <- CleanDataTable$Stress[CleanDataTable$Conditions==Condition][1]
          Temp <- CleanDataTable$Temperature[CleanDataTable$Conditions==Condition][1]  # °C

          # TTF calculation
          TTF <- exp(A)*(I*0.001/S)^(-n)*exp((Ea*e)/(k*(273.15+Temp))+ Proba * Scale)

          # Dataframe creation
          ModelDataTable <- rbind(ModelDataTable, data.frame('TTF'=TTF,'Status'=1,'Probability'=Proba,'Conditions'=Condition,'Stress'=I,'Temperature'=Temp))
      }

      # Drawing of the residual plots
      plot(nlsResiduals(Model))
      # Display of fit results
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
}


ErrorEstimation <- function(ExpDataTable, ModelDataTable, ConfidenceValue=0.95)
# Generation of confidence intervals
{
    # list of conditions
    ListConditions <- levels(ExpDataTable$Conditions)
    # DataFrame initialisation
    ConfidenceDataTable <- data.frame()

    if (length(ListConditions) != 0){

          for (i in seq_along(ListConditions)){

              NbData <- length(ExpDataTable$TTF[ExpDataTable$Conditions == ListConditions[i]])
              if (NbData > 30) {
                  mZP_Value <- qnorm((1 - ConfidenceValue) / 2) # Normal case. Valid if sample size > 30.
              } else {
                  mZP_Value <- qt((1 - ConfidenceValue) / 2, df=(NbData -1) ) # t-test statistic for low sample size
              }
              CDF <- pnorm(ModelDataTable$Probability[ModelDataTable$Conditions == ListConditions[i]])
              sef <- sqrt(CDF * (1 - CDF)/NbData) # TO BE CHECKED
              LowerLimit <- qnorm(CDF - sef * mZP_Value)
              HigherLimit <- qnorm(CDF + sef * mZP_Value)

              ConfidenceDataTable <- rbind(ConfidenceDataTable, data.frame('TTF'=ModelDataTable$TTF[ModelDataTable$Conditions == ListConditions[i]],
                                                                            'LowerLimit'=LowerLimit,'HigherLimit'=HigherLimit,'Conditions'=ListConditions[i]))
        }
    }
    return(ConfidenceDataTable)
}


CreateGraph <- function(ExpDataTable, ModelDataTable, ConfidenceDataTable, Title="", Scale="Lognormal", ErrorBands=TRUE, Save=TRUE)
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

    } else { # Lognormal
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

    # Save as png & pdf
    if (Save == TRUE){
        if (Title != ""){
            ggsave(filename=paste(Title,"png",sep="."),dpi=300)
            #ggsave(filename=paste(Title,"pdf",sep="."))
        } else {
            ggsave(filename="Chart.png",dpi=300)
            #ggsave(filename="Chart.pdf")
        }
    }
}

#' Electromigration data analysis
#'
#' Extract Black's parameters from a set of electromigration experiments.
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
#' BlackAnalysis()
#' BlackAnalysis(ErrorBand=FALSE)
#' @author Emmanuel Chery, \email{emmanuel.chery@@ams.com}
#' @import ggplot2 MASS scales grid nlstools
#' @export
BlackAnalysis <- function(ErrorBand=TRUE, ConfidenceValue=0.95, Save=TRUE)
{
    #rm(list=ls())
    ListFiles <- list.files(pattern="*exportfile.txt")
    DeviceID <- strsplit(ListFiles[1],split="_")[[1]][2]
    # case 1, there are one or several files available
    if (length(ListFiles) != 0){
          # Import the file(s) and create the 3 dataframes
          DataTable <- ReadDataAce(ListFiles,Scale="Lognormal")
          ModelDataTable <- BlackModelization(DataTable, DeviceID)
          ErrorDataTable <- ErrorEstimation(DataTable, ModelDataTable, ConfidenceValue)

    } else { # case 2, there are no files available
          print("You need to create the export files first!")
    }
    CreateGraph(DataTable,ModelDataTable,ErrorDataTable,DeviceID,Scale="Lognormal",ErrorBand,Save)
    #return(DataTable)
}