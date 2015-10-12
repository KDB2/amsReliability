# Script collection for ams AG process reliability team.
# Allow to display and process electromigration data.
# Extraction of Black's parameters is performed.
# September 2015
# Emmanuel Chery
# Version 0.3




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


StackData <- function(DataTable1, DataTable2)
# Merge 2 DataTable
{
    NewDataTable <- merge(DataTable1, DataTable2, all=TRUE)
    NewDataTable <- NewDataTable[order(NewDataTable$"Conditions"),]
    return(NewDataTable)
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
    ListDevice <- read.delim("//fsup04/fntquap/Common/Qual/Process_Reliability/Process/0.18_um_Technology/0.18 FabB/BEOL/Elmig/ListDeviceName.txt")
    W <- ListDevice$Width[ListDevice$Device==DeviceID] # micrometers
    H <- ListDevice$Width[ListDevice$Device==DeviceID] # micrometers
    S <- W*H*1E-12 # m^2

    # Physical constants
    k <- 1.38E-23 # Boltzmann
    e <- 1.6E-19 # electron charge

    # Remove the units where status is 0
    DataTable <- DataTable[DataTable$Status==1,]

    # Black model / Log scale: use of log10 to avoid giving too much importance to data with a high TTF
    nls.control(maxiter = 100, tol = 1e-15, minFactor = 1/1024, printEval = FALSE, warnOnly = FALSE)
    Model <- nls(log10(TTF) ~ log10(exp(A)*(Stress*1E-3/S)^(-n)*exp((Ea*e)/(k*(Temperature+273.15))+Scale*Probability)), DataTable, start=list(A=30,n=1,Ea=0.7,Scale=0.3))#,trace = T)
    #Model <- nls(TTF ~ exp(A)*(Stress*1E-3/S)^(-n)*exp((Ea*e)/(k*(Temperature+273.15))+Scale*Probability), DataTable, start=list(A=30,n=1,Ea=0.7,Scale=0.3))
    # Parameters Extraction
    A <- coef(Model)[1]
    n <- coef(Model)[2]
    Ea <-coef(Model)[3]
    Scale <- coef(Model)[4]
    # Residual Sum of Squares
    RSS <- sum(resid(Model)^2)
    # Total Sum of Squares: TSS <- sum((TTF - mean(TTF))^2))
    TSS <- sum(sapply(split(DataTable[,1],DataTable$Conditions),function(x) sum((x-mean(x))^2)))
    Rsq <- 1-RSS/TSS # R-squared measure


    # Using the parameters and the conditions, theoretical distributions are created
    ListConditions <- levels(DataTable$Conditions)
    # Initialisation with first condition
    #####################################

    # Conditions
    Condition <- ListConditions[1]
    I <- DataTable$Stress[DataTable$Conditions==ListConditions[1]][1]
    Temp <- DataTable$Temperature[DataTable$Conditions==ListConditions[1]][1]  # °C

    # y axis points are calculated. (limits 0.01% -- 99.99%) Necessary to have nice confidence bands.
    Proba <- seq(qnorm(0.0001),qnorm(0.9999),0.05)
    # TTF calculation
    TTF <- exp(A)*(I*0.001/S)^(-n)*exp((Ea*e)/(k*(273.15+Temp))+ Proba * Scale)

    # Dataframe creation
    ModelDataTable <- data.frame('TTF'=TTF,'Status'=1,'Probability'=Proba,'Conditions'=Condition,'Stress'=I,'Temperature'=Temp)


    # Loop to create the DataFrame
    if ( length(ListConditions) > 1 ){
        for (i in 2:length(ListConditions)){

            # Conditions
            Condition <- ListConditions[i]
            I <- DataTable$Stress[DataTable$Conditions==ListConditions[i]][1]
            Temp <- DataTable$Temperature[DataTable$Conditions==ListConditions[i]][1]  # °C

            # TTF calculation
            TTF <- exp(A)*(I*0.001/S)^(-n)*exp((Ea*e)/(k*(273.15+Temp))+ Proba * Scale)

            # Dataframe creation
            NewData <- data.frame('TTF'=TTF,'Status'=1,'Probability'=Proba,'Conditions'=Condition,'Stress'=I,'Temperature'=Temp)

            #Stack in 1 global table
            ModelDataTable <- StackData(ModelDataTable,NewData)
        }
    }
    write.table(data.frame('A'=A,'n'=n,'Ea'=Ea,'Scale'=Scale,"RSS"=RSS,"Rsq=",Rsq),"fit.txt",quote=FALSE,sep="\t")
    print(paste("Ea=",Ea,"eV, n=",n,", A=",A," Scale=",Scale," RSS=",RSS," Rsq=",Rsq,sep=""))
    return(ModelDataTable)
}


ErrorEstimation <- function(ExpDataTable, ModelDataTable, ConfidenceValue=0.95)
# Genration of confidence intervals
{
    # list of conditions
    ListConditions <- levels(ExpDataTable$Conditions)

    if (length(ListConditions) != 0){
          # DataFrame initialisation
          NbData <- length(ExpDataTable$TTF[ExpDataTable$Conditions == ListConditions[1]])
          if (NbData > 30) {
              mZP_Value <- qnorm((1 - ConfidenceValue) / 2) # Normal case. Valid if sample size > 30.
          } else {
              mZP_Value <- qt((1 - ConfidenceValue) / 2, df=(NbData -1) ) # t-test statistic for low sample size
          }
          CDF <- pnorm(ModelDataTable$Probability[ModelDataTable$Conditions == ListConditions[1]])
          sef <- sqrt(CDF * (1 - CDF)/NbData) # TO BE CHECKED
          LowerLimit <- qnorm(CDF - sef * mZP_Value)
          HigherLimit <- qnorm(CDF + sef * mZP_Value)

          ConfidenceDataTable <- data.frame('TTF'=ModelDataTable$TTF[ModelDataTable$Conditions == ListConditions[1]],'LowerLimit'=LowerLimit,'HigherLimit'=HigherLimit,'Conditions'=ListConditions[1])

          if (length(ListConditions) > 1) {
              for (i in 2:length(ListConditions)){
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

                NewData <- data.frame('TTF'=ModelDataTable$TTF[ModelDataTable$Conditions == ListConditions[i]],'LowerLimit'=LowerLimit,'HigherLimit'=HigherLimit,'Conditions'=ListConditions[i])
                ConfidenceDataTable <- StackData(ConfidenceDataTable,NewData)
              }
          }
      }
    return(ConfidenceDataTable)
}


CreateGraph <- function(ExpDataTable, ModelDataTable, ConfidenceDataTable, Title="", Scale="Lognormale", ErrorBands=TRUE, Save=TRUE)
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
    if (Save == TRUE){
        if (Title != ""){
            ggsave(paste(Title,"png",sep="."))
        } else {
            ggsave("Chart.png")
        }
    }
}


BlackAnalysis <- function(Scale="Lognormale",ErrorBand=TRUE, ConfidenceValue=0.95, Save=TRUE)
# Main function calling the other. The one to use to open all the files.
# Open all the exportfiles from the workfolder
{
    #rm(list=ls())
    ListFiles <- list.files(pattern="*exportfile.txt")
    DeviceID <- strsplit(ListFiles[1],split="_")[[1]][2]
    # case 1, there are one or several files available
    if (length(ListFiles) != 0){
          # Import the first file to create the 3 dataframes
          DataTable <- ReadDataAce(ListFiles[1],Scale)

          # Let's now check if other files are available
          if (length(ListFiles) > 1){
                # loop to open all the files and stack them in the dataframe
                for (i in 2:length(ListFiles)){
                    NewDataTable <- ReadDataAce(ListFiles[i],Scale)
                    # Merging the tables
                    DataTable <- StackData(DataTable,NewDataTable)
                }
          }
          ModelDataTable <- BlackModelization(DataTable, DeviceID)
          ErrorDataTable <- ErrorEstimation(DataTable, ModelDataTable, ConfidenceValue)

    } else { # case 2, there are no files available
          print("You need to create the export files first!")
    }
    CreateGraph(DataTable,ModelDataTable,ErrorDataTable,DeviceID,Scale,ErrorBand,Save)
    #return(DataTable)
}
