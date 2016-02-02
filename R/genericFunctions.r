################################################################################
###                                                                          ###
###    INFORMATIONS                                                          ###
###    ---------------------------------                                     ###
###                                                                          ###
###       PACKAGE NAME        amsReliability                                 ###
###       MODULE NAME         genericFunctions.r                             ###
###       VERSION             0.9.1                                          ###
###                                                                          ###
###       AUTHOR              Emmanuel Chery                                 ###
###       MAIL                emmanuel.chery@ams.com                         ###
###       DATE                2016/01/13                                     ###
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
###       This module includes generic functions usable with every           ###
###    degradation mechanism.                                                ###
###                                                                          ###
###                                                                          ###
###    FUNCTIONS                                                             ###
###    ---------------------------------                                     ###
###                                                                          ###
###       CalculProbability         Standard deviation/weibit calculation    ###
###       Clean                     Remove TTF associated to bad devices     ###
###       CreateDataFrame           Place experimental data in a table       ###
###       CreateGraph               In charge of data representation         ###
###       ErrorEstimation           Calculation of confidence Intervals      ###
###       Ranking                   Calculation of fraction estimators       ###
###                                                                          ###
################################################################################


###### List of Constants  ######
k <- 1.38E-23 # Boltzmann
e <- 1.6E-19 # electron charge
################################


CalculLifeTime <- function(Model, Area, Stress, Temperature, Probability,  Law="BlackLaw")
# Calcul the lifetime of a device for a given condition (Temp/stress) at a given failure rate
# with a given model.
# For TDDB, area is transformed in m² whereas it is already given in m² for EM.
# Currently supports Black equation and a TDDB lifetime model.
{
    if (Law == "BlackLaw") {
        # Parameters Extraction
        A <- coef(Model)[1]
        n <- coef(Model)[2]
        Ea <-coef(Model)[3]
        Scale <- coef(Model)[4]

        TTF <- exp(A)*(Stress*0.001/Area)^(-n)*exp((Ea*e)/(k*(273.15+Temperature))+ Probability * Scale)

    } else if (Law == "TDDB"){
        # Parameters Extraction
        t0 <- coef(Model)[1]
        g <- coef(Model)[2]
        Ea <- coef(Model)[3]
        beta <- coef(Model)[4]

        TTF <- exp(t0)*exp(-g*Stress)*exp((Ea*e)/(k*(Temperature+273.15)))*(Area*1E-12)^(-1/beta)*exp(Probability/beta)

    }
    return(TTF)
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
    # Remove ghost levels: levels were no samples are listed anymore.
    CleanedTable <- droplevels(CleanedTable)
    return(CleanedTable)
}


CreateDataFrame <- function(TTF, Status, Condition, Stress, Temperature, Scale="Lognormal", Dimension = 1)
# Creation of the dataframe assembling the TTF, the status of the samples,
# the probability, the condition (stickers for charts),
# the stress condition and the temperature used durng the stress.
# The probability is calculated according to Lognormal or Weibull distribution.
# Data are given clean.
# Data(TTF,Status,Probability,Conditions,Stress,Temperature, Dimension)
{
    rk <- Ranking(TTF) # Fraction estimator calculation
    if (Scale=="Weibull") {
        Proba <- CalculProbability(rk,Scale="Weibull") # Probability calculation Weibull
    } else {
        Proba <- CalculProbability(rk,Scale="Lognormal") # Probability calculation Lognormal
    }
    # Generation of the final data frame
    DataTable <- data.frame('TTF'=TTF,'Status'=Status,'Probability'=Proba,'Conditions'=Condition, 'Stress'=Stress, 'Temperature'=Temperature, 'Dimension'=Dimension)
    return(DataTable)
}


CreateGraph <- function(ExpDataTable, ModelDataTable, ConfidenceDataTable, Title="", Scale="Lognormal", ErrorBands=TRUE, Save=TRUE)
# Use the table prepared with CreateDataFrame and create the probability plot.
# Default is Lonormale scale but Weibull is available as an option.
{
    # x scale limits calculation based on the data.
    lim <- range(ExpDataTable$TTF[ExpDataTable$Status==1]) # Min of the values is stored in [1] and max in  [2]
    lim.high <- 10^(ceiling(log(lim[2],10)))
    lim.low <- 10^(floor(log(lim[1],10)))
    # 3 decades minimum are needed for a good looking chart.
    # In case the distribution is only in 1 decade, we add a decade at both ends
    if ((log10(lim.high) - log10(lim.low)) == 1 ) {
        lim.high <- lim.high * 10
        lim.low <- lim.low / 10
    # if we have already tw0 decades, we add one decade in the area where the data are closer to the edge
    } else if ((log10(lim.high) - log10(lim.low)) == 2) {
        if ((log10(lim[1]) - log10(lim.low)) < (log10(lim.high) - log10(lim[2]))){
            lim.low <- lim.low / 10
        } else {
            lim.high <- lim.high * 10
        }
    }

    # Now that we have the limits, we create the graph labels for x axis.
    GraphLabels <- 10^(seq(log10(lim.low),log10(lim.high)))
    # Now we create the minor ticks
    ind.lim.high <- log10(lim.high)
    ind.lim.low <- log10(lim.low)
    MinorTicks <- rep(seq(1,9), ind.lim.high - ind.lim.low ) * rep(10^seq(ind.lim.low, ind.lim.high-1), each=9)

    # Function used to calculate the distance between ticks for logscale. See line 166:
    # minor_breaks=trans_breaks(faceplant1, faceplant2, n=length(MinorTicks)))
    faceplant1 <- function(x) {
        return (c(x[1]*10^.25, x[2]/10^.25))
    }

    faceplant2 <- function(x) {
        return (MinorTicks)
    }
    #############################

    # Label for y axis
    # Dynamique labels as a function of the minimal probability observed.
    # Minimal proba is 0.01 %

    # Case 0: Proba min is above 1%
    if (Scale == "Weibull"){ # Weibull requires 63% and details in low %
        ListeProba <- c(1,2,3,5,10,20,30,40,50,63,70,80,90,95,99)
    } else { # Lognormal scale is symetric.
        ListeProba <- c(1,5,10,20,30,40,50,60,70,80,90,95,99)
    }

    MinProba <- min(ExpDataTable$Probability)

    if (MinProba <= CalculProbability(1/100,Scale)){ # Case 1: lower than 1%
        ListeProba <- c(0.1,ListeProba, 99.9)
    }
    if (MinProba <= CalculProbability(0.1/100,Scale)){ # Case 2: lower than 0.1%
        ListeProba <- c(0.01,ListeProba, 99.99)
    }

    # Probability vector used to draw y axis.
    ProbaNorm <- CalculProbability(ListeProba/100,Scale)



    # We are only going to plot samples where status is '1' (experiment is finished).
    # Table is sorted & conditions stay togeteher.
    CleanExpTable <- ExpDataTable[ExpDataTable$Status==1,]
    CleanExpTable <- CleanExpTable[order(CleanExpTable$"Conditions"),]

    # Graph creation with CleanTable
    Graph <- ggplot(data=CleanExpTable, aes(x=TTF, y=Probability, colour=Conditions, shape=Conditions))
    # box around chart + background
    Graph <- Graph + theme_linedraw() + theme(panel.background = element_rect(fill="gray90", color="black"))
    # Definition of scales
    Graph <- Graph + scale_x_log10(limits = c(lim.low,lim.high),breaks = GraphLabels,labels = trans_format("log10", math_format(10^.x)), minor_breaks=trans_breaks(faceplant1, faceplant2, n=length(MinorTicks)))
    Graph <- Graph + scale_y_continuous(limits=range(ProbaNorm), breaks=ProbaNorm, labels=ListeProba)
    # Grid definitions
    Graph <- Graph + theme(panel.grid.major = element_line(colour="white", size=0.25, linetype=1))
    Graph <- Graph + theme(panel.grid.minor = element_line(linetype=2, colour="white", size = 0.25))
    Graph <- Graph + theme(panel.grid.minor.y = element_line(linetype=0, colour="white", size = 0.25))
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
    Graph <- Graph + theme(axis.text.x = element_text(face="bold", size=16, margin=margin(0.4,0,0,0, "cm")))
    Graph <- Graph + theme(axis.text.y = element_text(size=16, margin=margin(0,0.4,0,0.2, "cm")))
    Graph <- Graph + theme(axis.ticks.length = unit(-0.25, "cm"))#, axis.ticks.margin = unit(0.4, "cm")) #Depreciated see margin above.
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


CreateModelDataTable <- function(Model, ListConditions, Area, Law="BlackLaw", Scale="Lognormal")
# Return a theoretical lifetime for a given set of Conditions and a given model.
# Result is returned in a dataframe.
# Currently supported model are Black equation and a TDDB lifetime model.
{
    # Initialisation
    ModelDataTable <- data.frame()
    # y axis points are calculated. (limits 0.01% -- 99.99%) Necessary to have nice confidence bands.
    Proba <- seq(CalculProbability(0.0001, Scale), CalculProbability(0.9999, Scale), 0.05)

    # Extraction of temperature and stress conditions
    Temperature <- sapply(ListConditions,function(x){strsplit(x,split="[mAV]*/")[[1]][2]})
    Temperature <- as.numeric(sapply(Temperature,function(x){substr(x,1, nchar(x)-2)}))
    Stress <- as.numeric(sapply(ListConditions,function(x){strsplit(x,split="[mAV]*/")[[1]][1]}))

    for (i in seq_along(Temperature)){

        # TTF calculation
        TTF <- CalculLifeTime(Model, Area, Stress[i], Temperature[i], Proba, Law)
        # Dataframe creation
        ModelDataTable <- rbind(ModelDataTable, data.frame('TTF'=TTF,'Status'=1,'Probability'=Proba,'Conditions'=ListConditions[i],'Stress'=Stress[i],'Temperature'=Temperature[i], 'Area'=Area))

    }
    return(ModelDataTable)
}


ErrorEstimation <- function(ExpDataTable, ModelDataTable, ConfidenceValue=0.95, Scale="Lognormal")
# Generation of confidence intervals
# Based on Kaplan Meier estimator and Greenwood confidence intervals
{
    # list of conditions
    ListConditions <- levels(ExpDataTable$Conditions)
    # DataFrame initialisation
    ConfidenceDataTable <- data.frame()

    if (length(ListConditions) != 0){

          for (condition in ListConditions){

              NbData <- length(ExpDataTable$TTF[ExpDataTable$Conditions == condition & ExpDataTable$Status == 1])
              if (NbData > 30) {
                  mZP_Value <- qnorm((1 - ConfidenceValue) / 2) # Normal case. Valid if sample size > 30.
              } else {
                  mZP_Value <- qt((1 - ConfidenceValue) / 2, df=(NbData -1) ) # t-test statistic for low sample size
              }

              if (Scale == "Weibull"){
                  CDF <- 1-exp(-exp(ModelDataTable$Probability[ModelDataTable$Conditions == condition]))
                  sef <- sqrt(CDF * (1 - CDF)/NbData)
                  LowerLimit <- log(-log(1-(CDF - sef * mZP_Value)))
                  HigherLimit <- log(-log(1-(CDF + sef * mZP_Value)))

              } else {
                  CDF <- pnorm(ModelDataTable$Probability[ModelDataTable$Conditions == condition])
                  sef <- sqrt(CDF * (1 - CDF)/NbData)
                  LowerLimit <- qnorm(CDF - sef * mZP_Value)
                  HigherLimit <- qnorm(CDF + sef * mZP_Value)
              }

              ConfidenceDataTable <- rbind(ConfidenceDataTable, data.frame('TTF'=ModelDataTable$TTF[ModelDataTable$Conditions == condition],
                                                                            'LowerLimit'=LowerLimit,'HigherLimit'=HigherLimit,'Conditions'=condition))
        }
    }
    return(ConfidenceDataTable)
}


FitDistribution <- function(DataTable,Scale="Lognormal")
# Extract simple distribution parameters (MTTF, scale) and return
# a ModelDataTable to plot the theoretical distribution
# Use fitdistr function
{
    # For each condtion we estimate a theoretical distribution
    ListConditions <- levels(DataTable$Conditions)

    # Initialisation of ModelDataTable
    ModelDataTable <- data.frame()

    for (ModelCondition in ListConditions){

        # Stress and Temperature stickers
        ModelStress <- DataTable$Stress[DataTable$Conditions==ModelCondition][1]
        ModelTemperature <- DataTable$Temperature[DataTable$Conditions==ModelCondition][1]

        # x axis limits are calculated
        lim <- range(DataTable$TTF[DataTable$Conditions==ModelCondition])
        lim.high <- 10^(ceiling(log(lim[2],10)))
        lim.low <- 10^(floor(log(lim[1],10)))
        # Generation of a vector for the calculation of the model. 200pts/decades
        x <- 10^seq(log(lim.low,10),log(lim.high,10),0.005)


        # Model calculation with the experimental TTF
        if (Scale=="Weibull") { # Weibull
              fit <- fitdistr(DataTable$TTF[DataTable$Conditions==ModelCondition & DataTable$Status==1],"weibull")
              fitShape <- fit$estimate[1]  # Beta
              fitScale <- fit$estimate[2]  # Characteristic time (t_63%)
              y <- CalculProbability(pweibull(x, fitShape, fitScale),"Weibull")
              # Display of Model parameters
              print(paste("Condition ",ModelCondition, " Beta= ", fitShape, " t63%=", fitScale,sep=""))

        } else { # Lognormale
              fit <- fitdistr(DataTable$TTF[DataTable$Conditions==ModelCondition & DataTable$Status==1],"lognormal")
              fitScale <- fit$estimate[1]  # meanlog
              fitShape <- fit$estimate[2]  # sdlog
              y <- CalculProbability(plnorm(x, fitScale, fitShape),"Lognormale")
              # Display of Model parameters
              print(paste("Condition ",ModelCondition, " Shape= ", fitShape, " MTTF=", exp(fitScale),sep=""))
        }

        # ModelDataTable creation
        ModelDataTable <- rbind(ModelDataTable, data.frame('TTF'=x,'Status'=1,'Probability'=y,'Conditions'=ModelCondition,'Stress'=ModelStress,'Temperature'=ModelTemperature) )
    }
    return(ModelDataTable)
}


FitResultsDisplay <- function(Model, DataTable, DeviceID)
# Given a model and an experimental dataset
# Return the parameters of the model and the residual error
# Save the information in a fit.txt file
{
    CleanDataTable <- DataTable[DataTable$Status==1,]
    # Residual Sum of Squares
    RSS <- sum(resid(Model)^2)
    # Total Sum of Squares: TSS <- sum((TTF - mean(TTF))^2))
    TSS <- sum(sapply(split(CleanDataTable[,1],CleanDataTable$Conditions),function(x) sum((x-mean(x))^2)))
    Rsq <- 1-RSS/TSS # R-squared measure

    # Drawing of the residual plots
    plot(nlsResiduals(Model))
    # Display of fit results
    cat(DeviceID,"\n")
    print(summary(Model))
    cat(paste("Residual squared sum: ",RSS,sep=""))
    # Save in a file
    capture.output(summary(Model),file="fit.txt")
    cat("Residual Squared sum:\t",file="fit.txt",append=TRUE)
    cat(RSS,file="fit.txt",append=TRUE)
    cat("\n \n",file="fit.txt",append=TRUE)
    cat("Experimental Data:",file="fit.txt",append=TRUE)
    cat("\n",file="fit.txt",append=TRUE)
    capture.output(DataTable,file="fit.txt",append=TRUE)
}


ModelFit <- function(dataTable, Law="BlackLaw")
# Perform a least square fit on an experimental dataset
# Return a model containing the parameters and the residuals.
{
    if (Law == "BlackLaw") {
        # Black model / Log scale: use of log10 to avoid giving too much importance to data with a high TTF
        Model <- nls(log10(TTF) ~ log10(exp(A)*(Stress*1E-3/Area)^(-n)*exp((Ea*e)/(k*(Temperature+273.15))+Scale*Probability)), dataTable,
                start=list(A=30,n=1,Ea=0.7,Scale=0.3),control= list(maxiter = 50, tol = 1e-7))#, minFactor = 1E-5, printEval = FALSE, warnOnly = FALSE))#,trace = T)
    } else if (Law == "TDDB"){
        # TDDB model / Log scale: use of log10 to avoid giving too much importance to data with a high TTF
        Model <- nls(log10(TTF) ~ log10(exp(t0)*exp(-g*Stress)*exp((Ea*e)/(k*(Temperature+273.15)))*(Area*1E-12)^(-1/beta)*exp(Probability/beta)), dataTable,
                start=list(t0=30,g=1,Ea=0.2,beta=1),control= list(maxiter = 50, tol = 1e-6))#, minFactor = 1E-5, printEval = FALSE, warnOnly = FALSE))#,trace = T)
    }
    return(Model)
}


Ranking <- function(TTF)
# Fraction estimator calculation
# rk(i)=(i-0.3)/(n+0.4)
# TTF is a vector.
{
    # ties.method="random" handles identical TTFs and provide a unique ID
    rk <- (rank(TTF, ties.method="random")-0.3)/(length(TTF)+0.4)
}


SelectFiles <- function()
# Allow graphical selection of multiple files.
# Return them as a list.
{
    # Create the Path.
    # initFile is used to enter the right directory.
    # initFile is removed before returning the list
    initFile <- list.files(pattern="bidule.txt")
    path2Current <- paste(getwd(), "/", initFile[1], sep="")

    # Filters for file selection
    Filters <- matrix(c("All files", "*", "Export Files", "*exportfile.txt", "Text", ".txt"),3, 2, byrow = TRUE)

    # Gui for file selection
    selection <- tk_choose.files(default = path2Current, caption = "Select files",
                            multi = TRUE, filters = Filters, index = 1)

    # Cleaning to remove the path and keep only the filename (last item)
    if (Sys.info()[['sysname']] == "Windows"){
        # List of file
        listFiles <- sapply(strsplit(selection,split="/"),function(x){x[length(x)]})
        # new working path
        newWD <- substr(selection[1], 1, nchar(selection[1])-nchar(listFiles[1]))
        setwd(newWD)
    } else {
        listFiles <- sapply(strsplit(selection[-1],split="/"),function(x){x[length(x)]})
    }

    return(listFiles)
}


SortConditions <- function(ListConditions)
# Sort a list of conditions to avoid 6mA being
# bigger as 14mA
# Return a list of Conditions sorted.
{
  Temperature <- sapply(ListConditions,function(x){strsplit(x,split="[mAV]*/")[[1]][2]})
  Temperature <- as.numeric(sapply(Temperature,function(x){substr(x,1, nchar(x)-2)}))
  Currents <- as.numeric(sapply(ListConditions,function(x){strsplit(x,split="[mAV]*/")[[1]][1]}))
  Table <- data.frame("Conditions"=ListConditions,"Current"=Currents,"Temperature"=Temperature)
  Table <-  Table[order(Table$Temperature),]
  ListCurrents <- levels(factor(Currents))

  SortedTable <- data.frame()
  for (current in ListCurrents){
    SortedTable <-  rbind(SortedTable,Table[Table$Current==current,])
  }
  return(as.character(SortedTable$Conditions))
}


OrderConditions <- function(DataTable)
# Order a list of conditions to avoid 6mA being
# bigger as 14mA.
# Return a vector of indice.
{
    SortedListConditions <- SortConditions(levels(DataTable$Conditions))
    VecIndices <- c()
    for (condition in SortedListConditions){
        VecIndices <- c(VecIndices, which(DataTable$Conditions == condition))
    }
    return(VecIndices)
}
