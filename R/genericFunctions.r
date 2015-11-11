################################################################################
###                                                                          ###
###    INFORMATIONS                                                          ###
###    ---------------------------------                                     ###
###                                                                          ###
###       PACKAGE NAME        amsReliability                                 ###
###       MODULE NAME         genericFunctions.r                             ###
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


FitDistribution <- function(DataTable,Scale="Lognormal")
# Extract simple distribution parameters (MTTF, scale) and return
# a ModelDataTable to plot the theoretical distribution
{
    # For each condtion we estimate a theoretical distribution
    ListConditions <- levels(DataTable$Conditions)

    # Initialisation of ModelDataTable
    ModelDataTable <- data.frame()

    for (i in seq_along(ListConditions)){

        # Condition, Stress and Temperature stickers
        ModelCondition <- ListConditions[i]
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

Ranking <- function(TTF)
# Fraction estimator calculation
# rk(i)=(i-0.3)/(n+0.4)
# TTF is a vector.
{
    # ties.method="random" handles identical TTFs and provide a unique ID
    rk <- (rank(TTF, ties.method="random")-0.3)/(length(TTF)+0.4)
}
