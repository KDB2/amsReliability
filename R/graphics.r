################################################################################
###                                                                          ###
###    INFORMATIONS                                                          ###
###    ---------------------------------                                     ###
###                                                                          ###
###       PACKAGE NAME        amsReliability                                 ###
###       MODULE NAME         graphics.r                                     ###
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
###       This module includes the graphics functions.                       ###
###                                                                          ###
###                                                                          ###
###    FUNCTIONS                                                             ###
###    ---------------------------------                                     ###
###                                                                          ###
###       CreateGraph               In charge of data representation         ###
###                                                                          ###
###                                                                          ###
################################################################################


#### Main graphics function ####

CreateGraph <- function(dataExp, ModelDataTable = NULL, ConfidenceDataTable = NULL,
                aesVec = c("TTF", "Probability", "Conditions"), title="", axisTitles = c("",""),
                scale.x = "Log", scale.y="Lognormal", save=TRUE)
# Use the table prepared with CreateDataFrame and create the probability plot.
# Default y scale is Lonormale scale but Weibull, Log and Lin (degradation charts (%)) are available as an option.
# Default x scale is log but linear (lin) is available in option.
# aes are given in vector form with max 4 names: aesVec = c(abs, ord, Conditions, Split)
{

    # Table is sorted & conditions stay togeteher.
    dataExp <- dataExp[order(dataExp$"Conditions"),]

    # Graph creation with dataExp
    if (length(aesVec) == 3){
        # Chart
        Graph <- ggplot(data=dataExp, aes(x=dataExp[[aesVec[1]]], y=dataExp[[aesVec[2]]],
                        colour=dataExp[[aesVec[3]]], shape=dataExp[[aesVec[3]]]))
        # Legend title
        Graph <- Graph + labs(colour = aesVec[3]) + labs(shape=aesVec[3])
    } else if (length(aesVec) == 4){
        # Chart
        Graph <- ggplot(data=dataExp, aes(x=dataExp[[aesVec[1]]], y=dataExp[[aesVec[2]]],
                colour=dataExp[[aesVec[3]]], shape=dataExp[[aesVec[4]]]))
        # Legend title
        Graph <- Graph + labs(colour = aesVec[3]) + labs(shape=aesVec[4])
    }

    # Definition of scales
    Graph <- CreateAxis.x(Graph, dataExp[[aesVec[1]]], scale.x)
    Graph <- CreateAxis.y(Graph, dataExp[[aesVec[2]]], scale.y)
    # Add default options
    Graph <- GraphBase(Graph, title)
    # Font size
    Graph <- AddThemeAxis(Graph, scale.x, scale.y)
    # x/y titles
    Graph <- Graph + xlab(axisTitles[1]) + ylab(axisTitles[2])
    # Grid
    if (scale.y == "Log" | scale.y == "Lin") {
        Graph <- AddGrid(Graph, minor.y = TRUE)
    } else { # Probability plots don't need a minor grid on y.
        Graph <- AddGrid(Graph, minor.y = FALSE)
    }

    # Add the theoretical model
    if (!is.null(ModelDataTable)){
        Graph <- Graph + geom_line(data=ModelDataTable, aes(x=ModelDataTable[[aesVec[1]]],
                            y=ModelDataTable[[aesVec[2]]], color=ModelDataTable[[aesVec[3]]], shape=NULL),
                            size=0.8)
    }
    # Add the confidence intervals
    if (!is.null(ConfidenceDataTable)){
        Graph <- Graph + geom_line(data=ConfidenceDataTable,
                        aes(x=ConfidenceDataTable[[aesVec[1]]], y=LowerLimit, color=ConfidenceDataTable[[aesVec[3]]],
                            shape=NULL), linetype="dashed", size=0.8)
        Graph <- Graph + geom_line(data=ConfidenceDataTable,
                        aes(x=ConfidenceDataTable[[aesVec[1]]], y=HigherLimit, color=ConfidenceDataTable[[aesVec[3]]],
                            shape=NULL), linetype="dashed",size=0.8)
    }

    print(Graph)

    # Save as png or pdf
    if (save == TRUE){
        GraphSave(title, extension ="png")
    }
}



#### Graphics utilities #####

GraphBase <- function(graph, title)
# Add default parameters to a graph
# background, legend, shape and color of points
{
    # box around chart + background
    graph <- graph + theme_linedraw() + theme(panel.background = element_rect(fill="gray90", color="black"))
    # Controled symbol list -- Max is 20 conditions on the chart.
    graph <- graph + scale_shape_manual(values=c(19,15,17,16,19,15,17,16,19,15,17,16,19,15,17,16,19,15,17,16))
    graph <- graph + scale_colour_manual(values = c("#d53e4f","#3288bd","#66a61e","#f46d43","#e6ab02","#8073ac","#a6761d","#666666","#bc80bd","#d53e4f","#3288bd","#66a61e","#f46d43","#e6ab02","#8073ac","#a6761d","#666666","#bc80bd","#d53e4f","#3288bd")) # "#5e4fa2" ,"#66c2a5", "#fec44f",
    # legend size
    graph <- graph + theme(legend.title = element_text(size=14, face="bold"))
    graph <- graph + theme(legend.text = element_text(size = 12))
    # Box around legend
    graph <- graph + theme(legend.background = element_rect())
    graph <- graph + theme(legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"))
    #Box around the conditions in legend
    graph <- graph + theme(legend.key = element_rect(fill="gray90", colour = "black", linetype=0))
    # Add a title
    graph <- graph + ggtitle(title)
    graph <- graph + theme(plot.title = element_text(face="bold", size=18, hjust=0.5))
    # Font size
    graph <- graph + theme(axis.title.x = element_text(face="bold", size=16))
    graph <- graph + theme(axis.title.y = element_text(face="bold", size=16))
    # Size of symbols
    graph <- graph + geom_point(size=4)

    return(graph)
}


AddGrid <- function(graph, minor.y = TRUE)
# Add Grid to graph
{
    # Grid definitions
    graph <- graph + theme(panel.grid.major = element_line(colour="white", size=0.25, linetype=1))
    graph <- graph + theme(panel.grid.minor = element_line(linetype=2, colour="white", size = 0.25))
    if (minor.y){
        graph <- graph + theme(panel.grid.minor.y = element_line(linetype=2, colour="white", size = 0.25))
    } else {
        graph <- graph + theme(panel.grid.minor.y = element_line(linetype=0, colour="white", size = 0.25))
    }
    return(graph)
}


AddThemeAxis <- function(graph, scale.x = "Log", scale.y="Lognormal" )
# Label/ticks size
# if the scale is logarithmic, font is in bold. (bug with 10^x presentation)
{
    if (scale.x == "Log"){
        graph <- graph + theme(axis.text.x = element_text(face="bold", size=16, margin=margin(0.4,0,0,0, "cm")))
    } else {
        graph <- graph + theme(axis.text.x = element_text(size=16, margin=margin(0.4,0,0,0, "cm")))
    }

    if (scale.y == "Log"){
        graph <- graph + theme(axis.text.y = element_text(face="bold", size=16, margin=margin(0,0.4,0,0.2, "cm")))
    } else {
        graph <- graph + theme(axis.text.y = element_text(size=16, margin=margin(0,0.4,0,0.2, "cm")))
    }

    graph <- graph + theme(axis.ticks.length = unit(-0.25, "cm"))

    return(graph)
}


ExtractLimits <- function(data, minDecades=3)
# Return the limits of the Data
#  Add necessary decades to follow the minimal number of decades requested.
{
    if (min(data) <= 0){ # if some values are negative, we simply return the range
        return(range(data))
    } else {
        lim <- range(data) # Min of the values is stored in [1] and max in  [2]
        lim.high <- 10^(ceiling(log(lim[2],10)))
        lim.low <- 10^(floor(log(lim[1],10)))

        nbDecades <- log10(lim.high) - log10(lim.low)
        diff <- minDecades - nbDecades

        # In case the number of minimal decades is not reached, we add decades at both ends
        if ( nbDecades < minDecades ) {
            # We estimate where the data are closer to the edge in order to add the aditional decade on this side.
            # (Odd case)
            if ((log10(lim[1]) - log10(lim.low)) < (log10(lim.high) - log10(lim[2]))){ # additional decade on the left.
                lim.low <- lim.low / 10^(diff %/% 2 + diff %% 2)
                lim.high <- lim.high * 10^(diff %/% 2)
            } else { # aditional decade on the right
                lim.low <- lim.low / 10^(diff %/% 2 )
                lim.high <- lim.high * 10^(diff %/% 2 + diff %% 2)
            }
        }
        return(c(lim.low, lim.high))
    }
}


GraphSave <- function(title, extension="png")
# Save a graph
{
    if (title != ""){
        ggsave(filename=paste(title, extension ,sep="."),dpi=300)
    } else {
        ggsave(filename=paste("Chart", extension, sep="."),dpi=300)
    }
}


GraphTargetLines <- function(graph, x=NULL, y=NULL, colorx="red", typex = 2, colory="red", typey = 2)
# Add target lines to a graph
{
    if (!is.null(x) & is.numeric(x)){
        graph <- graph + geom_vline(xintercept = x, color = colorx, linetype = typex)
    }

    if (!is.null(y) & is.numeric(x)){
        graph <- graph + geom_hline(xintercept = y, color = colory, linetype = typey)
    }

    return(graph)
}


CreateAxisLog <- function(graph, scaleLimits, axis = "x")
# Create a log axis.
# Limits <- c(lim.low, lim.high)
# axis: x or y
{
    lim.low <- scaleLimits[1]
    lim.high <- scaleLimits[2]
    graphLabels <- 10^(seq(log10(lim.low),log10(lim.high)))

    # Now we create the minor ticks
    minorTicks <- rep(seq(1,9), log10(lim.high) - log10(lim.low) ) * rep(10^seq(log10(lim.low), log10(lim.high)-1), each=9)

    # Function used to calculate the distance between ticks for logscale. See line 166:
    # minor_breaks=trans_breaks(faceplant1, faceplant2, n=length(MinorTicks)))
    faceplant1 <- function(x) {
        return (c(x[1]*10^.25, x[2]/10^.25))
    }

    faceplant2 <- function(x) {
        return (minorTicks)
    }

    if (axis == "x"){
        graph <- graph + scale_x_log10(limits = c(lim.low,lim.high), breaks = graphLabels, labels = trans_format("log10", math_format(10^.x)), minor_breaks=trans_breaks(faceplant1, faceplant2, n=length(minorTicks)))
        graph <- graph + annotation_logticks(sides='tb')
    } else if (axis == "y"){
        graph <- graph + scale_y_log10(limits = c(lim.low,lim.high), breaks = graphLabels, labels = trans_format("log10", math_format(10^.x)), minor_breaks=trans_breaks(faceplant1, faceplant2, n=length(minorTicks)))
        graph <- graph + annotation_logticks(sides='lr')
    }
    return(graph)
}


CreateAxisLin <- function(graph, scaleLimits, axis = "x")
# Create a lin axis.
# Limits <- c(lim.low, lim.high) # Not used in the current implementation
# axis: x or y
{
    lim.low <- scaleLimits[1]
    lim.high <- scaleLimits[2]

    if (axis == "x"){
        graph <- graph + scale_x_continuous(limits = NULL, breaks = waiver(), labels = waiver(), minor_breaks= waiver())

    } else if (axis == "y"){
        graph <- graph + scale_y_continuous(limits = NULL, breaks = waiver(), labels = waiver(), minor_breaks= waiver())
    }
    return(graph)
}


CreateAxisWeibull <- function(graph, minProba, axis = "y")
# Create a Weibull scale on axis y
# minProba defines the minimal probability being displayed.
# minProba is given in weibit
{
    minProba.ind <- min(0, floor(log10( (1-exp(-exp( minProba)))  *100)))
    listeProba <- c( 10^seq(minProba.ind,0),2,3,5,10,20,30,40,50,63,70,80,90,95, (100 - 10^seq(0 , minProba.ind)) )
    probaNorm <- CalculProbability(listeProba/100,"Weibull")

    if (axis == "x"){
        graph <- graph + scale_x_continuous(limits=range(probaNorm), breaks=probaNorm, labels=listeProba)
    } else if (axis == "y"){
        graph <- graph + scale_y_continuous(limits=range(probaNorm), breaks=probaNorm, labels=listeProba)
    }
    return(graph)
}


CreateAxisLognormal <- function(graph, minProba, axis = "y")
# Create a lognormal scale on axis y
# minProba defines the minimal probability being displayed.
# minProba is given in standard deviations.
{
    minProba.ind <- min(0, floor(log10(pnorm(minProba)*100)))
    listeProba <- c( 10^seq(minProba.ind,0),5,10,20,30,40,50,60,70,80,90,95, (100 - 10^seq(0 , minProba.ind)) )
    probaNorm <- CalculProbability(listeProba/100,"Lognormal")

    if (axis == "x"){
        graph <- graph + scale_x_continuous(limits=range(probaNorm), breaks=probaNorm, labels=listeProba)
    } else if (axis == "y"){
        graph <- graph + scale_y_continuous(limits=range(probaNorm), breaks=probaNorm, labels=listeProba)
    }

    return(graph)
}


CreateAxis.x <- function(graph, data, scale)
# Generic function to create an x axis.
# Scale can be Log (Logarithmic) or Lin (Linear).
{
    if (scale == "Log"){
        lim <- ExtractLimits(data, minDecades = 3)
        graph <- CreateAxisLog(graph, scaleLimits = lim ,axis="x")
    } else if (scale == "Lin"){
        lim <- ExtractLimits(data, minDecades = 1)
        graph <- CreateAxisLin(graph, scaleLimits = lim ,axis="x")
    }
    return(graph)
}


CreateAxis.y <- function(graph, data, scale)
# Generic function to create an y axis.
# Scale can be Log (Logarithmic) or Lin (Linear)
# or Weibull or Lognormal.
{
    if (scale == "Log"){
        lim <- ExtractLimits(data, minDecades = 3)
        graph <- CreateAxisLog(graph, scaleLimits = lim ,axis="y")
    } else if (scale == "Lin"){
        lim <- ExtractLimits(data, minDecades = 1)
        graph <- CreateAxisLin(graph, scaleLimits = lim ,axis="y")
    } else if (scale == "Weibull"){
        lim = min(data)
        graph <- CreateAxisWeibull(graph, minProba = lim ,axis="y")
    } else if (scale == "Lognormal"){
        lim = min(data)
        graph <- CreateAxisLognormal(graph, minProba = lim ,axis="y")
    }
    return(graph)
}
