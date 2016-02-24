################################################################################
###                                                                          ###
###    INFORMATIONS                                                          ###
###    ---------------------------------                                     ###
###                                                                          ###
###       PACKAGE NAME        amsReliability                                 ###
###       MODULE NAME         GUI.r                                          ###
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
###       This module is dedicated to the GUI function.                      ###
###    It allows starting all other functions by using a graphical menu.     ###
###                                                                          ###
###                                                                          ###
###    FUNCTIONS                                                             ###
###    ---------------------------------                                     ###
###                                                                          ###
###       amsReliability          Display a menu to select an analysis       ###
###                                                                          ###
################################################################################


#' Reliability data analysis
#'
#' This function displays a menu listing all available analyses.
#' User is guided with a GUI in order to select the analysis.
#'
#' @param None
#'
#' @return None
#'
#' @examples
#' RelAnalysis()
#' @author Emmanuel Chery, \email{emmanuel.chery@@ams.com}
#' @import ggplot2 MASS scales nlstools tcltk
#' @export
RelAnalysis <- function()
{
    # Main windows
    tt <- tktoplevel()
    tkwm.minsize(tt, 500, 300)
    tkwm.resizable(tt, TRUE, TRUE)
    tkwm.title(tt, paste("amsReliability ",packageVersion("amsReliability")))

    done <- tclVar(0)

    tkconfigure(tt, borderwidth= 20, width="70")

    # Main title
    titleFont <- tkfont.create(size=10,weight ="bold")
    titleWindow <- tklabel(tt, text= "amsReliability")
    tkconfigure(titleWindow, font=titleFont)

    # Menu
    topMenu <- tkmenu(tt)
    tkconfigure(tt, menu=topMenu)
    # File
    fileMenu <- tkmenu(topMenu, tearoff=FALSE)
    tkadd(fileMenu,"command",label="Quit",command=function() tkdestroy(tt))
    tkadd(topMenu,"cascade",label="File",menu=fileMenu)
    # About & Update Menu
    about <- tkmenu(topMenu, tearoff=FALSE)
    tkadd(about,"command",label="Update",
          command=function() {AutoUpdate(); tkdestroy(tt); RelAnalysis()}) # reload GUI after update.

    boxFont <- tkfont.create(size=10,weight ="bold")
    boxText <- tklabel(about, text="Welcome in amsReliability!")
    tkadd(about,"command",label="About",
          command=function()
                tkmessageBox(title="About",
                            message=paste("\n            amsReliability ", packageVersion("amsReliability"),"\n\n   \"Reliability Data Analysis Tool\"\n\n\n\nEmmanuel Chery\t(2015--2016)",sep=""),
                            icon="info"))
    tkadd(topMenu,"cascade",label="Help",menu=about)

    tkfocus(tt)


    # Main windows:
    ## Text
    titlePolice <- tkfont.create(size=14, weight="bold")
    subtitlePolice <- tkfont.create(size=10, weight="bold")
    welcomeText <- tklabel(tt, text="Welcome in amsReliability!")
    versionText <- tklabel(tt, text=paste("\n\t\t\t\tVersion ",packageVersion("amsReliability")," is currently loaded." ))
    tkconfigure(welcomeText,font=titlePolice)
    selectText <- tklabel(tt, text="\n\nPlease select the analysis you want to perform:\n")
    tkconfigure(selectText,font=subtitlePolice)
    # tkpack(welcomeText, versionText, selectText)
    tkgrid(welcomeText, columnspan=12)
    tkgrid(versionText, columnspan=12)
    tkgrid(selectText, columnspan=12)

    ## List of possibilities
    listAnalyses <- c("TDDB modelization and model parameter extraction","EM modelization and model parameter extraction", "EM exportfiles creation")
    maxChar <- max(sapply(listAnalyses, nchar))
    #myList<-tklistbox(tt, height= length(listAnalyses)+1, width= maxChar, selectmode="single",background="white")

    #tkgrid(myList,columnspan=12)
    #for (analysis in listAnalyses){
    #    tkinsert(myList,"end",analysis)
    #}
    # default selection
    #tkselection.set(myList,0)

    initMess <- tclVar("")
    myList <- ttkcombobox(tt, values=listAnalyses, textvariable=initMess, state="readonly", width=maxChar)
    tkgrid(myList, columnspan=18)

    ## Options
    ### text
    optionPolice <- tkfont.create(size=8, weight="bold")
    optionText <- tklabel(tt, text="\nOptions:")
    tkconfigure(optionText,font=optionPolice)
    tkgrid(optionText,row=6, column=1)
    ### Buttons
    error <- tkcheckbutton(tt)
    errorValue <- tclVar("0")
    tkconfigure(error,variable=errorValue)
    save <- tkcheckbutton(tt)
    saveValue <- tclVar("1")
    tkconfigure(save,variable=saveValue)
    tkgrid(tklabel(tt,text="\tDisplay error bands:"),error, columnspan=2)
    tkgrid(tklabel(tt,text="\tSave chart:"),save, columnspan=2)
    # Confidence level
    clevel <- tclVar("0.95")
    clevelEntry <-tkentry(tt,width="5",textvariable=clevel)
    tkgrid(tklabel(tt,text="\tConfidence level:"),clevelEntry,columnspan=2)

    # OnOK function
    OnOK <- function()
    {
        # Option values:
        errorband <- as.character(tclvalue(errorValue))
        save <- as.character(tclvalue(saveValue))
        if (errorband == "1"){
            errorband <- TRUE
        } else {
            errorband <- FALSE
        }

        if (save == "1"){
            save <- TRUE
        } else {
            save <- FALSE
        }
        # Confidence level
        confidence <- as.numeric(tclvalue(clevel))

        # List value:
        # userChoice <- as.numeric(tkcurselection(myList))+1 # Index starts at 1 now
        # ComboBox value:
        userChoice <- tclvalue(initMess)

        tkdestroy(tt)

        # Sanity check
        if (userChoice == ""){
            tkmessageBox(title = "Selection empty", message = "No analyses were selected!", icon = "error", type = "ok")
        }
        if (userChoice == listAnalyses[1]){
            OxideTDDB(ErrorBand = errorband, ConfidenceValue = confidence, Save = save)
        } else if (userChoice == listAnalyses[2]){
            BlackAnalysis(ErrorBand = errorband, ConfidenceValue = confidence, Save = save)
        } else if (userChoice == listAnalyses[3]){
            CreateExportFiles()
        }

        Sys.sleep(0.5)
        RelAnalysis()
    }


    ## Button OK & Quit
    OK.but <- tkbutton(tt,text="OK",command=OnOK)
    Quit.but <- tkbutton(tt,text="Quit",command=function() tkdestroy(tt))
    tkgrid(OK.but,row=9, column=11)
    tkgrid(Quit.but,row=9, column=12)

    # Watch for windows closure.
    tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)
    tkwait.variable(done)

    if(tclvalue(done)=="2") invisible()

    tkdestroy(tt)
}


AutoUpdate <- function()
{
    print("AutoUpdate ongoing")
    library(devtools)
    install_github("KDB2/amsReliability")
    detach("package:amsReliability", unload=TRUE)
    library("amsReliability")
    return()
}
