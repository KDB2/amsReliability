\name{amsReliability-package}
\alias{amsReliability-package}
\alias{amsReliability}
\docType{package}
\title{
\packageTitle{amsReliability}
}
\description{
\packageDescription{amsReliability}
}
\details{
\bold{CreateExportFiles.EM(ForceOverWrite=FALSE)}:

Use degradation and TCR files from Qualitau ACE and MIRA equipments to create respective exportfiles.  Times to failure are converted in seconds. If several conditions are present in one ACE experiment file, one exportfiles per condition will be created.
If exportfiles are already present, they are not replaced except if the ForceOverWrite option is set to TRUE.


\bold{BlackAnalysis(ErrorBand=TRUE,ConfidenceValue=0.95,Save=TRUE)}:

Using the exportfiles previously generated, perform the analysis of electromigration data. Extraction of Black parameters is performed. A chart is generated with the experimental data and the calculated model.
Confidence intervals can be disabled and the confidence value used can be choosen. If save is set to true, the chart is saved as PNG.
}
\author{
\packageAuthor{amsReliability}

Maintainer: \packageMaintainer{amsReliability}
}
\references{

}

\keyword{ package }
\seealso{
}
\examples{

}
