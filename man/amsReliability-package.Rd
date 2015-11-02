\name{amsReliability-package}
\alias{amsReliability-package}
\alias{amsReliability}
\docType{package}
\title{amsReliability package}
\description{
\packageDescription{amsReliability}
}
\details{

2 functions are currently available within this package.

\bold{CreateExportFiles}: Automatically creates the exportfiles for a set of electromigration experiments
by using degradation and TCR files provided by Qualitau ACE and MIRA equipments.
All files should be placed in the folder where the function is run.
The function matches file names to create the exportfiles but is able to detect
if several conditions are stored in one file. One exportfile per condition is created.
If an exportfile is already present, user is asked if old files should be kept, replaced or if data have to be merged.
Times to failure are converted in seconds if they are provided in hours.

\bold{BlackAnalysis}: Extract Black's parameters from a set of electromigration experiments.
The experimental data as well as the resulting model are displayed and
can be saved. Extracted parameters are saved in a fit.txt file.}
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
