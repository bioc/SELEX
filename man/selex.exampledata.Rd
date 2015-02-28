\name{selex.exampledata}
\alias{selex.exampledata}
\title{Extract example data files}
\description{A function used to extract the sample files embedded in this package. The package contains 2 FASTQ files (\code{R0.fastq.gz} and \code{ExdHox.R2.fastq.gz}, and an XML configuration file (\code{config.xml}).}
\usage{
selex.exampledata(outputFolder)
}
\arguments{
\item{outputFolder}{Location on disk where the files should be extracted.}
}
\value{Not applicable}
\seealso{
\code{\link{selex.defineSample}}
}
\examples{
#Initialize the SELEX package
#options(java.parameters="-Xmx1500M")
#library(SELEX) 

# Configure the current session
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)

# Extract sample data from package, including XML database
exampleFiles = selex.exampledata(workDir)
}
\keyword{misc}
\keyword{methods}
