\name{selex.mmSummary}
\alias{selex.mmSummary}
\title{Summarize Markov model properties}
\description{
This function returns sample, order, Markov model type, R^2 value, cross validation sample/length, and applied filters for either all computed Markov models or a specified sample in the current working directory.
}
\usage{
selex.mmSummary(sample=NULL, displayFilter=FALSE)
}
\arguments{
\item{sample}{A sample handle to the sample for which Markov model information should be returned.}
\item{displayFilter}{Display all applied sequence filter attributes.}
}
\value{
\code{selex.mmSummary} returns a data frame containing the sample, order, Markov model type, R^2 value, cross validation sample/length, and applied filters for all computed Markov models in the current working directory. If \code{sample} is provided, the above information is given for the specific sample only.
}
\seealso{
\code{\link{selex.mm}}, \code{\link{selex.summary}}
}
\examples{
\dontshow{
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)
sampleFiles = selex.exampledata(workDir)
selex.loadAnnotation(sampleFiles[3])

r0 = selex.sample(seqName="R0.libraries", sampleName="R0.barcodeGC", round=0)
r0.split = selex.split(sample=r0)

k = selex.kmax(sample=r0.split$test)
mm = selex.mm(sample=r0.split$train, order=NA, crossValidationSample=r0.split$test)
}
selex.mmSummary()
}
\keyword{misc}
\keyword{methods}
