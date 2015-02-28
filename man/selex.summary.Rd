\name{selex.summary}
\alias{selex.summary}
\title{Display all count table, Markov model, and information gain summaries}
\description{
This function returns the summaries of either all the calculated count tables, Markov models, and information gain values in the current working directory or those for a specific sample. 
}
\usage{
selex.summary(sample=NULL, displayFilter=FALSE)
}
\arguments{
\item{sample}{A sample handle to the sample for which summaries should be returned.}
\item{displayFilter}{Display all applied sequence filter attributes.}
}
\value{
\code{selex.summary} returns a list containing the following components:
\item{countSummary}{a data frame containing the sample, kmer length, offset, minimum/maximum count, total count, and applied filters for either all count tables or a specified sample in the current working directory.}
\item{markovModel}{a data frame containing the sample, order, Markov model type, R^2 value, cross validation sample/length, and applied filters for either all computed Markov models or a specified sample in the current working directory.}
\item{informationGain}{a data frame containing the sample, Kmer length, information gain value, Markov model/type, and applied filters for either all computed information gain values or a specified sample in the current working directory.}
}
\seealso{
\code{\link{selex.counts}}, \code{\link{selex.countSummary}}, \code{\link{selex.infogain}}, \code{\link{selex.infogainSummary}}, \code{\link{selex.mm}}, \code{\link{selex.mmSummary}}
}
\examples{
\dontshow{
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)
sampleFiles = selex.exampledata(workDir)
selex.loadAnnotation(sampleFiles[3])
r0 = selex.sample(seqName="R0.libraries", sampleName="R0.barcodeGC", round=0)
r2 = selex.sample(seqName='R2.libraries', sampleName='ExdHox.R2', round=2)
r0.split = selex.split(sample=r0)
selex.run(trainingSample=r0.split$train, crossValidationSample=r0.split$test, 
  infoGainSample=r2)
}
selex.summary()
}
\keyword{misc}
\keyword{methods}
