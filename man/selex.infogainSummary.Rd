\name{selex.infogainSummary}
\alias{selex.infogainSummary}
\title{Summarize available information gain values}
\description{
This function returns the sample, Kmer length, information gain value, Markov model/type, and applied filters for either all computed information gain values or a specified sample in the current working directory.
}
\usage{
selex.infogainSummary(sample=NULL, displayFilter=FALSE)
}
\arguments{
\item{sample}{A sample handle to the sample for which information gain statistics should be returned.}
\item{displayFilter}{Display all applied sequence filter attributes.}
}
\value{
\code{selex.infogainSummary} returns a data frame containing the sample, Kmer length, information gain value, Markov model/type, and applied filters for all computed information gain values in the current working directory. If \code{sample} is provided, the above information is given for the specific sample only.
}

\seealso{
\code{\link{selex.infogain}}, \code{\link{selex.summary}}
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
k = selex.kmax(sample=r0.split$test)
mm = selex.mm(sample=r0.split$train, order=NA, crossValidationSample=r0.split$test)
ig =  selex.infogain(sample=r2, k=8, mm)
}
selex.infogainSummary()
}
\keyword{misc}
\keyword{methods}
