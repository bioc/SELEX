\name{selex.countSummary}
\alias{selex.countSummary}
\title{Summarize available K-mer count tables}
\description{
This function returns the sample, kmer length, offset, minimum/maximum count, total count, and applied filters for either all count tables or a specified sample in the current working directory.
}
\usage{
selex.countSummary(sample=NULL, displayFilter=FALSE)
}
\arguments{
\item{sample}{A sample handle to the sample for which K-mer count statistics should be returned.}
\item{displayFilter}{Display all applied sequence filter attributes.}
}
\details{
Minimum count refers to the lowest count observed for a kmer of length k for a given sample. Total count is the sum of counts over all kmers of length k for a given sample.
}
\value{
\code{selex.countSummary} returns a data frame containing the sample, kmer length, offset, minimum/maximum count, total count, and applied filters for all count tables in the current working directory. If \code{sample} is provided, the above information is given for the specific sample only.
}
\seealso{
\code{\link{selex.counts}}, \code{\link{selex.summary}}
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
t1 =  selex.counts(sample=r2, k=2, offset=14, markovModel=NULL)
t2 =  selex.counts(sample=r2, k=4, markovModel=mm)
}
selex.countSummary()
}
\keyword{misc}
\keyword{methods}
