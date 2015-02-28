\name{selex.kmax}
\alias{selex.kmax}
\title{Calculate kmax for a dataset}
\description{
This function returns the longest oligonucleotide length k such that all DNA sequences of length k (`K-mers') are found at least a minimum count number of times for the given sample.
}
\usage{
selex.kmax(sample, threshold=100, seqfilter=NULL)
}
\arguments{
  \item{sample}{A sample handle.}
  \item{threshold}{The minimum count to be used.}
  \item{seqfilter}{A sequence filter object to include/exclude sequences that are read in from the FASTQ file.}
}
\details{
The kmax value is used to build the Markov model training and cross-validation datasets. While \code{\link{selex.mm}} contains a default kmax constructor, running \code{selex.kmax} can be useful in building analysis-specific Markov models. \cr \cr
\code{selex.kmax} discovers the kmax value by building K-mer count tables; after completion, the K-mer count tables can be viewed using \code{\link{selex.counts}}. When a new \code{seqfilter} object is provided, the kmax value is recomputed. See \code{\link{selex.seqfilter}} for more details.
}
\value{
\code{selex.kmax} returns the kmax value.
}
\seealso{
\code{\link{selex.counts}}, \code{\link{selex.mm}}, \code{\link{selex.sample}}, \code{\link{selex.seqfilter}}
}
\examples{
\dontshow{
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)
sampleFiles = selex.exampledata(workDir)
selex.loadAnnotation(sampleFiles[3])
selex.sampleSummary()
r0 = selex.sample(seqName="R0.libraries", sampleName="R0.barcodeGC", round=0)
}
kmax = selex.kmax(sample=r0, threshold=50)
}
\keyword{misc}
\keyword{manip}
\keyword{methods}
