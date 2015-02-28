\name{selex.sample}
\alias{selex.sample}
\title{
Create a sample handle
}
\description{
A function used to create a sample handle for conveniently referencing visible samples in other SELEX functions. 
}
\usage{
selex.sample(seqName, sampleName, round, index=NULL) 
}
\arguments{
  \item{seqName}{Sequencing run name}
  \item{sampleName}{Sample name}
  \item{round}{Sample round}
  \item{index}{Row number of desired sample when using \code{selex.sampleSummary}. When used, overrides other variables.}
}
\value{
\code{selex.sample} returns a sample handle object.
}
\details{
A sample is considered `visible' to the current SELEX session when its metadata has been loaded or defined and the sample's sequencing run info name, sample name, and round is displayed when \code{\link{selex.sampleSummary}} is called. Alternatively, a row number from \code{selex.sampleSummary} can be used to reference a sample.\cr \cr
A sample handle is an object that contains a reference to one visible sample. This object can be passed to SELEX functions instead of explicitly passing all sample information (such as variable region length, barcodes, round, FASTQ file location, etc). 
}
\seealso{
\code{\link{selex.defineSample}}, \code{\link{selex.getAttributes}}, \code{\link{selex.loadAnnotation}}, \code{\link{selex.sampleSummary}}
}
\examples{
\dontshow{
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)
sampleFiles = selex.exampledata(workDir)
selex.loadAnnotation(sampleFiles[3])
}
r0 = selex.sample(seqName="R0.libraries", sampleName="R0.barcodeGC", round=0)
r2 = selex.sample(index=3)
}
\keyword{misc}
\keyword{method}
