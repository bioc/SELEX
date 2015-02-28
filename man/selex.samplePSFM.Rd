\name{selex.samplePSFM}
\alias{selex.samplePSFM}
\title{Construct a diagnostic PSFM for a FASTQ file}
\description{A function used to calculate and return the PFSM (Position Specific Frequency Matrix) for a sample, regardless of sequence filtering.}
\usage{
selex.samplePSFM(sample)
}
\arguments{
  \item{sample}{A sample handle to the dataset on which the PSFM is calculated.}
}
\details{
The output can be used by the \code{seqLogo} package to create a sequence logo.
}
\value{
\code{selex.samplePSFM} returns a matrix containing the frequences for each base at every position.  
}
\seealso{
\code{\link{selex.kmerPSFM}}
}
\examples{
\dontshow{
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)
sampleFiles = selex.exampledata(workDir)
selex.loadAnnotation(sampleFiles[3])
}
# Display all currently loaded samples
selex.sampleSummary() 

r0 = selex.sample(seqName="R0.libraries", sampleName="R0.barcodeGC", round=0)

# Make PSFMs for the two visible FASTQ files:
psfm1 = selex.samplePSFM(r0)

# Can make sequence logos using the seqLogo package:
#library(seqLogo)
#seqLogo(makePWM(t(psfm1)))
}
\keyword{misc}
\keyword{methods}
