\name{selex.fastqPSFM}
\alias{selex.fastqPSFM}
\title{Construct a diagnostic PSFM for a FASTQ file}
\description{A function used to calculate and return the PFSM (Position Specific Frequency Matrix) for an entire FASTQ file, regardless of barcode or other sequence filtering.}
\usage{
selex.fastqPSFM(seqName)
}
\arguments{
\item{seqName}{A sequencing run name for the desired FASTQ file; this must match a sequencing run name of a sample currently visible to the SELEX session.}
}
\details{
The output can be used by the \code{seqLogo} package to create a sequence logo.
}
\value{
\code{selex.fastqPSFM} returns a matrix containing the frequences for each base at every position.  
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

# Make PSFMs for the two visible FASTQ files:
psfm1 = selex.fastqPSFM(seqName='R0.libraries')
psfm2 = selex.fastqPSFM(seqName='R2.libraries')

# Can make sequence logos using the seqLogo package:
#library(seqLogo)
#seqLogo(makePWM(t(psfm1)))
}
\keyword{misc}
\keyword{methods}
