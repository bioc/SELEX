\name{selex.kmerPSFM}
\alias{selex.kmerPSFM}
\title{Construct a PSFM from a K-mer table}
\description{A function used to calculate and return the PFSM (Position Specific Frequency Matrix) for a K-mer table of length \code{k} from \code{sample}. If an \code{offset} value is provided, K-mer counting takes place at a fixed position within the variable region.
}
\usage{
selex.kmerPSFM(sample, k, offset=NULL)
}
\arguments{
\item{sample}{A sample handle to the dataset on which K-mer counting should be perfomed.}
\item{k}{K-mer length(s) to be counted.}
\item{offset}{Location of window for which K-mers should be counted for. If not provided, K-mers are counted across all windows. \code{offset} starts from 0.}
}
\details{
A K-mer table will be constructed for the specified \code{sample}, length \code{k}, and \code{offset} if it does not already exist. 

The output can be used by the \code{seqLogo} package to create a sequence logo.
}
\value{
\code{selex.kmerPSFM} returns a matrix containing the frequences for each base at every position.  
}
\seealso{
\code{\link{selex.counts}}, \code{\link{selex.fastqPSFM}}
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
# Build the PSFM
psfm1 = selex.kmerPSFM(sample=r0, k=8, 0)

# Can make sequence logos using the seqLogo package:
#library(seqLogo)
#seqLogo(makePWM(t(psfm1)))
}
\keyword{misc}
\keyword{manip}
\keyword{methods}
