\name{selex.getAttributes}
\alias{selex.getAttributes}
\title{
Display sample handle properties
}
\description{
A function used to output sample handle's sequencing run name, sample name, round, left and right barcode, file path, and variable region length of the sample handle.
}
\usage{
selex.getAttributes(sample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample}{A sample handle.}
}
\value{
\code{selex.getAttributes} returns a data frame containing the sequencing run name, sample name, round, left and right barcode, file path, and variable region length of the sample handle.
}
\seealso{
\code{\link{selex.defineSample}}, \code{\link{selex.loadAnnotation}}, \code{\link{selex.sample}}, \code{\link{selex.sampleSummary}}
}
\examples{
\dontshow{	
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)
sampleFiles = selex.exampledata(workDir)
selex.loadAnnotation(sampleFiles[3])
}
# Create a sample handle
r0 = selex.sample(seqName="R0.libraries", sampleName="R0.barcodeGC", round=0)

# Use the sample handle to display sample properties
selex.getAttributes(r0)
}
\keyword{misc}
\keyword{methods}
