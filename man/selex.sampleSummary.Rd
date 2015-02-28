\name{selex.sampleSummary}
\alias{selex.sampleSummary}
\title{
Show samples visible to the current SELEX session
}
\description{
A function used to output the sequencing run names, sample names, rounds, left and right barcodes, and file paths of all samples that have been currently loaded from sequencing run info files or defined using \code{\link{selex.defineSample}}. 
}
\usage{
selex.sampleSummary()
}
\details{
\code{selex.sampleSummary} displays all the samples currently visible to the SELEX session. A sample is `visible' to the current SELEX session when its metadata has been loaded or defined and can be used to create a sample handle from its sequencing run info name, sample name, and round (see \code{selex.sample}). Alternatively, a sample handle can also be specified using its row number as returned by \code{selex.sampleSummary}.
}
\value{
\code{selex.sampleSummary} returns a data frame containing sequencing run names, sample names, rounds, left and right barcodes, and file paths.
}
\seealso{
\code{\link{selex.defineSample}}, \code{\link{selex.getAttributes}}, \code{\link{selex.loadAnnotation}}, \code{\link{selex.sample}}, \code{\link{selex.saveAnnotation}}
}
\examples{
\dontshow{
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)
sampleFiles = selex.exampledata(workDir)
selex.loadAnnotation(sampleFiles[3])
}

selex.sampleSummary()

}
\keyword{misc}
\keyword{methods}
