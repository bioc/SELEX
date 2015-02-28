\name{selex.getRound0}
\alias{selex.getRound0}
\title{Obtain round zero sample handle}
\description{A function used to return the sample handle of the round zero file associated with the input \code{sample}.}
\usage{
selex.getRound0(sample)
}
\arguments{
\item{sample}{A sample handle.}
}
\value{
\code{selex.getRound0} returns a sample handle to the corresponding round zero file. The latter needs to be defined in the annotation table. If not, an error will be generated.
}
\examples{
\dontshow{
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)
sampleFiles = selex.exampledata(workDir)
selex.loadAnnotation(sampleFiles[3])
}
#Show currently visible samples
selex.sampleSummary()

#Return the matched roound zero file
r2 = selex.sample(seqName='R2.libraries', sampleName='ExdHox.R2', round=2)
r0 = selex.getRound0(r2)
selex.getAttributes(r0)
}
\keyword{misc}
\keyword{methods}
