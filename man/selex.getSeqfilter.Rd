\name{selex.getSeqfilter}
\alias{selex.getSeqfilter}
\title{Display sequence filter attributes}
\description{Display all regular expressions used in a sequence filter object.
}
\usage{
selex.getSeqfilter(regex)
}
\arguments{
\item{regex}{A sequence filter object.}
}
\value{
A character string containing the regular expressions used in a sequence filter object.
}
\seealso{
\code{\link{selex.seqfilter}}
}
\examples{
\dontshow{
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)
regex = selex.seqfilter(kmerExcludeRegex="^A.{14}G", viewIncludeRegex="^T[ACTG]{14}T") 
}
selex.getSeqfilter(regex)
}
\keyword{misc}
\keyword{methods}
