\name{selex.jvmStatus}
\alias{selex.jvmStatus}
\title{
Display current JVM memory usage
}
\description{
A function that displays the current JVM memory usage. 
}
\usage{
selex.jvmStatus()
}
\details{
\code{selex.jvmStatus} is useful for verifying the proper installaztion and initialization of rJava. Setting \code{verbose=FALSE} in \code{\link{selex.config}} will disable terminal ouput and a call to \code{selex.jvmStatus} will display nothing. If the current JVM memory allocation is suboptimal, settings can be changed using the \code{memSize} option in \code{\link{selex.config}}.
}
\value{
\code{selex.jvmStatus} does not return a value. 
}
\seealso{
\code{\link{selex.config}}
}
\examples{
\dontshow{
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=TRUE, maxThreadNumber= 4)
}
selex.jvmStatus()
}
\keyword{misc}
\keyword{methods}
