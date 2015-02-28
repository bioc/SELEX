\name{selex.config}
\alias{selex.config}
\title{
Set SELEX system parameters
}
\description{A function used to set system parameters for the SELEX package. These parameters are stored in the .SELEX.properties file in the user's home directory. The existence of this file is checked for at runtime; if file does not exist, user is prompted to enter parameters. These parameters can be changed anytime by the user by invoking \code{selex.config}. \bold{It is important to set \code{java.parameters} before loading the SELEX library}. See `Details' for more information.
}
\usage{
selex.config(workingDir, verbose=FALSE, maxThreadNumber=NULL)
}

\arguments{
  \item{workingDir}{Physical location on disk for the current working directory. \bold{The full system path must be specified.}}
  \item{verbose}{Print more output to terminal. Useful for debugging.}
  \item{maxThreadNumber}{The maximum number of threads to be used while K-mer counting. If unspecified, 4 threads will be used.}
}
\details{
The working directory must be specified for every selex.run. This directory is used to store the output of all analyses performed by the selex package. If a certain calculation has already been performed in the given working directory, the previously computed values are returned. If the specified path does not exist, \code{selex.config} will attempt to create it. \cr \cr 
It is important to set the Java memory options before loading the SELEX library to allocate enough memory to the JVM. You can allocate memory by calling \code{options} and setting \code{java.parameters}. \code{-Xmx} must prefix the memory value; the code below allocates 1500 MB of RAM:\cr \cr 
\code{options(java.parameters="-Xmx1500M")} \cr \cr
If the memory option is unspecified, the default JVM memory setting will be used.
}
\value{Not applicable}
\seealso{
\code{\link{selex.jvmStatus}}, \code{\link{selex.setwd}}
}
\examples{
## Initialize SELEX
#options(java.parameters="-Xmx1500M")
#library(SELEX)

## Set working directory and verbose to true
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=TRUE, maxThreadNumber=4)
}
\keyword{misc}
\keyword{methods}
