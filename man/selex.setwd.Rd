\name{selex.setwd}
\alias{selex.setwd}
\title{Set or change the working directory}
\description{
A function that changes the current working directory.
}
\usage{
selex.setwd(path)
}
\arguments{
  \item{path}{Physical location on disk for the current working directory. \bold{The full system path must be specified.}}
}
\value{Not applicable}
\details{
The working directory must be specified for every selex.run. This directory is used to store the output of all analyses performed by the selex package. If a certain calculation has already been performed in the given working directory, the previously computed values are returned. If the specified path does not exist, \code{selex.config} will attempt to create it.
}
\seealso{
\code{\link{selex.config}}
}
\examples{
workDir = file.path(".", "SELEX_workspace")
# Change the current working directory
selex.setwd(path=workDir)

}
\keyword{misc}
\keyword{methods}
