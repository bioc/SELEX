\name{selex.split}
\alias{selex.split}
\title{Randomly split a dataset}
\description{
A function used to randomly split and store a dataset into smaller fractions as defined by \code{ratios}. This is useful in generating cross-validation datasets for Markov Model testing when none are available.
}
\usage{
selex.split(sample, ratios=NA)
}
\arguments{
  \item{sample}{A sample handle to the dataset to be split.}
  \item{ratios}{The fractions with with new datasets should be generated. If \code{NA}, ratios is set to \code{c(.5,.5)} and \code{sample} is evenly split into two datasets.}
}
\details{
\code{selex.split} saves the newly generated datasets in the current working directory and makes them visible to the current SELEX session. Naming convention of new samples: New sequence name =  <seqName>.split.x ; New sample name =  <sampleName>.split.x .
}
\value{Not applicable}
\seealso{
\code{\link{selex.defineSample}}, \code{\link{selex.loadAnnotation}}, \code{\link{selex.mm}}, \code{\link{selex.sampleSummary}}
}
\examples{
\dontshow{
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)
sampleFiles = selex.exampledata(workDir)
selex.loadAnnotation(sampleFiles[3])
r0 = selex.sample(seqName="R0.libraries", sampleName="R0.barcodeGC", round=0)
}
# Split the r0 sample into testing and training datasets
r0.split = selex.split(sample=r0)

# show information about subsamples
names(r0.split)
r0.split$info
selex.getAttributes(r0.split$train)

# Display all currently loaded samples
selex.sampleSummary() 


}
\keyword{misc}
\keyword{manip}
\keyword{methods}
