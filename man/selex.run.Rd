\name{selex.run}
\alias{selex.run}
\title{
Run a standard SELEX analysis
}
\description{
A function used to, in one shot, 
\itemize{
  \item Determine kmax on the \code{crossValidationSample} with the minimum count determined by \code{minCount}
  \item Build a Markov model on the \code{trainingSample} and test it on the \code{crossValidationSample} with kmax length K-mers used to determine model fit, and constructed using \code{mmMethod}
  \item Calculate information gain for \code{infoRange} K-mer lengths on the \code{infoGainSample}, using the Markov model order with the highest R^2 to predict previous round values.
}
}
\usage{
selex.run(trainingSample, crossValidationSample, minCount=100, infoGainSample, 
  infoRange=NULL, mmMethod="DIVISION", mmWithLeftFlank=FALSE)
}
\arguments{
  \item{trainingSample}{A sample handle to the training dataset.}
  \item{crossValidationSample}{A sample handle to the cross-validation dataset.}
  \item{minCount}{The minimum count to be used.}
  \item{infoGainSample}{A sample handle to the dataset on which to perform the information gain analysis.}
  \item{infoRange}{The range of K-mer lengths for which the information gain should be calculated. If \code{NULL}, the range is automatically set to start from the optimal Markov model order + 1 to the length of the variable region. This is the same as \code{k} in \code{selex.infogain}.}
  \item{mmMethod}{A character string indicating the algorithm used to evaluate the Markov model conditional probabilities. Can be either \code{"DIVISION"} (default) or \code{"TRANSITION"}.}
  \item{mmWithLeftFlank}{Predict expected counts by considering the sequences in the left flank of the variable region.}
			
}
\value{Not applicable}
\details{
Please see the individual functions or `References' for more details.
}
\references{
Slattery, M., Riley, T.R., Liu, P., Abe, N., Gomez-Alcala, P., Dror, I., Zhou, T., Rohs, R., Honig, B., Bussemaker, H.J.,and Mann, R.S. (2011)
\href{http://www.ncbi.nlm.nih.gov/pubmed/22153072}{{C}ofactor binding evokes latent differences in {D}{N}{A} binding specificity between {H}ox proteins.}
Cell 147:1270--1282.
   
Riley, T.R., Slattery, M., Abe, N., Rastogi, C., Liu, D., Mann, R.S., and Bussemaker, H.J. (2014)
\href{http://www.ncbi.nlm.nih.gov/pubmed/25151169}{{S}{E}{L}{E}{X}-seq: a method for characterizing the complete repertoire of binding site preferences for transcription factor complexes.}
Methods Mol. Biol. 1196:255--278.
}
\seealso{
\code{\link{selex.counts}}, \code{\link{selex.countSummary}}, \code{\link{selex.infogain}}, \code{\link{selex.infogainSummary}}, \code{\link{selex.mm}}, \code{\link{selex.mmSummary}}
}
\examples{
#Initialize the SELEX package
#options(java.parameters="-Xmx1500M")
#library(SELEX) 

# Configure the current session
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)

# Extract sample data from package, including XML database
sampleFiles = selex.exampledata(workDir)

# Load all sample files using XML database
selex.loadAnnotation(sampleFiles[3])

# Create sample handles
r0 = selex.sample(seqName="R0.libraries", sampleName="R0.barcodeGC", round=0)
r2 = selex.sample(seqName='R2.libraries', sampleName='ExdHox.R2', round=2)

# Split the r0 sample into testing and training datasets
r0.split = selex.split(sample=r0)

# Run entire analysis
selex.run(trainingSample=r0.split$train, crossValidationSample=r0.split$test,
  infoGainSample=r2)

# Display results
selex.mmSummary()[,c(1,2,3,4,5,6)]
selex.infogainSummary()[,c(1,2,3,4,5)]
}
\keyword{misc}
\keyword{manip}
\keyword{methods}
