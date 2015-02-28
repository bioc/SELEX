\name{selex.infogain}
\alias{selex.infogain}
\title{Compute or retrieve information gain between rounds}
\description{
A function used to compute and store the information gain for various K-mer lengths on \code{sample} using \code{markovModel} to predict prior probabilities.
}
\usage{
selex.infogain(sample, k=NULL, markovModel, seqfilter=NULL, 
  checkBarcode=TRUE)
}
\arguments{
  \item{sample}{A sample handle to the dataset on which to perform the analysis.}
  \item{k}{The range of K-mer lengths for which the information gain should be calculated. If \code{NA}, the range is automatically set to start from the optimal Markov model order + 1 to the length of the variable region.}
  \item{markovModel}{A Markov model handle.}
  \item{seqfilter}{A sequence filter object to include/exclude sequences that are read in from the FASTQ file.}  
  \item{checkBarcode}{Checks to see if the sample barcodes used to construct the Markov model match those in the current sample, and prevents computation if they do not match.}  
}
\details{
\code{\link{selex.infogainSummary}} is required to view the computed information gain values. When a new \code{seqfilter} object is provided, the information gain analysis is redone. See \code{\link{selex.seqfilter}} for more details. \cr \cr
See `References' for more on the computation of information gain values.
}
\value{
\code{selex.infogain} returns the highest information gain value.
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
\code{\link{selex.infogainSummary}}, \code{\link{selex.mm}}, \code{\link{selex.run}}
}
\examples{
\dontshow{
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)
sampleFiles = selex.exampledata(workDir)
selex.loadAnnotation(sampleFiles[3])
r0 = selex.sample(seqName="R0.libraries", sampleName="R0.barcodeGC", round=0)
r2 = selex.sample(seqName='R2.libraries', sampleName='ExdHox.R2', round=2)
r0.split = selex.split(sample=r0)
k = selex.kmax(sample=r0.split$test)
mm = selex.mm(sample=r0.split$train, order=NA, crossValidationSample=r0.split$test)
}

# Calculate information gain for a fixed range
ig1 =  selex.infogain(sample=r2, k=c(8:10), markovModel=mm)

# View the results
selex.infogainSummary()[,c(1,2,3,4,5)]

# Now calculate for the default range
ig2 = selex.infogain(sample=r2, markovModel=mm)

# View the results again
selex.infogainSummary()[,c(1,2,3,4,5)]
}
\keyword{misc}
\keyword{manip}
\keyword{methods}
