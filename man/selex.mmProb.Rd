\name{selex.mmProb}
\alias{selex.mmProb}
\title{Compute prior probability of sequence using Markov model}
\description{
A function used to calculate and return the prior observation probability of a DNA sequence \code{seqStr} as predicted by Markov model \code{markovModel}.
}
\usage{
selex.mmProb(seqStr, markovModel)
}
\arguments{
  \item{seqStr}{A character string representing the DNA sequence to be evaluated.}
  \item{markovModel}{A Markov model handle.}
}
\value{
\code{selex.mm} returns the prior observation probability of \code{seqStr}.
}
\seealso{
\code{\link{selex.mmProb}}
}
\examples{
\dontshow{
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)
sampleFiles = selex.exampledata(workDir)
selex.loadAnnotation(sampleFiles[3])
r0 = selex.sample(seqName="R0.libraries", sampleName="R0.barcodeGC", round=0)
r0.split = selex.split(sample=r0)
k = selex.kmax(sample=r0.split$test)
}
# Build the Markov model on the training dataset
mm = selex.mm(sample=r0.split$train, order=NA, crossValidationSample=r0.split$test)

# Compute expected Markov model probability value for a random 12-mer
selex.mmProb(seqStr="ATTGCAGACTTG", markovModel=mm)

}
\references{
Slattery, M., Riley, T.R., Liu, P., Abe, N., Gomez-Alcala, P., Dror, I., Zhou, T., Rohs, R., Honig, B., Bussemaker, H.J.,and Mann, R.S. (2011)
\href{http://www.ncbi.nlm.nih.gov/pubmed/22153072}{{C}ofactor binding evokes latent differences in {D}{N}{A} binding specificity between {H}ox proteins.}
Cell 147:1270--1282.
   
Riley, T.R., Slattery, M., Abe, N., Rastogi, C., Liu, D., Mann, R.S., and Bussemaker, H.J. (2014)
\href{http://www.ncbi.nlm.nih.gov/pubmed/25151169}{{S}{E}{L}{E}{X}-seq: a method for characterizing the complete repertoire of binding site preferences for transcription factor complexes.}
Methods Mol. Biol. 1196:255--278.
}
\keyword{misc}
\keyword{manip}
\keyword{methods}
