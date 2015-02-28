\name{selex.affinities}
\alias{selex.affinities}
\title{
Construct a K-mer affinity table
}
\description{
A function used to calculate and return the affinities and affinity standard errors of K-mers of length \code{k} for a given dataset in addition to all the output provided by \code{\link{selex.counts}}. A Markov model is necessary for evaluation.
}
\usage{
selex.affinities(sample, k, minCount=100, top=-1, numSort=TRUE, offset=NULL,
	markovModel=NULL, seqfilter=NULL)
}
\arguments{
  \item{sample}{A sample handle to the dataset on which K-mer counting should be perfomed.}
  \item{k}{K-mer length(s) to be counted.}
  \item{minCount}{The minimum number of counts for a K-mer to be output.}
  \item{top}{Give the first N K-mers (by count).}
  \item{numSort}{Sort K-mers in descending order by count. If \code{FALSE}, K-mers are sorted alphabetically.}
  \item{offset}{Location of window for which K-mers should be counted for. If not provided, K-mers are counted across all windows.}
  \item{markovModel}{Markov model handle to use to predict previous round probabilities and expected counts.}
  \item{seqfilter}{A sequence filter object to include/exclude sequences that are read in from the FASTQ file.}  
}
\details{
When a new \code{seqfilter} object is provided, K-mer counting and affinity table construction is redone. See \code{\link{selex.seqfilter}} for more details. \cr \cr

See `References' for more details regarding K-mer counting and affinity calculation.
}
\value{
\code{selex.affinities} returns a data frame containing the K-mer sequence, observed counts, predicted prior observation probability, predicted prior observed counts, affinities, and standard errors. 
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
\code{\link{selex.counts}}, \code{\link{selex.infogain}}, \code{\link{selex.kmax}}, \code{\link{selex.mm}} 
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
r2Aff = selex.affinities(sample=r2, k=10, markovModel=mm)
}
\keyword{misc}
\keyword{manip}
\keyword{methods}
