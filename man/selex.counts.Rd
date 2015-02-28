\name{selex.counts}
\alias{selex.counts}
\title{
Construct or retrieve a K-mer count table
}
\description{
A function used to count and return the number of instances K-mers of length \code{k} appear within the \code{sample}'s variable regions. If an \code{offset} value is provided, K-mer counting takes place at a fixed position within the variable region of the read. If a Markov model is supplied, the expected count and the probability of observing the K-mer are also returned.
}
\usage{
selex.counts(sample, k, minCount=100, top=-1, numSort=TRUE, offset=NULL,
  markovModel=NULL, forceCalculation=FALSE, seqfilter=NULL, outputPath = "")
}
\arguments{
  \item{sample}{A sample handle to the dataset on which K-mer counting should be perfomed.}
  \item{k}{K-mer length(s) to be counted.}
  \item{minCount}{The minimum number of counts for a K-mer to be output.}
  \item{top}{Give the first N K-mers (by count).}
  \item{numSort}{Sort K-mers in descending order by count. If \code{FALSE}, K-mers are sorted in ascending order.}
  \item{offset}{Location of window for which K-mers should be counted for. If not provided, K-mers are counted across all windows.}
  \item{markovModel}{Markov model handle to use to predict previous round probabilities and expected counts.}
  \item{forceCalculation}{Forces K-mer counting to be performed again, even if a previous result exists.}
  \item{seqfilter}{A sequence filter object to include/exclude sequences that are read in from the FASTQ file.}  
  \item{outputPath}{Prints the computed K-mer table to a plain text file. This is useful when the number of unique K-mers in the dataset exceeds R's memory limit.}
}
\details{
The \code{offset} feature counts K-mers of length \code{k} \code{offset} bp away from the 5' end in the variable region. For example, if we have 16-mer variable regions and wish to count K-mers of length 12 at an offset of 3 bp, we are looking at the K-mers found only at the position indicated by the bolded nucleotides in the variable region: \cr \cr
5'   NNN\bold{NNNNNNNNNNNN}N   3' \cr \cr
Minimum count refers to the lowest count observed for a kmer of length k for a given sample. Total count is the sum of counts over all kmers of length k for a given sample. These statistics can be viewed for all K-mer lengths and samples counting was performed on using \code{\link{selex.countSummary}}. When a new \code{seqfilter} object is provided, K-mer counting is redone. See \code{\link{selex.seqfilter}} for more details.\cr \cr
See `References' for more details regarding the K-mer counting process.
}
\value{
\code{selex.counts} returns a data frame containing the K-mer sequence and observed counts for a given sample if a Markov model has not been supplied. \cr \cr
If a Markov model is supplied, a data frame containing K-mer sequence, observed counts, predicted previous round probability, and predicted previous round expected counts is returned. \cr \cr
If the number of unique K-mers exceeds R's memory limit, \code{selex.counts} will cause R to crash when returning a data frame containing the K-mers. The \code{outputPath} option can be used to avoid such a situation, as the Java code will directly write the table to a plain text file at the specified location instead.
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
\code{\link{selex.affinities}}, \code{\link{selex.countSummary}}, \code{\link{selex.infogain}}, \code{\link{selex.kmax}}, \code{\link{selex.mm}}, \code{\link{selex.run}}
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
# Kmer counting for a specific length on a given dataset
t1 = selex.counts(sample=r2, k=8, minCount=1)

# Kmer counting with an offset
t2 = selex.counts(sample=r2, k=2, offset=14, markovModel=NULL)

# Kmer counting with a Markov model (produces expected counts)
t3 = selex.counts(sample=r2, k=4, markovModel=mm)

# Display all available kmer statistics
selex.countSummary()
}
\keyword{misc}
\keyword{manip}
\keyword{methods}
