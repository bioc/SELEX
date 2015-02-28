\name{selex.mm}
\alias{selex.mm}
\title{Build or retrieve a Markov model}
\description{
A function used to compute and store a Markov model built by training and cross-validating on specified files. It returns a Markov model handle to conveniently reference the calculated model in other SELEX functions.
}
\usage{
selex.mm(sample, order=NA, crossValidationSample=NULL, Kmax= NULL, 
  seqfilter=NULL, mmMethod="DIVISION", mmWithLeftFlank=FALSE)
}
\arguments{
  \item{sample}{A sample handle to the training dataset.}
  \item{order}{The order of the Markov model to be built. If \code{NA}, builds Markov models of all orders from \code{0} to \code{Kmax-1}, selecting the one with highest R^2.}
  \item{crossValidationSample}{A sample handle to the cross-validation dataset. If \code{NULL}, no R^2 value will be provided.}
  \item{Kmax}{The K-mer length to determine model fit on. If \code{NA}, automatically finds kmax with a minimum count set to 100. See \code{\link{selex.kmax}} and `Details'.}
  \item{seqfilter}{A sequence filter object to include/exclude sequences that are read in from the FASTQ file.}  
  \item{mmMethod}{A character string indicating the algorithm used to evaluate the Markov model conditional probabilities. Can be either \code{"DIVISION"} (default) or \code{"TRANSITION"}. See `Details'.}
  \item{mmWithLeftFlank}{Predict expected counts by considering the sequences in the left flank of the variable region. See `Details'.}
}
\details{
\code{selex.mm} builds an N-th order Markov model by training the model on the \code{sample} dataset and tests its accuracy by attempting to predict the K-mer counts of the cross-validation dataset. This K-mer length is determined either by setting \code{Kmax} or by an internal call to \code{selex.kmax}. If a cross-validation dataset does not exist, \code{\link{selex.split}} can be used to partition a dataset into testing and training datasets. \cr \cr
The Markov model uses conditional probabilities to predict the probability of observing a given K-mer. For example, let us consider the following 6-mer sequence: \cr \cr 
AGGCTA \cr \cr
If a fourth-order Markov model were to predict the prior observation probability of the above sequence, the following would have to be evaluated:
\deqn{p(AGGCTA) = p(AGG)~p(C|AGG)~p(T|GGC)~p(A|GCT)}{p(AGGCTA) = p(AGG) p(C|AGG) p(T|GGC) p(A|GCT)}
These conditional probabilities can be evaluated using one of two algorithms. The \code{TRANSITION} algorithm relies on K-mer counts of length equal to the order of the Markov model. For the example above, 
\deqn{p(C|AGG) = \frac{count(AGGC)}{count(AGGA)+count(AGGC)+count(AGGG)+count(AGGT)}}{p(C|AGG) = count(AGGC) / [count(AGGA) + count(AGGC) + count(AGGG) + count(AGGT)]}
The \code{DIVISION} algorithm relies on K-mer frequencies for lengths equal to Markov model order and order-1. For the example above,
\deqn{p(C|AGG) = \frac{P(AGGC)}{P(AGG)} = \frac{\frac{count(AGGC)}{total 4-mer counts}}{\frac{(count(AGG)}{total 3-mer counts}} }{ 
p(C|AGG) = P(AGGC) / P(AGG) = [count(AGGC) / (total 4-mer counts)] / [count(AGG) / (total 3-mer counts)] }
The the flanking sequences in the left flank of the variable region can be taken into consideration when predicting prior observation probabilities by setting \code{mmWithLeftFlank} to \code{TRUE}. If the left barcode of the sequence in the example above is TGG, P(AGG) in the prior probability becomes
\deqn{p(AGG) = p(A|TGG)~p(G|GGA)~p(G|GAG).}{p(AGG) = p(A|TGG) p(G|GGA) p(G|GAG).}
After computation, the Markov model is stored in the working directory. When a new \code{seqfilter} object is provided, the Markov model is reconstructed. See \code{\link{selex.seqfilter}} for more details. See `References' for more details regarding the model construction process. \cr \cr
\code{selex.mmSummary} can be used to view the R^2 values for all orders that have been computed for the Markov model. If \code{crossValidationSample} is \code{NULL}, the resulting Markov model will not be displayed by \code{selex.mmSummary}.

}
\value{
\code{selex.mm} returns a Markov model handle.
}
\seealso{
\code{\link{selex.counts}}, \code{\link{selex.infogain}}, \code{\link{selex.kmax}}, \code{\link{selex.mmSummary}}, \code{\link{selex.run}}, \code{\link{selex.split}}
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
mm = selex.mm(sample=r0.split$train, order=NA, crossValidationSample=r0.split$test)

# View Markov model R^2 values
selex.mmSummary()
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
