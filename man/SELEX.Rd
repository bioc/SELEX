\name{SELEX}
\alias{SELEX-package}
\alias{Selex}
\alias{selex}
\docType{package}
\title{
	SELEX Package
}
\description{
Functions to assist in discovering transcription factor DNA binding specificities from SELEX-seq experimental data according to the Slattery \emph{et al.} paper. For a more comprehensive example, please look at the vignette. Sample data used in the Slattery, et. al. is stored in the \code{extdata} folder for the package, and can be accessed using either the base R function \code{\link{system.file}} or the package function \code{\link{selex.exampledata}}.

Functions available:
\tabular{ll}{
\code{\link{selex.affinities}} \tab Construct a K-mer affinity table \cr
\code{\link{selex.config}} \tab Set SELEX system parameters \cr
\code{\link{selex.counts}} \tab Construct or retrieve a K-mer count table \cr
\code{\link{selex.countSummary}} \tab Summarize available K-mer count tables\cr
\code{\link{selex.defineSample}} \tab Define annotation for an individual sample \cr
\code{\link{selex.exampledata}} \tab Extract example data files \cr
\code{\link{selex.fastqPSFM}} \tab Construct a diagnostic PSFM for a FASTQ file \cr
\code{\link{selex.getAttributes}} \tab Display sample handle attributes \cr
\code{\link{selex.getRound0}} \tab Obtain round zero sample handle \cr
\code{\link{selex.getSeqfilter}} \tab Display sequence filter attributes \cr
\code{\link{selex.infogain}} \tab Compute or retrieve information gain between rounds \cr
\code{\link{selex.infogainSummary}} \tab Summarize available information gain values \cr
\code{\link{selex.jvmStatus}} \tab Display current JVM memory usage \cr
\code{\link{selex.kmax}} \tab Calculate kmax for a dataset \cr
\code{\link{selex.kmerPSFM}} \tab Construct a PSFM from a K-mer table \cr
\code{\link{selex.loadAnnotation}} \tab Load a sample annotation file \cr
\code{\link{selex.mm}} \tab Build or retrieve a Markov model \cr
\code{\link{selex.mmProb}} \tab Compute prior probability of sequence using Markov model\cr
\code{\link{selex.mmSummary}} \tab Summarize Markov model properties \cr
\code{\link{selex.revcomp}} \tab Create forward-reverse complement data pairs \cr 
\code{\link{selex.run}} \tab Run a standard SELEX analysis \cr
\code{\link{selex.sample}} \tab Create a sample handle \cr
\code{\link{selex.sampleSummary}} \tab Show samples visible to the current SELEX session\cr
\code{\link{selex.saveAnnotation}} \tab Save sample annotations to file \cr
\code{\link{selex.seqfilter}} \tab Create a sequence filter \cr
\code{\link{selex.setwd}} \tab Set or change the working directory \cr
\code{\link{selex.split}} \tab Randomly split a dataset \cr
\code{\link{selex.summary}} \tab Display all count table, Markov model, and information gain summaries \cr
}
}
\details{
\tabular{ll}{
Package: \tab SELEX\cr
Type: \tab Package\cr
Version: \tab .99\cr
Date: \tab 2014-11-3\cr
License: \tab GPL\cr
}
}
\author{
Chaitanya Rastogi, Dahong Liu, and Harmen Bussemaker \cr \cr
Maintainer: Harmen Bussemaker \email{hjb2004@columbia.edu}
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

# Load & display all sample files using XML database
selex.loadAnnotation(sampleFiles[3])
selex.sampleSummary()

# Create sample handles
r0 = selex.sample(seqName="R0.libraries", sampleName="R0.barcodeGC", round=0)
r2 = selex.sample(seqName='R2.libraries', sampleName='ExdHox.R2', round=2)

# Split the r0 sample into testing and training sets
r0.split = selex.split(sample=r0)
r0.split

# Display all currently loaded samples
selex.sampleSummary() 

# Find kmax on the test dataset
k = selex.kmax(sample=r0.split$test)

# Build the Markov model on the training dataset
mm = selex.mm(sample=r0.split$train, order=NA, crossValidationSample=r0.split$test)
# See Markov model R^2 values
selex.mmSummary()

# Kmer counting with an offset
t1 =  selex.counts(sample=r2, k=2, offset=14, markovModel=NULL)
# Kmer counting with a Markov model (produces expected counts)
t2 =  selex.counts(sample=r2, k=4, markovModel=mm)
# Display all available kmer statistics
selex.countSummary()

# Calculate information gain
ig =  selex.infogain(sample=r2, k=8, mm)
# View information gain results
selex.infogainSummary()

# Perform the default analysis
selex.run(trainingSample=r0.split$train, crossValidationSample=r0.split$test, 
  infoGainSample=r2)

# View all stats
selex.summary()
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
