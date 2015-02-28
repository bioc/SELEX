\name{selex.revcomp}
\alias{selex.revcomp}
\title{
Create forward-reverse complement data pairs 
}
\description{
A function used to find and return the reverse complement of K-mers and the values associated with them. It is useful for calculating forward/reverse complement symmetrized values. 
}
\usage{
selex.revcomp(kmer,value)
}
\arguments{
  \item{kmer}{A string array representing K-mers.}
  \item{value}{An array of associated values.}
}
\seealso{
\code{\link{selex.affinities}}, \code{\link{selex.counts}}
}
\details{
\code{selex.revcomp} finds and returns the reverse complement and associated value of every input K-mer, if it exists. If a reverse complement does not exist for a given K-mer, it is removed from the output. For example, consider the following K-mer and value arrays:
\tabular{ccc}{
ACGT \tab \tab .34 \cr
GCTA \tab \tab .22 \cr
CGAC \tab \tab .98 \cr
TAGC \tab \tab .19 
}
The output of \code{selex.revcomp} will be:
\tabular{ccccccc}{
ACGT \tab \tab .34 \tab \tab ACGT \tab \tab .34 \cr
GCTA \tab \tab .22 \tab \tab TAGC \tab \tab .19 \cr
TAGC \tab \tab .19 \tab \tab GCTA \tab \tab .22 
}
}
\value{
\code{selex.revcomp} returns a data frame containing the original K-mers and values, along with their reverse complements and associated values.  
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
# Find round 2 affinities
r2Aff = selex.affinities(sample=r2, k=10, markovModel=mm)

# Find the reverse complement affinities and standard errors
Aff = selex.revcomp(kmer=r2Aff$Kmer, value=r2Aff$Affinity)
SE = selex.revcomp(kmer=r2Aff$Kmer, value=r2Aff$SE)

# Find the forward/reverse complement symmetrized Affinity and SE values
symAff = (Aff$Value+Aff$Reverse.Complement.Values)/2
symSE = sqrt((SE$Value^2+SE$Reverse.Complement.Values^2)/2)

# Final Result
final = data.frame(Kmer=Aff$Kmer, Affinity=Aff$Value, 
  SymmetrizedAffinity=symAff/max(symAff), SE=SE$Value, 
  SymmetrizedSE=symSE/max(symAff))
final = final[order(-final$SymmetrizedAffinity),]
}
\keyword{misc}
\keyword{manip}
\keyword{methods}
