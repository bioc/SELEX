\name{selex.seqfilter}
\alias{selex.seqfilter}
\title{Create a sequence filter}
\description{A function used to create a sequence filter object to conveniently and precisely include or exclude sequences from being counted or displayed. The filters are formed using Java regular expressions and can be used by a variety of functions within the package.}
\usage{
selex.seqfilter(variableRegionIncludeRegex=NULL, 
  variableRegionExcludeRegex=NULL, variableRegionGroupRegex=NULL, 
  kmerIncludeRegex=NULL, kmerExcludeRegex=NULL, kmerIncludeOnly=NULL,
  viewIncludeRegex=NULL, viewExcludeRegex=NULL, viewIncludeOnly=NULL)
}
\arguments{
\item{variableRegionIncludeRegex}{Include reads with variable regions containing this regular expression.}
\item{variableRegionExcludeRegex}{Exclude reads with variable regions containing this regular expression.}
\item{variableRegionGroupRegex}{Select subsequences of variable regions matching this regular expression.}
\item{kmerIncludeRegex}{Perform K-mer counting on variable regions containing this regular expression.}
\item{kmerExcludeRegex}{Perform K-mer counting on variable regions  \bold{not} containing this regular expression.}
\item{kmerIncludeOnly}{Perform K-mer counting on variable regions \bold{exactly} matching this regular expression.}
\item{viewIncludeRegex}{Display K-mers containing this reguar expression.}
\item{viewExcludeRegex}{Display K-mers \bold{not} containing this regular expression.}
\item{viewIncludeOnly}{Display K-mers \bold{exactly} matching this regular expression.}
}
\details{
The filters described by \code{selex.seqfilter} are used to filter sequences in the different stages of the K-mer counting process: read filtering, variable region filtering, and K-mer filtering. 
\tabular{ccccc}{
Read Filtering  \tab \tab  Variable Region Filtering   \tab \tab K-mer Filtering \cr
\code{variableRegionIncludeRegex}   \tab \tab \code{kmerIncludeRegex}  \tab \tab \code{viewIncludeRegex} \cr
\code{variableRegionExcludeRegex}   \tab \tab \code{kmerExcludeRegex}  \tab \tab \code{viewExcludeRegex} \cr
\code{variableRegionGroupRegex}     \tab \tab \code{kmerIncludeOnly}   \tab \tab \code{viewIncludeOnly}
}
Read filtering includes or excludes reads from the FASTQ file, acting as additional filters to those used to extract the variable regions. For example, consider an experimental design where the left barcode is TGG, right right barcode is TTAGC, and the variable region length is 10. FASTQ reads will be rejected unless they have the correct format; the sequences below represent hypothetical FASTQ reads:
\tabular{ccl}{
5' TGG NNNNNNNNNN TTAGC 3' \tab \tab template \cr 
5' TCG ATCAGTGGAC TTAGC 3' \tab \tab \bold{fails} (left match failed)\cr
5' TGG NAGGTCAGAC TTAGC 3' \tab \tab \bold{fails} (indeterminate base in variable region)\cr
5' TGG ATCAGTGGAC TTAGC 3' \tab \tab \bold{passes}
}
The read filter options then act as additional filters on the 10-bp variable region. \code{variableRegionGroupRegex} has the added functionality of selecting substrings from the variable region itself. Using the same example, \code{variableRegionGroupRegex} could be used to select 5-mers regions flanked by AA on the left and the right (or AA NNNNN AA):
\tabular{ccl}{
5' TCG ATCAGTGGAC TTAGC 3' \tab \tab \bold{fails} template (left match failed)\cr
5' TGG ATCAGTGGAC TTAGC 3' \tab \tab \bold{passes} template, \bold{fails} filter (no match)\cr
5' TGG TAAGTGCCAA TTAGC 3' \tab \tab \bold{passes} template and filter
}
When the \code{variableRegionGroupRegex} filter matches, \emph{only the subsequence} will be used in future counting. In the above example, this would be GTGCC. \cr \cr
After the variable regions have been extracted from the FASTQ file, the next step involves K-mer counting. Variable region filtering comes into play here, allowing or preventing K-mer counting on these sequences. Lastly, K-mer filtering determines what K-mers are returned or displayed in tables. \cr \cr
Any function utilizing sequence filters will recompute results if, for a given sample, new values for the read filtering or variable region filter options are provided.
}
\value{
\code{selex.seqfilter} returns a sequence filter object.
}
\seealso{
\code{\link{selex.affinities}}, \code{\link{selex.counts}}, \code{\link{selex.getSeqfilter}}, \code{\link{selex.infogain}}, \code{\link{selex.kmax}}, \code{\link{selex.mm}} 
}
\examples{
\dontshow{
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)
sampleFiles = selex.exampledata(workDir)
selex.loadAnnotation(sampleFiles[3])
r0 = selex.sample(seqName="R0.libraries", sampleName="R0.barcodeGC", round=0)
}

# Raw K-mer counts
my.counts1 = selex.counts(sample=r0, k=16, top=100)

# Include reads whose variable regions begin with TGTA
regex = selex.seqfilter(variableRegionIncludeRegex="^TGTA") 
my.counts2 = selex.counts(sample=r0, k=16, top=100, seqfilter=regex)

# Exclude reads whose variables regions begin with TGT
regex = selex.seqfilter(variableRegionExcludeRegex="^TGT")
my.counts3 = selex.counts(sample=r0, k=16, top=100, seqfilter=regex)

# Extract 13-bp substring from reads whose variable regions begin with TGT
regex = selex.seqfilter(variableRegionGroupRegex="^TGT([ACGT]{13})")
my.counts4 = selex.counts(sample=r0, k=13, top=100, seqfilter=regex)

# Extract 5-bp substring from reads whose variable regions begin with TGT
regex = selex.seqfilter(variableRegionGroupRegex="^TGT([ACGT]{5})")
my.counts5 = selex.counts(sample=r0, k=5, top=100, seqfilter=regex)

# Select variable regions beginning with A and ending with G
regex = selex.seqfilter(kmerIncludeRegex="^A.{14}G") 
my.counts6 = selex.counts(sample=r0, k=16, top=100, seqfilter=regex)

# Exclude variable regions beginning with A and ending with G 
regex = selex.seqfilter(kmerExcludeRegex="^A.{14}G") 
my.counts7 = selex.counts(sample=r0, k=16, top=100, seqfilter=regex)

# Exclude variable regions beginning with A and ending with G, and display
# 16-mers that start and end with T
regex = selex.seqfilter(kmerExcludeRegex="^A.{14}G", 
  viewIncludeRegex="^T[ACTG]{14}T") 
my.counts8 = selex.counts(sample=r0, k=16, top=100, seqfilter=regex)

# Exclude variable regions beginning with A and ending with G, and display
# 16-mers that do not start and end with T
regex = selex.seqfilter(kmerExcludeRegex="^A.{14}G", 
  viewExcludeRegex="^T[ACTG]{14}T") 
my.counts9 = selex.counts(sample=r0, k=16, top=100, seqfilter=regex)

# Only count variable regions containing TGTAAAATCAGTGCTG or TGTAAGTGGACTCTCG
regex = selex.seqfilter(kmerIncludeOnly=c('TGTAAAATCAGTGCTG', 
  'TGTAAGTGGACTCTCG')) 
my.counts10 = selex.counts(sample=r0, k=16, top=100, seqfilter=regex)

# Only display results for the K-mers TGTAAAATCAGTGCTG and TGTAAGTGGACTCTCG
regex = selex.seqfilter(viewIncludeOnly=c('TGTAAAATCAGTGCTG', 
  'TGTAAGTGGACTCTCG')) 
my.counts11 = selex.counts(sample=r0, k=16, top=100, seqfilter=regex)
}
\keyword{misc}
\keyword{methods}
