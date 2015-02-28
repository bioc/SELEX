\name{selex.defineSample}
\alias{selex.defineSample}
\title{Define annotation for an individual sample}
\description{
A function used to manually load SELEX sample metadata and make it visible to the current session, which can then be used to create a sample handle (see \code{\link{selex.sample}}. It is functionally identical to \code{selex.loadAnnotation}, except constrained to take manual inputs for a single sample. 
}
\usage{
selex.defineSample(seqName, seqFile=NULL, sampleName, round, varLength, 
  leftBarcode, rightBarcode, leftFlank=NULL, rightFlank=NULL, 
  round0SeqName=NULL, round0SampleName=NULL)
}
\arguments{
  \item{seqName}{Desired sequencing run name}
  \item{seqFile}{Path to the FASTQ file containing the sample. Can be \code{NULL} if \code{seqFile} has been previously specified for the \code{seqName}.}
  \item{sampleName}{Desired sample name}
  \item{round}{Sample round}
  \item{varLength}{Length of the variable region}
  \item{leftBarcode}{Left barcode sequence}
  \item{rightBarcode}{Right barcode sequence}
  \item{leftFlank}{Left flank sequence}
  \item{rightFlank}{Right flank sequence}
  \item{round0SeqName}{The sequencing run name of the round 0 data associated with this sample}
  \item{round0SampleName}{The sample name of the round 0 data associated with this sample}
}
\details{
\code{selex.defineSample} should be used for rapid testing or when it is not worthwhile to generate a generate a sample annotation file. Unlike \code{selex.loadAnnotation}, the unique name requirement is relaxed, requiring only unique \code{seqName}, \code{seqFile}, and \code{round} combinations. Only one \code{seqFile} can be specified for a given \code{seqName}; for example, the following will throw an error: \cr \cr
\code{selex.defineSample('Seq1', 'Seq1.fastq.gz', round=1, rightBarcode='CCAGCTG', ...)} \cr
\code{selex.defineSample('Seq1', 'Seq2.fastq.gz', round=2, rightBarcode='CCAGCTG', ...)} \cr \cr
However, \code{seqFile} can be omitted if a new sample is being specified with the same \code{seqName}: \cr \cr
\code{selex.defineSample('Seq1', 'Seq1.fastq.gz', round=1, rightBarcode='CCAGCTG', ...)} \cr
\code{selex.defineSample('Seq1', NULL, round=1, rightBarcode='CCACGTC', ...)} \cr \cr
The sequencing run name and sample name of the round 0 file associated with a later-round sample can be provided to keep samples and their random pools in order.
}
\value{Not applicable}
\seealso{
\code{\link{selex.getAttributes}}, \code{\link{selex.loadAnnotation}}, \code{\link{selex.sample}}, \code{\link{selex.sampleSummary}}, \code{\link{selex.saveAnnotation}}
}
\examples{
\dontshow{
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir, verbose=FALSE, maxThreadNumber= 4)
exampleFiles = selex.exampledata(workDir)
}
selex.defineSample(seqName='R0.libraries', seqFile=exampleFiles[1],
                sampleName='R0.barcodeGC', round = 0, varLength = 16, 
                leftBarcode = 'TGG', rightBarcode= 'CCAGCTG')
}

\keyword{misc}
\keyword{methods}
