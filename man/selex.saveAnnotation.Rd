\name{selex.saveAnnotation}
\alias{selex.saveAnnotation}
\title{Save sample annotations to file}
\description{
A function used to save sample metadata visible to the current SELEX session in a sample annotation file.
}
\usage{
selex.saveAnnotation(filePath)
}
\arguments{
  \item{filePath}{Location on disk to create the sample annotation file. \bold{The full system path must be specified.}}
}
\details{
A sample annotation file is an XML file that acts as a database storing metadata for different SELEX experiments. \code{selex.saveAnnotation} provides a convenient way to permanently store manually entered sample information using \code{\link{selex.defineSample}}. For more information on the XML format used to store the information, see \code{\link{selex.loadAnnotation}}.
}
\value{Not applicable}
\seealso{
\code{\link{selex.defineSample}}, \code{\link{selex.loadAnnotation}}, \code{\link{selex.sampleSummary}}
}
\examples{
\dontshow{
workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)
sampleFiles = selex.exampledata(workDir)
}
selex.defineSample(seqName='R0.libraries', seqFile=sampleFiles[1],
                sampleName='R0.barcodeGC', round=0, varLength=16, 
                leftBarcode='TGG', rightBarcode='CCAGCTG')

selex.saveAnnotation(paste0(workDir, "sample_annotations.xml"))
}
\keyword{misc}
\keyword{methods}
