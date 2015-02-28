\name{selex.loadAnnotation}
\alias{selex.loadAnnotation}
\title{Load a sample annotation file}
\description{
A function used to load sample metadata contained within a sample annotation file and make it visible to the current SELEX session. These samples can then be used to create sample handles (see \code{\link{selex.sample}}).
}
\usage{
selex.loadAnnotation(config_path, data_folder=NULL)
}
\arguments{
  \item{config_path}{Location on disk to the sample annotation file.}
  \item{data_folder}{Location on disk where FASTQ sample files are stored. This is either an absolute path, or relative to the location of the annotation file. If unspecified, it uses the parent folder of the annotation file.}
}
\value{Not applicable}
\details{
A sample annotation file is an XML file that acts as a database storing metadata for different SELEX experiments. Here, a SELEX experiment refers to a single SELEX round that has been sequenced. Such a database allows the user to explicitly store all relevant information in a structured manner for easy future access. \cr \cr
A sample annotation file is provided below. Every annotation file can contain multiple \code{SequencingRunInfo} instances; every instance within an annotation file must contain a unique name. If multiple annotation files are used in a given SELEX session, all such names must be unique. For example, the following annotation files A and B \cr \cr
File A \cr
\code{<SequencingRunInfo name="exdUbx.Run1">} and \cr
\code{<SequencingRunInfo name="exdUbx.Run2">} \cr \cr
File B \cr
\code{<SequencingRunInfo name="exdUbx.Run1">} \cr \cr
have a legal naming system if \bold{either} File A or File B is used in a single SELEX session, but have an invalid naming system if \bold{both} are used. In general, it is a good idea to ensure that every \code{SequencingRunInfo} name is unique. Every \code{SequencingRunInfo} instance references a single FASTQ file. The user has the option of providing additional metadata regarding the FASTQ file. \cr \cr 
A \code{SequencingRunInfo} instance can contain multiple \code{Sample}s. Every \code{Sample} name within a \code{SequencingRunInfo} instance must contain a unique name \bold{and} round combination. For example, \cr \cr
\code{<Sample name="exdLab", round="0">} and \cr
\code{<Sample name="exdLab", round="1">} \cr \cr
is a valid name \code{Sample} name combinations while \cr \cr
\code{<Sample name="exdLab", round="0">} and \cr
\code{<Sample name="exdLab", round="0">} \cr \cr
is not. Non Round 0 \code{Sample}s have the option of referencing a Round 0 file, working as a checking mechanism to prevent the wrong Round 0 sample from being used to analyze a later round sample. \cr \cr
Once samples have been loaded into the current SELEX session, a sample handle can be generated using the \code{SequencingRunInfo} name, \code{Sample} name, and \code{Round} number. Sample handles make it easier to reference individual \code{Sample}s while running an analysis. See \code{\link{selex.sample}} for more information.
}
\note{
Sample annotation files are structured as follows: \cr \cr
<?xml version="1.0" encoding="UTF-8"?> \cr \cr
  <SELEXSequencingConfig  xmlns:xsi="http://www.columbia.edu/"  xsi:noNamespaceSchemaLocation="selex.xsd"> \cr \cr
  <SequencingRunInfo name="exdUbx.exdScr.0">	 <!-- information needed for differentiating multiple sequencing info instances --> \cr
	 \cr
		<DataFile>/Users/Documents/Data/Run1.fastq.gz</DataFile> <!-- absolute or relative path -->
		\cr <SequencingPlatform>Illumina</SequencingPlatform>  	      <!-- #optional --> 
		\cr <ResearcherName>John Smith</ResearcherName>  	            <!-- #optional --> 
		\cr <ResearcherEmail>jsmith@columbia.edu</ResearcherEmail>  	<!-- #optional --> 
		\cr <SequencingFacilityName>Columbia University Genome Center</SequencingFacilityName>  	<!-- #optional -->
		\cr <SequencingFacilityEmail>cugc@columbia.edu</SequencingFacilityEmail>  	          <!-- #optional --> 
		\cr <Description>Ubx/Scr Round 0 Probes</Description>  	<!-- #optional -->
		\cr <Notes>Our first SELEX Run</Notes>  	              <!-- #optional --> 
			 \cr \cr
		<Sample name="barcodeCCAGCTG.v1" round="0"> \cr
			<Protein>Probes</Protein> \cr
			<Concentration></Concentration>		<!-- #optional --> 
			\cr <VariableRegionLength>16</VariableRegionLength> \cr
			<LeftFlank>GTTCAGAGTTCTACAGTCCGACGATCTGG</LeftFlank>   \cr
			<RightFlank>CCAGCTGTCGTATGCCGTCTTCTGCTTG</RightFlank>  \cr
			<LeftBarcode>TGG</LeftBarcode> \cr
			<RightBarcode>CCAGCTG</RightBarcode> \cr
			<Round0></Round0> \cr
			<Notes></Notes>  	                <!-- #optional -->
 \cr </Sample> \cr
	 \cr
		<Sample name="barcodeCCACGTC.v1" round="0"> \cr
			<Protein>Probes</Protein> \cr
			<Concentration></Concentration>		
			\cr <VariableRegionLength>16</VariableRegionLength> \cr
			<LeftFlank>GTTCAGAGTTCTACAGTCCGACGATCTGG</LeftFlank>   \cr
			<RightFlank>CCACGTCTCGTATGCCGTCTTCTGCTTG</RightFlank>  \cr
			<LeftBarcode>TGG</LeftBarcode> \cr
			<RightBarcode>CCACGTC</RightBarcode> \cr
			<Round0></Round0> \cr
			<Notes></Notes>                   
		\cr </Sample> \cr
 \cr
	</SequencingRunInfo> \cr
 \cr <!-- #New FASTQ file below -->
 \cr \cr <SequencingRunInfo name="exdUbx.exdScr.L.2"> \cr
	\cr
		<DataFile>/Users/Documents/Data/Run2.fastq.gz</DataFile>
  	\cr <SequencingPlatform>Illumina</SequencingPlatform>  	      
		\cr <ResearcherName>John Smith</ResearcherName>  	             
		\cr <ResearcherEmail>jsmith@columbia.edu</ResearcherEmail>  	
		\cr <SequencingFacilityName>Columbia University Genome Center</SequencingFacilityName> 
		\cr <SequencingFacilityEmail>cugc@columbia.edu</SequencingFacilityEmail>  	             
		\cr <Description>Ubx/Scr Round 2</Description>  	   
		\cr <Notes>Our first SELEX Run</Notes>  	        
		\cr \cr 	
		<Sample name="barcodeCCAGCTG.v1.low" round="2"> \cr
			<Protein>hmExdUbx</Protein> \cr
			<Concentration>low</Concentration> \cr
			<VariableRegionLength>16</VariableRegionLength> \cr
			<LeftFlank>GTTCAGAGTTCTACAGTCCGACGATCTGG</LeftFlank> \cr  
			<RightFlank>CCAGCTGTCGTATGCCGTCTTCTGCTTG</RightFlank> \cr
			<LeftBarcode>TGG</LeftBarcode> \cr
			<RightBarcode>CCAGCTG</RightBarcode> \cr
			<Round0 sequencingName="exdUbx.exdScr.0" sampleName="barcodeCCAGCTG.v1"/> \cr
			<Notes></Notes> \cr
		</Sample>\cr
	\cr
	</SequencingRunInfo>\cr \cr
</SELEXSequencingConfig>
}
\seealso{
\code{\link{selex.defineSample}}, \code{\link{selex.getAttributes}}, \code{\link{selex.sample}}, \code{\link{selex.sampleSummary}}, \code{\link{selex.saveAnnotation}}
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

# Create a sample handle
r0 = selex.sample(seqName="R0.libraries", sampleName="R0.barcodeGC", round=0)

# Use the sample handle to display sample properties
selex.getAttributes(r0)
}
\keyword{misc}
\keyword{methods}
