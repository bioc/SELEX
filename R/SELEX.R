selex.getAttributes <- function(sample) {
  seqRun = sample$getSequencingName()
  sampName = sample$getSampleName()
  round = sample$getSampleRound()
  output = data.frame(SequencingRun=seqRun,Sample=sampName,Round=round, stringsAsFactors=FALSE)
  
  info = J("main/SELEX")$getExperimentReferenceInfo(sample)
  a=sapply(info[1],.jevalArray)
  b=sapply(info[2],.jevalArray)
  a = c(names(output),a)
  output = data.frame(output,t(b))
  names(output) = a
  
  rm(info)
  return(output)
}

selex.affinities <-function(sample,k,minCount=100,top=-1, numSort=TRUE,offset=NULL,markovModel=NULL, seqfilter=NULL) {
  
  if(is.null(seqfilter))
  {
    seqfilter=selex.seqfilter()
  }
  round = selex.getAttributes(sample)$Round
  table = selex.counts(sample,k,minCount, top, numSort,offset,markovModel, FALSE, seqfilter)
  affinities = table$ObservedCount/table$ExpectedCount
  affinities = (affinities/max(affinities))^(1/round)
  SE = affinities*sqrt(2/table$ObservedCount)
  table = cbind(table,Affinity=affinities,SE=SE)
  return(table)
}

selex.revcomp <- function (kmer, value) {
  rc = as.character(reverseComplement(DNAStringSet(kmer)))
  idx = match(kmer, rc)
  rcValues = value[idx]
  data = data.frame(Kmer=kmer, Value=value, Reverse.Complement=rc, Reverse.Complement.Values=rcValues)
  data = data[complete.cases(data), ]
  return(data)
}

#done
selex.config <- function(workingDir=NULL, verbose=NULL, maxThreadNumber=NULL)
{

	if(!is.null(verbose))
	{
	  options("bioconductorSELEX.verbose"=verbose)
	  selex.verbose(verbose)
	}

	if(!is.null(maxThreadNumber) && maxThreadNumber>0)
	{
      options("bioconductorSELEX.maxthreads"=maxThreadNumber)
      J("base/Util")$setMaxThreadNumber(selex.getInt(maxThreadNumber))
    }

    if( !(is.null(workingDir) || is.na(workingDir)))
	{
	  options("bioconductorSELEX.workdir"=workingDir)
      selex.setwd(workingDir)
	}
}


#done #checked
selex.loadAnnotation <- function(config_path=NULL, data_folder=NULL)
{
	if(!is.null(config_path))
	{
    if(is.null(data_folder))
    {
      data_folder = "";
    }
		J("main/SELEX")$loadConfigFile(config_path, data_folder)
	}
}

#done
selex.defineSample <- function(seqName, seqFile=NULL, sampleName, round, varLength, 
                             leftBarcode, rightBarcode, leftFlank=NULL, rightFlank=NULL, 
                             round0SeqName = NULL, round0SampleName = NULL )
{
  if(is.null(seqName))
  {
    seqName=""
  }
  if(is.null(seqFile))
  {
    seqFile=""
  }
  
  if(is.null(round0SeqName))
  {
    round0SeqName="";
  }
  if(is.null(round0SampleName))
  {
    round0SampleName ="";
  }
  if(is.null(leftFlank))
  {
    leftFlank=leftBarcode;
  }
  if(is.null(rightFlank))
  {
    rightFlank=rightBarcode;
  }
	J("main/SELEX")$addSequenceInfo(seqName,seqFile,sampleName,
        	 selex.getInt(round),selex.getInt(varLength), 
           leftBarcode,rightBarcode,leftFlank,rightFlank,
        	 round0SeqName, round0SampleName )
}

#done
# should be the first method to call
selex.setwd <- function(path)
{
	J("main/SELEX")$setWorkingDirectory(path)
	#selex.workingDirectory = sprintf('%s/',path)
}

#done
#selex.run <- function()
#{
#	J("main/SELEX")$runSELEX()
#}

#done
selex.verbose <- function(trueFalse)
{
	trueFalseStr = "FALSE"
	if(trueFalse)
	{
		trueFalseStr = "TRUE"
	}
	J("base/DebugLog")$verbose(trueFalseStr)
}

#done
selex.getInt <- function(k)
{
	.jnew("java/lang/Integer",sprintf("%d",k))
}

#done
selex.getBoolean <- function(trueFalse)
{
	if(trueFalse)
	{
		.jnew("java/lang/Boolean", "true")
	}
	else
	{
		.jnew("java/lang/Boolean", "false")
	}
}
 

#
selex.unloadSamples <-function()
{
  J("main/SELEX")$unloadSamples()
}


#
selex.saveAnnotation <-function(filePath)
{
  J("main/SELEX")$exportConfigFile(filePath)
}

# 
selex.getRound0 <- function(sample)
{
  J("main/SELEX")$getRound0(sample)
}

#done #checked
selex.sample <- function(seqName, sampleName, round, index=NULL)
{
  if(is.null(index))
  {
    J("main/SELEX")$getExperimentReference(seqName,sampleName, selex.getInt(round))
  }
  else
  {
    row = selex.sampleSummary()[index,]
    J("main/SELEX")$getExperimentReference(row$seqName, row$sampleName, selex.getInt(strtoi(row$rounds) ))
  }
}

#done #checked
selex.sampleSummary <- function()
{
	results = J("main/SELEX")$showSamples()
	
	seqName = sapply(results[1],.jevalArray)
	sampleName = sapply(results[2],.jevalArray)
	rounds = sapply(results[3],.jevalArray)
	leftBarcode = sapply(results[4],.jevalArray)
	rightBarcode = sapply(results[5],.jevalArray)
	leftFlank = sapply(results[6],.jevalArray)
	rightFlank = sapply(results[7],.jevalArray)
	seqFile = sapply(results[8],.jevalArray)
	
	rm(results)
  
	mat = data.frame(seqName, sampleName, rounds, leftBarcode, rightBarcode, leftFlank, rightFlank, seqFile, stringsAsFactors=FALSE)
	mat = mat[order(mat[,1],mat[,2],mat[,3]), ]

	mat
}

#done #checked
selex.kmax <- function(sample,  threshold=100, seqfilter=NULL)
{
  
  if(is.null(seqfilter))
  {
    seqfilter=selex.seqfilter()
  }
  	kmax=0
	dataSet=sample
	sampleName = J("main/SELEX")$getSampleID(sample)
  
	mmOption = .jnew("base/MarkovModelOption")
	mmOption$setBuildMarkovModel(selex.getBoolean(TRUE));
  
	for (i in 1:16 ) 
	{
		kmax = i-1
		write(sprintf("Counting [%s][ K = %d ]", sampleName,i), stdout())
	    flush.console()
	
		result = J("main/SELEX")$doMinimalCounting(dataSet, selex.getInt(i), 
			mmOption, selex.getBoolean(FALSE), selex.getInt(threshold), seqfilter);
		
		if(result==0)
			break
		
	}
	
	kmax.key = sprintf("%s.kmax",sampleName )
	
	filter = seqfilter$getVariableRegionRegexFormattedString()
	ipsStats = J("main/SELEX")$getInputDataSetStats(dataSet,filter)
	write(sprintf("[ sample id : %s,  filter:  %s ]",
                  J("main/SELEX")$getSampleID(sample), filter), stdout())
  if(!exists("ipsStats"))
  {
     return -1
  }
	ipsStats$setKmax(selex.getInt(kmax));
	
	J("main/SELEX")$saveStats()
	write(sprintf("[ %s = %d ]",kmax.key,kmax), stdout())
	flush.console()
	
	kmax
}

#done #checked
selex.split <- function(sample, ratios=NA)
{
	if(is.null(ratios) || is.na(ratios))
	{
		ratios = c(0.5 , 0.5);
	}
	else
	{
		ratios = ratios / sum(ratios);
	}
	ratios
	
	ratioArray <- .jnew("java/util/ArrayList")
	for(i in 1:length(ratios))
	{
		ratioArray$add(.jnew("java/lang/Double", ratios[i]))
	}
	results =  J("main/SELEX")$splitDataSet(sample,"split",ratioArray)
	
	seqName=sapply(results[1],.jevalArray)
	sampleName=sapply(results[2],.jevalArray)
	round=sapply(results[3],.jevalArray)
	nReads=sapply(results[4],.jevalArray)

	train = selex.sample(seqName=seqName[1], sampleName=sampleName[1], round=round[1])
	test  = selex.sample(seqName=seqName[2], sampleName=sampleName[2], round=round[2])
	info = data.frame(seqName, sampleName, round, nReads, stringsAsFactors=FALSE)

	list(test=test, train=train, info=info)
}

#done #checked
selex.mm <- function( sample, order=NA, crossValidationSample=NULL, Kmax= NULL, 
                            seqfilter=NULL, mmMethod='DIVISION',
				 mmWithLeftFlank = FALSE)
{
	if(mmWithLeftFlank)
	{
	   mmMethod = paste(mmMethod,'WITH_LEFT_FLANK',sep='|') ;
	}

	if(is.null(seqfilter))
	{
	  seqfilter=selex.seqfilter()
	}
  
	trainingDataSet=sample 
	validationDataSet=crossValidationSample
	
	trainingSampleName = J("main/SELEX")$getSampleID(sample)
  if(!is.null(crossValidationSample))
	{
    validationSampleName = J("main/SELEX")$getSampleID(crossValidationSample)
  }
	
  filter = seqfilter$getVariableRegionRegexFormattedString()
	#
	kmax =0;
	if(is.null(Kmax) && !is.null(crossValidationSample))
	{
		kmax.key = sprintf("%s.kmax",validationSampleName )
		
		ipsStats = J("main/SELEX")$getInputDataSetStats(crossValidationSample,filter); 
    
		if(is.null(ipsStats))
		{
			sprintf("%s not found.", kmax.key)
			kmax = selex.kmax(sample=crossValidationSample, 100, seqfilter)
		}else {
		  kmax = ipsStats$getKmax()
		  if(is.null(kmax))
		  {
		    sprintf("%s not found.", kmax.key)
		    kmax = selex.kmax(sample=crossValidationSample,100,  seqfilter)
		  } else{
        kmax = strtoi(kmax)
      }
    }
		write(sprintf("[ %s = %d ]", kmax.key, kmax), stdout())
	} else {
		write(  sprintf("Overwriting Kmax = %d", Kmax) , stdout())
		kmax = Kmax;
	}
	
	indices = NA
	
	if( is.null(order) || is.na(order))
	{
		indices = 1:kmax
	} else {
		indices = order+1
	}
	
	maxR = 0
	markovLength = 0
	bestMM = NULL

	for (i in indices) 
	{
		write(sprintf("Counting [%s][ K = %d ]",trainingSampleName, i), stdout())
		flush.console()

		
		if(is.null(crossValidationSample))
		{
		    mm = J("main/SELEX")$trainMarkovModel(trainingDataSet, selex.getInt(i), 
                                              seqfilter, mmMethod)
		}else {
		    mm = J("main/SELEX")$trainMarkovModel(trainingDataSet, validationDataSet, 
		                selex.getInt(i), selex.getInt(kmax), seqfilter, mmMethod)
		}
				
		R = mm$getMarkovR2()
		if( (maxR < R) || is.null(bestMM))
		{
			maxR=R;
			markovLength=i;
			bestMM = mm
		}
	
	}

	write(sprintf("[ markovLength = %d ]",markovLength), stdout())
	write(sprintf("[ maxR = %f ]",maxR), stdout())
	write(sprintf("[ Model = %s ]", bestMM$toString() ), stdout())
	flush.console()

	bestMM
	
}

#done #checked

selex.infogain <- function(sample, k=NULL, markovModel, seqfilter=NULL, checkBarcode=TRUE)
{
  
  if(is.null(seqfilter))
  {
    seqfilter=selex.seqfilter()
  }
	dataSet=sample
  
	maxIG = 0
	
	minStart = markovModel$getMarkovLength()
	indices = NA
	
	if( is.null(k) || is.na(k))
	{
		indices = minStart:16
	} else {
		indices = k
	}
	
	for(i in indices)
	{
		write(sprintf("Counting [InfoGain][ K = %d ]",i), stdout())
		flush.console()
		IG = J("main/SELEX")$calculateInformationGain(dataSet, selex.getInt(i), markovModel, seqfilter,
                                                  selex.getBoolean(checkBarcode))	
		J("main/SELEX")$saveStats()
		
		if(maxIG<IG)
		{
			maxIG = IG
		}
	}
	
	maxIG
}


#done #checked
selex.counts <- function( sample, k, minCount=100, top=-1, numSort=TRUE, offset =NULL,
                          markovModel=NULL, forceCalculation =FALSE, seqfilter=NULL, outputPath = "" )
{	
  
  if(is.null(seqfilter))
  {
    seqfilter=selex.seqfilter()
  }
  	sampleName = J("main/SELEX")$getSampleID(sample)
	write(sprintf("Counting [%s][ K = %d ]", sampleName, k), stdout())
	flush.console()

	if(is.null(offset))
	{
	  offset = -1;
	}
  
	
	mmOption = .jnew("base/MarkovModelOption")
  
	result = J("main/SELEX")$doMinimalCounting( sample ,selex.getInt(k), mmOption, 
         selex.getBoolean(FALSE) ,selex.getInt(offset), selex.getBoolean(FALSE), seqfilter);
	
	countStats= J("main/SELEX")$getKMerCountStats(sample, selex.getInt(k),selex.getInt(offset),seqfilter)
	minCountResult = countStats$getLowestCount()
	write(sprintf("[ Lowest Count =  %s ]", minCountResult), stdout())
	
	flush.console()
	
	J("main/SELEX")$saveStats()
	
	selex.loadCountTable(sample, k, offset, minCount, top, numSort, markovModel=markovModel, seqfilter, outputPath)
	
}


#done #checked 
#Draw diagram with seqLogo: seqLogo(makePWM(t(mat)))
selex.kmerPSFM <- function(sample, k, offset=NULL)
{
  if(is.null(offset))
  {
    offset=-1;
  }
  stats = J("main/SELEX")$getPSFM(sample, selex.getInt(k), selex.getInt(offset))
  p1=sapply(stats[1],.jevalArray)
  p2=sapply(stats[2],.jevalArray)
  p3=sapply(stats[3],.jevalArray)
  p4=sapply(stats[4],.jevalArray)
  
  mat = data.frame(p1, p2, p3, p4, stringsAsFactors=FALSE)
  colnames(mat)=c("A","C","G","T")
  
  rm(stats)
  mat
}


#done #checked 
#Draw diagram with seqLogo: seqLogo(makePWM(t(mat)))
selex.samplePSFM <- function(sample)
{
  stats = J("main/SELEX")$getSamplePSFM(sample)
  p1=sapply(stats[1],.jevalArray)
  p2=sapply(stats[2],.jevalArray)
  p3=sapply(stats[3],.jevalArray)
  p4=sapply(stats[4],.jevalArray)
  
  mat = data.frame(p1, p2, p3, p4, stringsAsFactors=FALSE)
  colnames(mat)=c("A","C","G","T")
  
  rm(stats)
  mat
}

#done #checked 
#Draw diagram with seqLogo: seqLogo(makePWM(t(mat)))
selex.fastqPSFM <- function(seqName)
{
  #stats = J("main/SELEX")$getFastqPSFM(sample)
  stats = J("main/SELEX")$getFastqPSFM(seqName)
  p1=sapply(stats[1],.jevalArray)
  p2=sapply(stats[2],.jevalArray)
  p3=sapply(stats[3],.jevalArray)
  p4=sapply(stats[4],.jevalArray)
  
  mat = data.frame(p1, p2, p3, p4, stringsAsFactors=FALSE)
  colnames(mat)=c("A","C","G","T")
  
  rm(stats)
  mat
}



#done #checked
selex.loadCountTable <- function(sample, k, offset=NULL, minCount=100, top=-1, sort=FALSE, 
                             markovModel=NULL, seqfilter=NULL, outputPath = "")
{
  
  if(is.null(seqfilter))
  {
    seqfilter=selex.seqfilter()
  }
	useMarkov = TRUE
	if(is.null(markovModel))
	{
		useMarkov = FALSE
		markovModel = .jnew("base/MarkovModelInfo")
	}
  
  if(is.null(offset))
  {
    offset=-1;
  }
  if(top!=-1) #overwriting minCount
  {
    minCount=1;
  }
  
	stats = J("main/SELEX")$getKmerCount(sample, selex.getInt(k),selex.getInt(offset),
	               selex.getInt(minCount), selex.getInt(top), markovModel, seqfilter, outputPath)
	
  if(outputPath !="")
  {
    write(sprintf("Counts saved in file: %s", outputPath ), stdout())
    
    flush.console()
  }
  else
  {
    seqs=sapply(stats[1],.jevalArray)
    counts=sapply(stats[2],.jevalArray)
    if(useMarkov)
    {
      probs=sapply(stats[3],.jevalArray)
      exps=sapply(stats[4],.jevalArray)
      head(exps)
      mat = data.frame(seqs, counts, probs, exps, stringsAsFactors=FALSE)
      colnames(mat)=c("Kmer","ObservedCount","Probability","ExpectedCount")
    } else {
      mat = data.frame(seqs, counts, stringsAsFactors=FALSE)
      colnames(mat)=c("Kmer","ObservedCount")
    }
    if(sort)
    {
      mat = mat[order(-mat[,2]), ]
      rownames(mat) = NULL
    }
    rm(stats)
    mat
  }
	
}



#done #checked
selex.mmProb <- function(seqStr, markovModel)
{
  J("main/SELEX")$getProbability(seqStr, markovModel)
}


#done
selex.infogainSummary <- function(sample=NULL, displayFilter=FALSE)
{
	stats = J("main/SELEX")$getInformationGain()
	samples=sapply(stats[1],.jevalArray)
	ks=sapply(stats[2],.jevalArray)
	rs=sapply(stats[3],.jevalArray)
	markov=sapply(stats[4],.jevalArray)
	filters=sapply(stats[5],.jevalArray)
	methods=sapply(stats[6],.jevalArray)
	mat = data.frame(samples, ks, rs, markov,methods, filters, stringsAsFactors=FALSE)
	colnames(mat)=c("Sample","K","InformationGain","MarkovModel","MarkovModelType","Filter")
	mat = mat[order(mat[,1],mat[,5], mat[,2]), ]
	rownames(mat) = NULL
  
  if(!is.null(sample))
  {
    mat = mat[mat[,1] ==  J("main/SELEX")$getSampleID(sample) ,]
  }
	rm(stats)
	if(displayFilter==FALSE){
    mat = mat[,1:ncol(mat)-1]
	}
	return(mat)
}

#done
selex.mmSummary <- function(sample=NULL, displayFilter=FALSE)
{
	stats = J("main/SELEX")$getMarkovPR()
	samples=sapply(stats[1],.jevalArray)
	ks=sapply(stats[2],.jevalArray)
	ks=ks-1
	rs=sapply(stats[3],.jevalArray)
	cvS=sapply(stats[4],.jevalArray)
	cvL=sapply(stats[5],.jevalArray)
	filterL=sapply(stats[6],.jevalArray)
	methods=sapply(stats[7],.jevalArray)
  
	mat = data.frame(samples, ks, methods, rs,cvS, cvL, filterL, stringsAsFactors=FALSE)
	colnames(mat)=c("Sample","Order","MarkovModelType", "R","CVSample","CVLength", "Filter")
	mat = mat[order(mat[,1], mat[,3], mat[,2]), ]
	rownames(mat) = NULL
	
	if(!is.null(sample))
	{
	  mat = mat[mat[,1] ==  J("main/SELEX")$getSampleID(sample) ,]
	}
	rm(stats)
	if(displayFilter==FALSE){
	  mat = mat[,1:ncol(mat)-1]
	}
	return(mat)
}

#done
selex.countSummary <- function(sample=NULL, displayFilter=FALSE)
{
	stats = J("main/SELEX")$getCountStats()
	samples=sapply(stats[1],.jevalArray)
	ks=		  sapply(stats[2],.jevalArray)
	offsets=	sapply(stats[3],.jevalArray)
	offsets[offsets==-1] = NA
	mincounts= 	sapply(stats[4],.jevalArray)
	maxcounts= 	sapply(stats[5],.jevalArray)
	nums= 	sapply(stats[6],.jevalArray)
	filters= 	sapply(stats[7],.jevalArray)
	mat = data.frame(samples, ks, offsets, mincounts, maxcounts, nums,filters, stringsAsFactors=FALSE)
	colnames(mat)=c("Sample","K","offset","MinCount","MaxCount","TotalCounts","Filter")
	mat = mat[order(mat[,1],mat[,2],mat[,3]), ]
	rownames(mat) = NULL
	
	if(!is.null(sample))
	{
    sampleID=J("main/SELEX")$getSampleID(sample)
    print(sampleID) 
	  mat = mat[mat[,1] == sampleID  ,]
	}
	rm(stats)
	if(displayFilter==FALSE){
    mat = mat[,1:ncol(mat)-1]
	}
	return(mat)
}

#done
selex.summary <- function(sample=NULL, displayFilter=FALSE)
{
	informationGain = selex.infogainSummary(sample, displayFilter)
	markovModel     = selex.mmSummary(sample, displayFilter)
	countSummary    = selex.countSummary(sample, displayFilter)
	
	if(is.null(sample)){
		list(informationGain = informationGain,
				 markovModel   = markovModel,
				 countSummary  = countSummary
				 )	
	}else{
    list(informationGain = informationGain,
				 markovModel   = markovModel,
				 countSummary  = countSummary
				 )
	}
}

#done
selex.run <- function(trainingSample, crossValidationSample, minCount=100, 
                      infoGainSample, infoRange=NULL, mmMethod='DIVISION', mmWithLeftFlank = FALSE)
{
	k = selex.kmax(sample=crossValidationSample,threshold =minCount)
	mmOpt =   selex.mm( trainingSample, order=NA,  crossValidationSample = crossValidationSample, 
                            mmMethod=mmMethod, mmWithLeftFlank= mmWithLeftFlank)
	IG =  selex.infogain(infoGainSample, infoRange, mmOpt)
	print(selex.mmSummary(trainingSample))
	print(selex.infogainSummary(infoGainSample))
}

#
selex.jvmStatus <- function()
{
	stats = J("base/Util")$getMemeoryUsage()
  
	a=sapply(stats[1],.jevalArray)
	b=sapply(stats[2],.jevalArray)
	mat = data.frame(a,b, stringsAsFactors=FALSE)
	colnames(mat)=c("JVM","Status")
	rownames(mat) = NULL
	
	rm(stats)
  mat
}

selex.seqfilter <- function( 
  variableRegionIncludeRegex=NULL,
  variableRegionExcludeRegex=NULL,
  variableRegionGroupRegex=NULL,
  kmerIncludeRegex=NULL,
  kmerExcludeRegex=NULL,
  kmerIncludeOnly= NULL,  #c('AAAAA','AACCT')   #5. if present, select k-mers that appears in the list
  #after k-mer counting options
  viewIncludeRegex=NULL,
  viewExcludeRegex=NULL,
  viewIncludeOnly= NULL  #c('AAAAA','AACCT')   #8. shows only k-mer that in the list 
)
{
  seqfilter = .jnew("base/RegexOption")
  if(!is.null(variableRegionIncludeRegex))
  {
    seqfilter$setVariableRegionIncludeRegex(variableRegionIncludeRegex)
  }
  if(!is.null(variableRegionExcludeRegex))
  {
    seqfilter$setVariableRegionExcludeRegex(variableRegionExcludeRegex)
  }
  if(!is.null(variableRegionGroupRegex))
  {
    seqfilter$setVariableRegionGroupRegex(variableRegionGroupRegex)
  }
  
  if(!is.null(kmerIncludeRegex))
  {
    seqfilter$setKmerIncludeRegex(kmerIncludeRegex)
  }
  if(!is.null(kmerExcludeRegex))
  {
    seqfilter$setKmerExcludeRegex(kmerExcludeRegex)
  }
  if(!is.null(kmerIncludeOnly))
  {
    for(kmer in kmerIncludeOnly) {
      seqfilter$getKmerIncludeOnly()$add(kmer)
    }
  }
  
  if(!is.null(viewIncludeRegex))
  {
    seqfilter$setViewIncludeRegex(viewIncludeRegex)
  }
  if(!is.null(viewExcludeRegex))
  {
    seqfilter$setViewExcludeRegex(viewExcludeRegex)
  }
  if(!is.null(viewIncludeOnly))
  {
    for(kmer in viewIncludeOnly) {
      seqfilter$getViewIncludeOnly()$add(kmer)
    }
  }

  seqfilter
}

selex.getSeqfilter<- function(regex)
{
  regex$toString()
}

selex.exampledata<- function(outputFolder)
{
  files = c("R0.fastq.gz", "R2.fastq.gz", "config.xml")
  
  fullNames = c()
  for(f1 in files)
  { 
    srcf1 = system.file("extdata", f1, package="SELEX")
    status = file.copy(srcf1, outputFolder, overwrite = TRUE)
    write(sprintf("Extracting example data file to : %s/%s [%s]", outputFolder, f1, status ), stdout())
    fullNames = c(fullNames, sprintf("%s/%s", outputFolder, f1))
  }
  fullNames
}

