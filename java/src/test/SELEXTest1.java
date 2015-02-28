package test;

import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;

import main.SELEX;
import main.SimpleKmerCount;
import base.ArraysMergerHeap;
import base.CountObject;
import base.DebugLog;
import base.MarkovModelInfo;
import base.MarkovModelOption;
import base.RegexOption;
import base.Sequence;
import config.ExperimentReference;

public class SELEXTest1
{
	/**
	 * Setup the configuration folder , data folder and temp folder.
	 * 
	 * sample_data and sample_config can be found under SELEX/dahong/ on github.
	 */
	public static String dataFolder = "./sample_data/";
	public static String configFolder = "./sample_config/";
	public static String tempFolder = "./tmp/";
	
	public static void main(String[] args)
	{	
		/*
		 * Uncomment one of these test cases to run 
		 */
		
		testKmerCounts();			//k-mer counting test
		//testKmerCounts2();		//more k-mer counting tests
		//testKmerCounts3();
		//testKmerCounts4();		//exception example
		//testKmerCountsWithOffset();	//k-mer counts with offset setting
		//testKmerCountsMplex();

		//k-mer counting with various types of filters
		//testKmerCountsWithFilters();		
		//testKmerCountsWithKmerFilters();
		//testKmerCountsWithGroupFilters();
		//testKmerCountsWithViewFilters();

		//k-mer counting with k =  20
		//testKmerCounts20_1();
		//testKmerCounts20_2();
		
		//testGetRound0();          //Test getRound0 and k-mer counts
		//testMarkov();				//Markov model generation
		//testMarkovWithFilters();	//Markov model generation with filters
		//testGetStatus();			//show sample info and kmer count stats
		
		//testSplitting();			//sample splitting and testing
				
		//testKmerCountsPSFM();		//PFSM on k-mers
		//testFastqPSFM();			//PSFM on fastq file

		//testGetProbability();		//get probability for a sequence
		
//		testKmerCountsReversedComplementAndProp();
//		testKmerCountsReversedComplementAndMarkov();

	}
	
	/**
	 * Test getRound0() function and k-mer counting.
	 * 
	 * The Round0 link can be setup from either the XML configuration file,
	 * 		or on the fly in the SELEX.addSequenceInfo() method.
	 */
	public static void testGetRound0()
	{
		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path,dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		
		SELEX.loadConfigFile(config_path,dataFolder);
		DebugLog.log(Arrays.toString(SELEX.showSamples()));

		SELEX.addSequenceInfo("seqX",configFolder + "/fastq_small1",
							  "sampleX",  0, 16, "TGG","CCAGCTG","XXXXXXXX","ZZZZ",
							  "exdUbx.exdScr.0", "barcodeCCACGTC.v1");
		//GAGGTCGATCTACCTT
		
		SELEX.exportConfigFile(configFolder + "/config-Seq99.xml");
		
		ExperimentReference info=     
			SELEX.getExperimentReference("exdUbx.exdScr.L.2","barcodeCCACGTC.v1.low", 2);
		ExperimentReference info_round0=  
			SELEX.getRound0(info);
		
		if(info_round0!=null)
		{
			int mincount = 0;
			int top = 10;
			int k=10;
			SELEX.doMinimalCounting(info_round0, k,  null, null, false, -1,false,-1, null);
			Object[] counts = SELEX.getKmerCount(info_round0, k, -1, mincount, top, null, null);
			
			DebugLog.log("=================COUNTS== "+SELEX.getSampleID(info_round0)+" ===============");
			
			print2DArray(counts);
		}
		
		ExperimentReference sampleX= SELEX.getExperimentReference("seqX","sampleX", 0);
		ExperimentReference sampleX_round0=  SELEX.getRound0(sampleX);
		if(sampleX_round0!=null)
		{
			int mincount = 0;
			int top = 10;
			int k=10;
			SELEX.doMinimalCounting(sampleX_round0, k,  null, null, false, -1,false,-1, null);;
			Object[] counts = SELEX.getKmerCount(sampleX_round0, k, -1, mincount, top, null, null);
			DebugLog.log("=================COUNTS== "+SELEX.getSampleID(sampleX_round0)+" ===============");

			print2DArray(counts);
		}
	}

	/**
	 * Tests the Markov model generation with and without the left flank
	 * Uses them to get expected counts.
	 */
	public static void testMarkov()
	{
		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		DebugLog.log(Arrays.toString(SELEX.showSamples()));
		
		ExperimentReference training= SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);
		ExperimentReference testing=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCAGCTG.v1", 0);
		ExperimentReference info=     SELEX.getExperimentReference("exdUbx.exdScr.L.2","barcodeCCACGTC.v1.low", 2);

		//SELEX.kmax(testing);
		
		int kmax = 8;

		MarkovModelInfo mm=SELEX.trainMarkovModel(training, testing, 6,  kmax, null, 
				MarkovModelOption.DIVISION+"|"+MarkovModelOption.WITH_LEFT_FLANK);
		
		DebugLog.log(mm.getMarkovR2());
		
		double IG = SELEX.calculateInformationGain(info, 12, mm, null);
		DebugLog.log("IG:"+IG);
		
		int mincount = 0;
		int top = 15;
		Object[] counts = SELEX.getKmerCount(testing, kmax, -1, mincount, top, mm, null);
		DebugLog.log("=================COUNTS== Top ===============");

		print2DArray(counts);

		// alternative: Markov model without left flank information
		mm=SELEX.trainMarkovModel(training, testing, 6,  kmax, null, 
				MarkovModelOption.DIVISION);
		DebugLog.log(mm.getMarkovR2());
		
		IG = SELEX.calculateInformationGain(info, 12, mm, null);
		DebugLog.log("IG:"+IG);
		
		counts = SELEX.getKmerCount(testing, kmax, -1, mincount, top, mm, null);
		DebugLog.log("=================COUNTS== Top ===============");

		print2DArray(counts);

	}
	
	/**
	 * Loads samples and displays the info.
	 * 
	 * Shows Kmer count stats.
	 */
	public static void testGetStatus()
	{
		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);

		DebugLog.log("========Samples==========");
		Object[] newDS =SELEX.showSamples();
		print2DArray(newDS);

		DebugLog.log("========getCountStats==========");
		newDS =SELEX.getCountStats();
		print2DArray(newDS);
		
	}
	
	public static void testSetConfig0()
	{
		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		DebugLog.log(Arrays.toString(SELEX.showSamples()));
		
	}
	
	public static void testGetProbability()
	{
		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		DebugLog.log(Arrays.toString(SELEX.showSamples()));
		
		ExperimentReference training= SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);
		ExperimentReference testing=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCAGCTG.v1", 0);
		ExperimentReference info=     SELEX.getExperimentReference("exdUbx.exdScr.L.2","barcodeCCACGTC.v1.low", 2);

		SELEX.kmax(testing);

		MarkovModelInfo mm=SELEX.trainMarkovModel(training, testing, 6,  8, null, MarkovModelOption.DIVISION);
		DebugLog.log(mm.getMarkovR2());
		
		double IG = SELEX.calculateInformationGain(info, 12, mm, null);
		DebugLog.log(IG);
		
		System.out.println(SELEX.getProbability("AAAAAAAAAA", mm));
	}
	
	
	
	public static void testMarkovWithFilters()
	{
		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		DebugLog.log(Arrays.toString(SELEX.showSamples()));
		
		ExperimentReference training= SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);
		ExperimentReference testing=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCAGCTG.v1", 0);
		ExperimentReference info=     SELEX.getExperimentReference("exdUbx.exdScr.L.2","barcodeCCACGTC.v1.low", 2);

		SELEX.kmax(testing);

		MarkovModelInfo mm=SELEX.trainMarkovModel(training, testing, 6,  8, null, MarkovModelOption.DIVISION);
		DebugLog.log("================> R2:"+mm.getMarkovR2());
		
		double IG = SELEX.calculateInformationGain(info, 12, mm, null);
		DebugLog.log("================> IG:"+IG);

		RegexOption regexOption = new RegexOption();
		regexOption.setVariableRegionExcludeRegex("^TTT");
		MarkovModelInfo mm2=SELEX.trainMarkovModel(training, testing, 6,  8, regexOption, MarkovModelOption.DIVISION);
		DebugLog.log("================> R2:"+mm2.getMarkovR2());
		
		double IG2 = SELEX.calculateInformationGain(info, 12, mm2, null);
		DebugLog.log("================> IG2:"+IG2);
		
		
	}
	
	/**
	 * K-mer counting test and retrieves the top result.
	 */
	private static void testKmerCounts()
	{
		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);

		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		int k=16;
		int offset = -1;
		int mincount = -1;
		int top = 10;
		SELEX.doMinimalCounting(dataSet, k,  null, null, false, offset,false,-1, null);
		
		Object[] counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, null);
		DebugLog.log("=================COUNTS==Top===============");
		print2DArray(counts);
		
	}
	

	private static void testKmerCountsMplex()
	{
		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.mplex.xml";
		SELEX.loadConfigFile(config_path, dataFolder);

		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);
																							 
		int k=8;
		int offset = -1;
		int mincount = -1;
		int top = 10;
		SELEX.doMinimalCounting(dataSet, k,  null, null, false, offset,false,-1, null);
		
		Object[] counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, null);
		DebugLog.log("=================COUNTS==Top===============");
		print2DArray(counts);
		
		
	}
	
	private static void print2DArray(Object[] counts)
	{
		int len=Array.getLength(counts[0]);

		for(int i=0;i<len;i++)
		{
			StringBuffer sb=new StringBuffer();
			for(Object ds:counts)
			{
				sb.append(Array.get(ds, i)+"\t");
			}
			DebugLog.log(sb.toString());
		}
	}
	
	private static void testKmerCountsReversedComplementAndProp()
	{

		SimpleKmerCount.includeReversed = true;
		
		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);

		ExperimentReference training= SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);
		
		MarkovModelInfo mm=SELEX.trainMarkovModel(training, 6 , null, MarkovModelOption.DIVISION);
		DebugLog.log(mm.getMarkovR2());
		
		float[] marKovProbabilities = null;
		String morkovObjPath = null;
		String markovCountPath = null;
		int markovModelLength = 0;
		int[] markovCounts = null;
		long markovTotalCount = 0;
		try
		{
			morkovObjPath = mm.getMarkovObjPath();
			markovCountPath = mm.getMarkovCountsPath();
			markovModelLength = mm.getMarkovLength();
			markovTotalCount = mm.getMarkovLengthTotalCount();
			
			// load markovProbabilities
			DebugLog.log("Reading Markov Prob file:" + morkovObjPath);
			FileInputStream fos = new FileInputStream(morkovObjPath);
			ObjectInputStream oos = new ObjectInputStream(fos);
			marKovProbabilities = (float[]) oos.readObject();

			DebugLog.log("Reading Markov Count file:" + markovCountPath);
			ArraysMergerHeap.MinHeapNode node2 = new ArraysMergerHeap.MinHeapNode(
					markovCountPath);
			CountObject obj = null;
			markovCounts = new int[1 << (2 * (markovModelLength))];

			while ((obj = node2.peek()) != null)
			{
				Sequence seq = (Sequence) obj.getKey();
				markovCounts[(int) seq.getValue()] = obj.getCount();
				node2.pop();
			}
			
		}
		catch(Exception ex)
		{
			DebugLog.log(ex);
			throw new RuntimeException(ex);
		}
		
		long size =  (1L << 32) -1;
		DebugLog.log("size = "+ size);
		//size =  100;
		int len=16;
		double sum = 0;
		//PrintWriter p=new PrintWriter("16-mer-prob.txt");
		for(long i=0;i<size;i++)
		{
			Sequence seq=new Sequence(i, len);
			
			double modelProb = SimpleKmerCount.getPredictedCount(
					markovModelLength, markovTotalCount, seq,
					markovCounts, marKovProbabilities);

			sum += modelProb;
			
			//DebugLog.log(seq.getString() + " --> "+ modelProb);
			if(i%20000000 ==0)
			{
				DebugLog.log(i);
			}
		}
		DebugLog.log("sum = "+ sum);
		
	}
	

	private static void testKmerCountsReversedComplementAndMarkov()
	{

		SimpleKmerCount.includeReversed = true;
		
		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);

		ExperimentReference training= SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);
		ExperimentReference testing=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCAGCTG.v1", 0);
		
		int kmax = 8;
		MarkovModelInfo mm=SELEX.trainMarkovModel(training,testing, 6 , kmax, null, MarkovModelOption.DIVISION);
		DebugLog.log(mm.getMarkovR2());
		
	}
	
	
	private static void testKmerCounts2()
	{

		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		
		//
		//exdUbx.exdScr.0',   sampleName='barcodeCCACGTC.v1
		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		obj.setTempFolder(tempFolder);
		int k=16;
		int offset = -1;
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		
		int top = -1;
		int minCount = 3;
		
		Object[] counts = SELEX.getKmerCount(dataSet, k, offset, minCount, top, null, null);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);
		
		k = 15;
		minCount = 4;
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		counts = SELEX.getKmerCount(dataSet, k, offset, minCount, top, null, null);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);

		k = 13;
		top=20;
		minCount = -1;
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		counts = SELEX.getKmerCount(dataSet, k, offset, minCount, top, null, null);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);
	}
	
	/**
	 * Outputs the result to files.
	 */
	private static void testKmerCounts3()
	{

		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		
		//
		//exdUbx.exdScr.0',   sampleName='barcodeCCACGTC.v1
		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		obj.setTempFolder(tempFolder);
		int k=16;
		int offset = -1;
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		
		int top = -1;
		int minCount = 1;
		
		Object[] counts = SELEX.getKmerCount(dataSet, k, offset, minCount, top, null, null,"/tmp/bb1_list");
		
		top = 1000;
		minCount = 1;
		counts = SELEX.getKmerCount(dataSet, k, offset, minCount, top, null, null,"/tmp/bb2_list");
		
		top = -1;
		minCount = 3;
		counts = SELEX.getKmerCount(dataSet, k, offset, minCount, top, null, null,"/tmp/bb3_list");
		
	}
	
	private static void testKmerCounts4()
	{

		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		
		//
		//exdUbx.exdScr.0',   sampleName='barcodeCCACGTC.v1
		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		obj.setTempFolder(tempFolder);
		int k=20;
		int offset = -1;
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		
		int top = -1;
		int minCount = 1;
		
		Object[] counts = SELEX.getKmerCount(dataSet, k, offset, minCount, top, null, null,"/tmp/bb1_list");
		
		top = 1000;
		minCount = 1;
		counts = SELEX.getKmerCount(dataSet, k, offset, minCount, top, null, null,"/tmp/bb2_list");
		
		top = -1;
		minCount = 3;
		counts = SELEX.getKmerCount(dataSet, k, offset, minCount, top, null, null,"/tmp/bb3_list");
		
	}
	
	private static void testKmerCountsWithOffset()
	{
		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);

		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		obj.setTempFolder(tempFolder);
		int k=8;
		int offset = 8;
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		
		Object[] counts = null;
		counts = SELEX.getKmerCount(dataSet, k, offset, 0, 20, null, null);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);
	}
	
	/**
	 * K-mer counting with VariableRegion filters
	 */
	private static void testKmerCountsWithFilters()
	{

		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		
		//
		//exdUbx.exdScr.0',   sampleName='barcodeCCACGTC.v1
		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		int k=16;
		int offset = -1;
		int top = 30;
		
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		
		Object[] counts = SELEX.getKmerCount(dataSet, k, offset, 0, top, null, null);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);

		RegexOption regeOption =new RegexOption();
		regeOption.setVariableRegionIncludeRegex("^A.{14}G");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, 0, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);
		

		regeOption =new RegexOption();
		regeOption.setVariableRegionExcludeRegex("^A.{14}G");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, 0, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);
	}

	/**
	 * K-mer counting with Kmer filters
	 */
	private static void testKmerCountsWithKmerFilters()
	{

		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		
		//
		//exdUbx.exdScr.0',   sampleName='barcodeCCACGTC.v1
		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		int k=16;
		int offset = -1;
		int top = 100;
		int mincount = -1;
		
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		
		Object[] counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, null);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);

		RegexOption regeOption =new RegexOption();
		regeOption.setKmerIncludeRegex("^A.{14}G");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);
		

		regeOption =new RegexOption();
		regeOption.setKmerExcludeRegex("^A.{14}G");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);
	}
	
	
	/**
	 * Retrieves K-mer count results with view filters
	 */
	private static void testKmerCountsWithViewFilters()
	{

		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		
		//
		//exdUbx.exdScr.0',   sampleName='barcodeCCACGTC.v1
		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		int k=16;
		int offset = -1;
		int top = -1;
		int mincount = 3;
		
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		
		Object[] counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, null);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);

		RegexOption regeOption =new RegexOption();
		regeOption.setViewIncludeRegex("^A.{14}G");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);
		
		regeOption =new RegexOption();
		regeOption.setViewExcludeRegex("^A.{14}G");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);
		
		//////////////////

		regeOption =new RegexOption();
		regeOption.getKmerIncludeOnly().add("GATTAAATATTTTTTG");
		regeOption.getKmerIncludeOnly().add("GATGTAGTCTTTTTTT");
		regeOption.getKmerIncludeOnly().add("GATATTAAGCCTATGT");
		regeOption.getKmerIncludeOnly().add("CGGTTAATCAATTTGT");
		regeOption.getKmerIncludeOnly().add("ATTTTATATGTTTTTT");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, regeOption);
		DebugLog.log("=================COUNTS==Kmer selection Top ===============");
		print2DArray(counts);
		
		//////////////////////

		regeOption =new RegexOption();
		regeOption.getViewIncludeOnly().add("GATTAAATATTTTTTG");
		regeOption.getViewIncludeOnly().add("GATGTAGTCTTTTTTT");
		regeOption.getViewIncludeOnly().add("GATATTAAGCCTATGT");
		regeOption.getViewIncludeOnly().add("CGGTTAATCAATTTGT");
		regeOption.getViewIncludeOnly().add("ATTTTATATGTTTTTT");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);
	}
	
	private static void testKmerCounts20_1()
	{

		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq100.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		
		//
		//exdUbx.exdScr.0',   sampleName='barcodeCCACGTC.v1
		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.L.2.xx","xxx", 2);

		SimpleKmerCount obj=new SimpleKmerCount();
		int k=20;
		int offset = -1;
		int top = 100;
		int mincount = -1;
		
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		
		Object[] counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, null);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);
		
	}
	
	private static void testKmerCounts20_2()
	{

		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq101.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		
		//
		//exdUbx.exdScr.0',   sampleName='barcodeCCACGTC.v1
		ExperimentReference dataSet=  SELEX.getExperimentReference("2013-1897_35","R3_35", 3);

		SimpleKmerCount obj=new SimpleKmerCount();
		int k=16;
		int offset = -1;
		int top = 100;
		int mincount = -1;
		
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		
		Object[] counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, null);
		DebugLog.log("=================COUNTS==Top ===============");

		print2DArray(counts);
	}
	
	/**
	 * Test k-mer counting with group filters
	 */
	private static void testKmerCountsWithGroupFilters()
	{

		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		
		//
		//exdUbx.exdScr.0',   sampleName='barcodeCCACGTC.v1
		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		int k=16;
		int offset = -1;
		int top = 100;
		int mincount = 1;

		RegexOption regeOption =new RegexOption();
		regeOption.setVariableRegionIncludeRegex("^TGT");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);

		SELEX.getKmerCount(dataSet, k, offset, mincount, -1, null, regeOption, "/tmp/output");
		
		Object[] counts =null ;
		
		DebugLog.log("=================COUNTS==Top ===============");

		regeOption =new RegexOption();
		k=13;
		regeOption.setVariableRegionGroupRegex("^TGT([ACGT]{13})");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, 0, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);

		regeOption =new RegexOption();
		k=5;
		regeOption.setVariableRegionGroupRegex("^TGT([ACGT]{5})");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, 0, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		print2DArray(counts);
	}
	
	private static void testKmerCountsPSFM()
	{
		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		
		//
		//exdUbx.exdScr.0',   sampleName='barcodeCCACGTC.v1
		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		obj.setTempFolder(tempFolder);
		SELEX.doMinimalCounting(dataSet, 2,  obj, null, false, null,true,-1, null);
		
		Object[] counts = SELEX.getKmerCount(dataSet, 2, -1, 1, -1, null,null);
		DebugLog.log("=================COUNTS=================");
		print2DArray(counts);
		
		Object[] psfm = SELEX.getPSFM( dataSet, 2,-1);

		DebugLog.log("=================PSFM=================");
		print2DArray(psfm);
		
	}
	
	private static void testFastqPSFM()
	{

		SELEX.setWorkingDirectory(tempFolder);
		
		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path, dataFolder);
		
		//
		//exdUbx.exdScr.0',   sampleName='barcodeCCACGTC.v1
		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		obj.setTempFolder(tempFolder);
		SELEX.doMinimalCounting(dataSet, 2,  obj, null, false, null,true,-1, null);
		
		Object[] psfm = SELEX.getFastqPSFM(dataSet);

		DebugLog.log("=================FASTQ PSFM=================");
		print2DArray(psfm);
		

		psfm = SELEX.getSamplePSFM(dataSet);
		DebugLog.log("=================Sample PSFM=================");
		print2DArray(psfm);
		
	}
	
	private static void testSplitting()
	{

		SELEX.setWorkingDirectory(tempFolder);

		String 
		config_path=configFolder + "/config-Seq1.xml";
		SELEX.loadConfigFile(config_path,dataFolder);
		config_path=configFolder + "/config-Seq2.xml";
		SELEX.loadConfigFile(config_path,dataFolder);
		

		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);
		
		ArrayList<Double> ratios = new  ArrayList<Double>();
		ratios.add(0.5);
		ratios.add(0.5);
		Object[] newDS = SELEX.splitDataSet( dataSet, "test", ratios);

		print2DArray(newDS);
		
		DebugLog.log("=================SAMPLES=================");
		Object[] samples =SELEX.showSamples();

		print2DArray(samples);
		
		ExperimentReference training= SELEX.getExperimentReference("exdUbx.exdScr.0.test.1","barcodeCCACGTC.v1.test.1", 0);
		ExperimentReference testing=  SELEX.getExperimentReference("exdUbx.exdScr.0.test.2","barcodeCCACGTC.v1.test.2", 0);
		
		int k = SELEX.kmax(testing);
		
		MarkovModelInfo mm=SELEX.trainMarkovModel(training, testing, 5 ,  k, null, null);
		DebugLog.log(mm.getMarkovR2());
	}
	
//	private static void testGetMarkModelInfo()
//	{
//		Properties stats =loadStats(false) ;
//		MarkovModelInfo mm = getMarkModelInfo("training", 8 , stats);
//		System.out.println(mm);
//		Object[] objs= getKmerCount("testing", 8, 100, mm);
//		for(Object obj:objs)
//		{
//			System.out.println("=======================================");
//			int len= Array.getLength(obj);
//			for(int i=0;i<Math.min(len, 3);i++)
//			{
//				Object obj2= Array.get(obj, i);
//				System.out.println(obj2);
//			}
//
//		}
//	}
}
