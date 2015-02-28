package test;

import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Properties;
import java.util.TreeSet;

import main.SELEX;
import main.SimpleKmerCount;
import base.ArraysMergerHeap;
import base.CountObject;
import base.DebugLog;
import base.MarkovModelInfo;
import base.MarkovModelOption;
import base.RegexOption;
import base.SELEXConfigReader;
import base.Sequence;
import config.ExperimentReference;
import config.KmerCountStats;

public class SELEXTest2
{

	
	public static void main(String[] args)
	{
		//testKmerCountsWithErrorOffset();
		
		//testKmerCounts2();
		//testKmerCounts3();
		//testKmerCounts3();
		
		//testSetConfig();System.exit(0);
		
		//testKmerCountsReversedComplementAndProp();
		//testKmerCountsReversedComplementAndMarkov();
		
		//testSplitting();
		//testSplitting2();System.exit(0);
		
		//testSplitting3_Markov();
		
		//testGetStatus();
		
		//testKmerCounts();
		//testKmerCounts2();
		
		//testKmerCountsWithFilters();
		//testKmerCountsWithKmerFilters();
		//testKmerCountsWithGroupFilters();
		//testKmerCountsWithViewFilters();

		//testKmerCounts20_1();
		//testKmerCounts20_2();
		
		
		//testKmerCountsPSFM();
		//testFastqPSFM();
		
		//testMarkovWithFilters();
		
		testFASTAKmerCounts();
	}
	
	public static void testGetStatus()
	{
		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);

		DebugLog.log("========Samples==========");
		Object[] newDS =SELEX.showSamples();
		for(Object ds:newDS)
		{
			DebugLog.log("==================");
			for(int i=0;i<Array.getLength(ds);i++)
			{
				DebugLog.log(Array.get(ds, i));
			}
		}
		

		DebugLog.log("========getCountStats==========");
		newDS =SELEX.getCountStats();
		for(Object ds:newDS)
		{
			DebugLog.log("==================");
			for(int i=0;i<Array.getLength(ds);i++)
			{
				DebugLog.log(Array.get(ds, i));
			}
		}
		
	}
	
	public static void testSetConfig0()
	{
		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);
		DebugLog.log(Arrays.toString(SELEX.showSamples()));
		
	}
	
	public static void testSetConfig()
	{
		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);
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
		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);
		DebugLog.log(Arrays.toString(SELEX.showSamples()));
		
		ExperimentReference training= SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);
		ExperimentReference testing=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCAGCTG.v1", 0);
		ExperimentReference info=     SELEX.getExperimentReference("exdUbx.exdScr.L.2","barcodeCCACGTC.v1.low", 2);

		SELEX.kmax(testing);

		MarkovModelInfo mm=SELEX.trainMarkovModel(training, testing, 6,  8, null, MarkovModelOption.DIVISION);
		DebugLog.log(mm.getMarkovR2());
		
		double IG = SELEX.calculateInformationGain(info, 12, mm, null);
		DebugLog.log("IG:"+IG);
		

		String filter="TTTT";
		MarkovModelInfo mm2=SELEX.trainMarkovModel(training, testing, 6,  8, null, MarkovModelOption.DIVISION);
		DebugLog.log(mm2.getMarkovR2());
		
		double IG2 = SELEX.calculateInformationGain(info, 12, mm2, null);
		DebugLog.log("IG2:"+IG2);
		
		
	}
	
	private static void testFASTAKmerCounts()
	{
		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1-fasta.small.xml";
		SELEX.loadConfigFile(config_path);

		ExperimentReference dataSet=  SELEX.getExperimentReference("sequence1","sample2", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		obj.setTempFolder("./tmp/");
		int k=5;
		int offset = -1;
		int mincount = -1;
		int top = 100;
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		
		Object[] counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, null);
		DebugLog.log("=================COUNTS==No Top===============");
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
	
	
	private static void testKmerCounts()
	{
		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);

		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		obj.setTempFolder("./tmp/");
		int k=16;
		int offset = -1;
		int mincount = -1;
		int top = 100;
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		
		Object[] counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, null);
		DebugLog.log("=================COUNTS==No Top===============");
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
		
		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);

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
		
		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);

		ExperimentReference training= SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);
		ExperimentReference testing=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCAGCTG.v1", 0);
		
		int kmax = 8;
		MarkovModelInfo mm=SELEX.trainMarkovModel(training,testing, 6 , kmax, null, MarkovModelOption.DIVISION);
		DebugLog.log(mm.getMarkovR2());
		
	}
	
	
	private static void testKmerCounts2()
	{

		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);
		
		//
		//exdUbx.exdScr.0',   sampleName='barcodeCCACGTC.v1
		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		obj.setTempFolder("./tmp/");
		int k=16;
		int offset = -1;
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		
		int top = -1;
		int minCount = 3;
		
		Object[] counts = SELEX.getKmerCount(dataSet, k, offset, minCount, top, null, null);
		DebugLog.log("=================COUNTS==Top ===============");
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
		
		k = 15;
		minCount = 4;
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		counts = SELEX.getKmerCount(dataSet, k, offset, minCount, top, null, null);
		DebugLog.log("=================COUNTS==Top ===============");
		len=Array.getLength(counts[0]);

		for(int i=0;i<len;i++)
		{
			StringBuffer sb=new StringBuffer();
			for(Object ds:counts)
			{
				sb.append(Array.get(ds, i)+"\t");
			}
			DebugLog.log(sb.toString());
		}
		

		k = 13;
		top=20;
		minCount = -1;
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		counts = SELEX.getKmerCount(dataSet, k, offset, minCount, top, null, null);
		DebugLog.log("=================COUNTS==Top ===============");
		len=Array.getLength(counts[0]);

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
	
	

	private static void testKmerCounts3()
	{

		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);
		
		//
		//exdUbx.exdScr.0',   sampleName='barcodeCCACGTC.v1
		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		obj.setTempFolder("./tmp/");
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

		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);
		
		//
		//exdUbx.exdScr.0',   sampleName='barcodeCCACGTC.v1
		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		obj.setTempFolder("./tmp/");
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
	
	private static void testKmerCountsWithErrorOffset()
	{
		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);

		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		obj.setTempFolder("./tmp/");
		int k=8;
		int offset = 8;
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, null);
		
		Object[] counts = null;
		int len=0;
		
		counts = SELEX.getKmerCount(dataSet, k, offset, 0, 20, null, null);
		DebugLog.log("=================COUNTS==Top ===============");
		len=Array.getLength(counts[0]);

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
	
	private static void testKmerCountsWithFilters()
	{

		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);
		
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

		RegexOption regeOption =new RegexOption();
		regeOption.setVariableRegionIncludeRegex("^A.{14}G");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, 0, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		len=Array.getLength(counts[0]);

		for(int i=0;i<len;i++)
		{
			StringBuffer sb=new StringBuffer();
			for(Object ds:counts)
			{
				sb.append(Array.get(ds, i)+"\t");
			}
			DebugLog.log(sb.toString());
		}
		

		regeOption =new RegexOption();
		regeOption.setVariableRegionExcludeRegex("^A.{14}G");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, 0, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		len=Array.getLength(counts[0]);

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
	
	private static void testKmerCountsWithKmerFilters()
	{

		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);
		
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

		RegexOption regeOption =new RegexOption();
		regeOption.setKmerIncludeRegex("^A.{14}G");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		len=Array.getLength(counts[0]);

		for(int i=0;i<len;i++)
		{
			StringBuffer sb=new StringBuffer();
			for(Object ds:counts)
			{
				sb.append(Array.get(ds, i)+"\t");
			}
			DebugLog.log(sb.toString());
		}
		

		regeOption =new RegexOption();
		regeOption.setKmerExcludeRegex("^A.{14}G");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		len=Array.getLength(counts[0]);

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
	
	private static void testKmerCountsWithViewFilters()
	{

		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);
		
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

		RegexOption regeOption =new RegexOption();
		regeOption.setViewIncludeRegex("^A.{14}G");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		len=Array.getLength(counts[0]);

		for(int i=0;i<len;i++)
		{
			StringBuffer sb=new StringBuffer();
			for(Object ds:counts)
			{
				sb.append(Array.get(ds, i)+"\t");
			}
			DebugLog.log(sb.toString());
		}
		
		regeOption =new RegexOption();
		regeOption.setViewExcludeRegex("^A.{14}G");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, mincount, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		len=Array.getLength(counts[0]);

		for(int i=0;i<len;i++)
		{
			StringBuffer sb=new StringBuffer();
			for(Object ds:counts)
			{
				sb.append(Array.get(ds, i)+"\t");
			}
			DebugLog.log(sb.toString());
		}
		
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
		len=Array.getLength(counts[0]);

		for(int i=0;i<len;i++)
		{
			StringBuffer sb=new StringBuffer();
			for(Object ds:counts)
			{
				sb.append(Array.get(ds, i)+"\t");
			}
			DebugLog.log(sb.toString());
		}
		
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
		len=Array.getLength(counts[0]);

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
	
	private static void testKmerCounts20_1()
	{

		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq100.xml";
		SELEX.loadConfigFile(config_path);
		
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
	
	private static void testKmerCounts20_2()
	{

		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq101.xml";
		SELEX.loadConfigFile(config_path);
		
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
	
	private static void testKmerCountsWithGroupFilters()
	{

		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);
		
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
		int len = 0;
		
		DebugLog.log("=================COUNTS==Top ===============");

		regeOption =new RegexOption();
		k=13;
		regeOption.setVariableRegionGroupRegex("^TGT([ACGT]{13})");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, 0, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		len=Array.getLength(counts[0]);

		for(int i=0;i<len;i++)
		{
			StringBuffer sb=new StringBuffer();
			for(Object ds:counts)
			{
				sb.append(Array.get(ds, i)+"\t");
			}
			DebugLog.log(sb.toString());
		}
		

		regeOption =new RegexOption();
		k=5;
		regeOption.setVariableRegionGroupRegex("^TGT([ACGT]{5})");
		SELEX.doMinimalCounting(dataSet, k,  obj, null, false, offset,false,-1, regeOption);
		
		counts = SELEX.getKmerCount(dataSet, k, offset, 0, top, null, regeOption);
		DebugLog.log("=================COUNTS==Top ===============");
		len=Array.getLength(counts[0]);

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
	
	private static void testKmerCountsPSFM()
	{
		String tempFolder = "./tmp";
		SELEX.setWorkingDirectory(tempFolder);
		SELEX.addSequenceInfo("seqX","/home/hamburger/workspace/SELEX/fastq_small1","sampleX",  0,16, "TGG","CCAGCTG",null,null);
		//GAGGTCGATCTACCTT
		ExperimentReference dataSet=  SELEX.getExperimentReference("seqX","sampleX", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		obj.setTempFolder(tempFolder);
		SELEX.doMinimalCounting(dataSet, 2,  obj, null, false, null,true,-1, null);
		
		Object[] counts = SELEX.getKmerCount(dataSet, 2, -1, 1, -1, null,null);
		DebugLog.log("=================COUNTS=================");
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
		
		Object[] psfm = SELEX.getPSFM( dataSet, 2,-1);

		DebugLog.log("=================PSFM=================");
		len=Array.getLength(psfm[0]);

		for(int i=0;i<len;i++)
		{
			StringBuffer sb=new StringBuffer();
			for(Object ds:psfm)
			{
				sb.append(Array.get(ds, i)+"\t");
			}
			DebugLog.log(sb.toString());
		}
		
	}
	
	private static void testFastqPSFM()
	{
		String tempFolder = "./tmp";
		SELEX.setWorkingDirectory(tempFolder);
		SELEX.addSequenceInfo("seqX","/home/hamburger/workspace/SELEX/fastq_small1","sampleX",  0,16, "TGG","CCAGCTG",null,null);
		//GAGGTCGATCTACCTT
		ExperimentReference dataSet=  SELEX.getExperimentReference("seqX","sampleX", 0);

		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		obj.setTempFolder(tempFolder);
		SELEX.doMinimalCounting(dataSet, 2,  obj, null, false, null,true,-1, null);
		
		Object[] psfm = SELEX.getFastqPSFM(dataSet);

		DebugLog.log("=================PSFM=================");
		int len=Array.getLength(psfm[0]);

		for(int i=0;i<len;i++)
		{
			StringBuffer sb=new StringBuffer();
			for(Object ds:psfm)
			{
				sb.append(Array.get(ds, i)+"\t");
			}
			DebugLog.log(sb.toString());
		}
		
	}
	
	private static void testSplitting()
	{

		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);
		
		//
		SimpleKmerCount obj=new SimpleKmerCount();

		ExperimentReference dataSet=  SELEX.getExperimentReference("exdUbx.exdScr.0","barcodeCCACGTC.v1", 0);
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		ArrayList<Double> ratios = new  ArrayList<Double>();
		ratios.add(0.5);
		ratios.add(0.5);
		Object[] newDS = SELEX.splitDataSet( dataSet, "xxx", ratios);
		for(Object ds:newDS)
		{
			for(int i=0;i<Array.getLength(ds);i++)
			{
				DebugLog.log(Array.get(ds, i));
			}
		}
		
		
		DebugLog.log("=================SAMPLES=================");
		Object[] samples =SELEX.showSamples();
		for(Object ds:samples)
		{
			for(int i=0;i<Array.getLength(ds);i++)
			{
				DebugLog.log(Array.get(ds, i));
			}
		}
	}
	
	public static void testSplitting2()
	{

		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);

		DebugLog.log("=================SAMPLES=================");
		Object[] samples =SELEX.showSamples();
		for(Object ds:samples)
		{
			for(int i=0;i<Array.getLength(ds);i++)
			{
				DebugLog.log(Array.get(ds, i));
			}
		}
		
		ExperimentReference training= SELEX.getExperimentReference("exdUbx.exdScr.0.xxx.1","barcodeCCACGTC.v1.xxx.1", 0);
		ExperimentReference testing=  SELEX.getExperimentReference("exdUbx.exdScr.0.xxx.2","barcodeCCACGTC.v1.xxx.2", 0);
		
		int k = SELEX.kmax(testing);

		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), training);
		
		MarkovModelInfo mm=SELEX.trainMarkovModel(training, testing, 5 ,  k, null, null);
		DebugLog.log(mm.getMarkovR2());
		
	}
	

	public static void testSplitting3_Markov()
	{

		SELEX.setWorkingDirectory("./tmp/");
		
		String 
		config_path="/home/hamburger/workspace/SELEX/config-Seq1.xml";
		SELEX.loadConfigFile(config_path);
		config_path="/home/hamburger/workspace/SELEX/config-Seq2.xml";
		SELEX.loadConfigFile(config_path);

		DebugLog.log("=================SAMPLES=================");
		Object[] samples =SELEX.showSamples();
		for(Object ds:samples)
		{
			for(int i=0;i<Array.getLength(ds);i++)
			{
				DebugLog.log(Array.get(ds, i));
			}
		}
		
		ExperimentReference training= SELEX.getExperimentReference("exdUbx.exdScr.0.xxx.1","barcodeCCACGTC.v1.xxx.1", 0);
		ExperimentReference testing=  SELEX.getExperimentReference("exdUbx.exdScr.0.xxx.2","barcodeCCACGTC.v1.xxx.2", 0);
		
		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), training);

		obj.setTempFolder("./tmp/");
		SELEX.doMinimalCounting(training, 2,  obj,  null, false,-1, false, -1, null);
		
//		
//		MarkovModelInfo mm=SELEX.trainMarkovModel(training, 5,  obj);
//		DebugLog.log(mm.getMarkovR2());
//		
//		obj.initTraining(SELEX.getConfigReader(), testing);
//		double IG = SELEX.calculateInformationGain(testing, 9, obj, mm);
//		DebugLog.log(IG);
		
		
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
