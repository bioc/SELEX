package test;

import main.CountFileReader;

public class CountFileReaderTest
{

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception
	{
		//test1();
		test2();
	}
	

	public static void test2() throws Exception
	{
		/*
		String filePath="/Users/homliu/Documents/workspace/SELEX/tmp/backup/exdUbx.exdScr.0.barcodeCCACGTC.v1.0.15.dat_D41D8CD98F00B204E9800998ECF8427E";
		String outputfilePath="/Users/homliu/Documents/workspace/SELEX/15_1.txt";
		CountFileReader.main(new String[]{filePath, outputfilePath});
		*/
		String filePath="/Users/homliu/Documents/workspace/SELEX/tmp/exdUbx.exdScr.L.2.barcodeCCACGTC.v1.low.2.14.dat";
		String outputfilePath="/Users/homliu/Documents/workspace/SELEX/15_2.txt";
		CountFileReader.main(new String[]{filePath});
	}
	
	public static void test1() throws Exception
	{
		String filePath="/Users/homliu/Documents/workspace/SELEX/tmp/training.8.dat";
		Long totalReads = 62760474L;
		int k=4;
		String folder = "/Users/homliu/Documents/workspace/SELEX/tmp/training.8.split";
		CountFileReader.split(filePath, totalReads, k, folder);
	}

}
