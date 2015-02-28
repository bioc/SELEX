package test;

import base.Sequence;

import java.util.Arrays;
import java.util.Map;

public class SequenceTest
{

	public static void main(String args[])
	{
		//test3();
		
		//test1();
		//test2();
		//test("ttggccaattggccaattggccaattggccaat"); //wrong. too long
		
		test_reversed();
        //testGenerateFeaturesMatrix();
        //testFeatures();

	}
	
	public static void test4()
	{

		Sequence[] seqs =new Sequence[]{
				new Sequence("TC"),
				new Sequence("AA"),
				new Sequence("GG"),
				new Sequence("GC"),
				new Sequence("CG"),
				new Sequence("CC"),
				new Sequence("CT"),
				new Sequence("TA"),
				new Sequence("TT"),
				new Sequence("TG"),
				new Sequence("AG"),
				new Sequence("AC"),
				new Sequence("AT"),
				new Sequence("CA"),
				new Sequence("GT"),
				new Sequence("GA")
		};
		Arrays.sort(seqs);
		System.out.println(Arrays.toString(seqs));
		
	}
	public static void test3()
	{
		long temp = 454761243;
		Sequence test = new Sequence(temp,16);
		System.out.println(test.getString());
		//ACGTACGTACGTACGT

		test = new Sequence(temp,25);
		System.out.println(test.getString());
		///AAAAAAAAAACGTACGTACGTACGT


		Sequence test2 = new Sequence(temp, 4);
		System.out.println(
				test2.getString()
				);
		//CGTACGTACGTACGT
	}
	public static void test1()
	{
		test("aattggccaaaattggccaaaattggccaa");
		test("ttggccaa");
		test("ttggccaattggccaattggccaattggccaa");
	}
	
	public static void test2()
	{
		String seqStr = "aaaa";
		Sequence seq3 =new Sequence(seqStr, 0, seqStr.length());
		System.out.println("==============test2=================");
		System.out.println("value:"+seq3.getValue());
		System.out.println("string:"+seq3.getString());
		System.out.println("length:"+seq3.getLength());
	}
	
	public static void test(String seqStr)
	{
		//sequence --> number
		//String seqStr="aattggccaaaattggccaaaattggccaa";
		Sequence seq1= new Sequence(seqStr);
		System.out.println("==============1=================");
		System.out.println("value:"+seq1.getValue());
		System.out.println("string:"+seq1.getString());
		System.out.println("length:"+seq1.getLength());
		
		//number --> sequence
		Sequence seq2= new Sequence(seq1.getValue(),seqStr.length());
		System.out.println("==============2=================");
		System.out.println("value:"+seq2.getValue());
		System.out.println("string:"+seq2.getString());
		System.out.println("length:"+seq2.getLength());
		
		Sequence seq3 =new Sequence(seq1.getString(), 0, seqStr.length());
		System.out.println("==============3=================");
		System.out.println("value:"+seq3.getValue());
		System.out.println("string:"+seq3.getString());
		System.out.println("length:"+seq3.getLength());
	}
	

	public static void test_reversed()
	{
		String seqStr = "AAAACGT";
		Sequence seq3 =new Sequence(seqStr, 0, seqStr.length());
		System.out.println("==============test2=================");
		System.out.println("value:"+seq3.getValue());
		System.out.println("string:"+seq3.getString());
		System.out.println("length:"+seq3.getLength());

		Sequence seq4 =seq3.getReverseComplement();
		System.out.println("==============test2=================");
		System.out.println("value:"+seq4.getValue());
		System.out.println("string:"+seq4.getString());
		System.out.println("length:"+seq4.getLength());
		
	}

    public static void testFeatures() {
        String seqStr = "CAGGACGTAAGCCTAACGTA";
        Sequence seq3 =new Sequence(seqStr, 0, seqStr.length());
        byte order = 6;

        Map<String,Integer> result = Sequence.generateFeaturesMatrix(order);

        long start = System.currentTimeMillis();
        for(int i=0;i<1000000; i++) {
            result = Sequence.getFeatures(seq3, order, result);
        }
        long end = System.currentTimeMillis();

        for(String key : result.keySet()) {
            System.out.println(key + ":" +  result.get(key));
        }
        System.out.println("Time taken to generate 1000,000 runs of " + order + "th order features for "
                + seqStr.length() + "-char string: " + (end-start)/1000 + " seconds");
    }

    public static void testGenerateFeaturesMatrix() {
        long start = System.currentTimeMillis();
        Map<String,Integer> features = Sequence.generateFeaturesMatrix((byte)10);
        long end = System.currentTimeMillis();
        /*
        for(String key : features.keySet()) {
            //System.out.println(key + ":" +  features.get(key));
        }
        */

        System.out.println("Time taken to generate 10th order feature matrix is " + (end-start)/1000 + " seconds");
    }

}
