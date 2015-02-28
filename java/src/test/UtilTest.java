package test;

import config.Sample;
import base.Util;

public class UtilTest
{
	/**
	 * @param args
	 */
	public static void main(String[] args)
	{
		cloneSample();
		
//		System.out.println(Util.getMD5(str));
	}
	
	public static void cloneSample()
	{
		Sample s =new Sample();
		s.setLeftBarcode("CCC");
		s.setRightBarcode("ZZZZ");
		s.setVariableRegionLength(10);
		System.out.println(s);
		
		Sample s2 = Util.cloneSample(s);
		System.out.println(s2);
	}
}
