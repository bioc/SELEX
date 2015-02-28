package test;

import java.util.Arrays;
import java.util.LinkedList;

import base.CountObject;
import base.Sequence;
import base.Util;

public class RadixSortTest
{

	/**
	 * @param args
	 */
	public static void main(String[] args)
	{
		//System.out.println( Math.ceil(5.0/5.0));
		//System.out.println( Math.ceil(30.0/5.0));
		test();
	}

	public static void test()
	{
		CountObject[] list = new CountObject[]{
				new CountObject(new Sequence("ATTTT"),1),
				new CountObject(new Sequence("CTTAA"),1),
				new CountObject(new Sequence("GGGAT"),1),
				new CountObject(new Sequence("AGAAT"),1),
				new CountObject(new Sequence("ATTTG"),1),
				new CountObject(new Sequence("AATGG"),1),
		};
		
		Util.radixSort(list, list.length);
		
		System.out.println(Arrays.toString(list));
	}
	
	
}
