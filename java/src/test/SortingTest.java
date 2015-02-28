package test;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;

import base.CountObject;
import base.Sequence;

public class SortingTest
{

	/**
	 * @param args
	 */
	public static void main(String[] args)
	{
		
		int SIZE = 5000000;
		CountObject[] array = new CountObject[SIZE];
		Random r=new Random();
		int len = 16;
		for(int i=0;i<SIZE;i++)
		{
			array[i]=new CountObject(new Sequence(r.nextLong(),len),1);
		}

		System.out.println("Sorting...");
		long t1= System.currentTimeMillis();

		Arrays.sort(array, 0, SIZE ,new Comparator<CountObject>()
		{
			@Override
			public int compare(CountObject arg0, CountObject arg1)
			{
//				if (arg0 == null)
//					return 1;
//				if (arg1 == null)
//					return -1;
				return arg0.compareToByKey(arg1);
			}

		});
		long t2= System.currentTimeMillis();
		System.out.println(((t2-t1))+" miliseconds");
	}

}
