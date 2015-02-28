package test;

import java.util.Arrays;
import java.util.Random;

public class Sorting
{
	
	public static void main(String[] args)
	{
		long t1=System.currentTimeMillis();
		int size = 30*1000*1000;
		Random r=new Random();
		long largeNums[] = new long[size];
		for(int i=0;i<size;i++)
		{
			largeNums[i]=r.nextLong();
		}
		
		Arrays.sort(largeNums);
		long t2=System.currentTimeMillis();
		System.out.println((t2-t1)/1000+" seconds");
	}

}
