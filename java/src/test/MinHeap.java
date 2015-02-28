package test;

import java.util.Arrays;

import base.DebugLog;

/**
 * This class is not used in SELEX. It serves as a reference for implementing ArraysMergerHeap
 *
 */
public class MinHeap
{
	public MinHeap(int[] array)
	{
		this.array = array;
		this.len= array.length;
	}
	
	private int[] array;
	private int len;
	
	public void minHeapify(int i)
	{
		int l=2*i;
		int r=l+1;
		int smallest = i;
		if(l<len && array[l] < array[i])
		{
			smallest=l;
		}
		if(r<len && array[r] < array[smallest])
		{
			smallest=r;
		}
		
		if(smallest!=i)
		{
			exchange(smallest,i);
			minHeapify(smallest);
		}
	}
	
	public void buildMinHeap()
	{
		for(int i=(len-1)/2;i>=0;i--)
		{
			//DebugLog.log(Arrays.toString(array));
			minHeapify(i);
		}
	}
	
	public void sort()
	{
		//DebugLog.log(Arrays.toString(array));
		buildMinHeap();
		DebugLog.log(Arrays.toString(array));
		for(int i=len-1;i>=0;i--)
		{
			DebugLog.log(array[0]);
			exchange(0,i);
			len--;
			minHeapify(0);
		}
	}
	
	private void exchange(int i,int j)
	{
		int temp=array[i];
		array[i]=array[j];
		array[j] =temp;
	}
}
