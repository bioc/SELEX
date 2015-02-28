package test;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import base.DebugLog;

public class ExecutorTest
{

	/**
	 * @param args
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws InterruptedException
	{

		ExecutorService service = Executors.newFixedThreadPool(2);
		// now submit our jobs
		for(int i=0;i<5;i++)
		{
			final int localI =i;
			service.submit(new Runnable() {
			    public void run() {
			        System.out.println("Thread started:"+localI);
			        try
					{
						Thread.sleep(3000);
					} catch (InterruptedException e)
					{
						// TODO Auto-generated catch block
						DebugLog.log(e);
					}
			        System.out.println("Thread ended:"+localI);
			    }
			});
		}
		
		service.shutdown();
		// now wait for the jobs to finish
		service.awaitTermination(10, TimeUnit.SECONDS);
	}

}

