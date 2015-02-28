package test;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class RegTest
{

	/**
	 * @param args
	 */
	public static void main(String[] args)
	{

	      Pattern r = Pattern.compile(
	    		"TGA[T|C].{5,6}TGA[T|C]|" +
	      		"TGA[T|C].{5,6}[A|G]TCA|" +
	      		"[A|G]TCA.{5,6}TGA[T|C]|" +
	      		"[A|G]TCA.{5,6}[A|G]TCA|" +
	      		"TAATTA.{1,4}TAATTA");

	      //String line = "TGATDDDDDTGAT";
	      //String line = "DGACDDDDDDTGAT";
	      //String line = "AAAAACGATCAA";
	      //String line ="TGGAAAAATCAATAAATCACCACGTCTCGTATGCCG";
	      String line ="TGGTTTTGATTAATGATCACCAGCTGTCGTATGCCG";
	      // Now create matcher object.
	      Matcher m = r.matcher(line);
	      if (m.find( )) {
	         System.out.println("Found value: " + m.group(0) );
	      } else {
	         System.out.println("NO MATCH");
	      }
	}

}
