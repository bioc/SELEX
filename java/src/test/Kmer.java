package test;

/**
 * This class is not used in the project.
 * @author homliu
 *
 */
public class Kmer implements Comparable<Kmer>
{
	String base;
	int start;
	int length;
	
	public Kmer(String base, int start, int length)
	{
		super();
		this.base = base;
		this.start = start;
		this.length = length;
	}
	
	public String toString()
	{
		return this.base.substring(this.start, this.start+this.length);
	}

	private Integer hashCode = null;
	@Override
	public int hashCode()
	{
		if(hashCode ==null)
		{
			hashCode =this.toString().hashCode();
		}
		return hashCode;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		Kmer other = (Kmer) obj;
		if(this.length != other.length)
			return false;
		for(int i=0;i<this.length;i++)
		{
			if(base.charAt(i+this.start)!=other.base.charAt(i+other.start))
				return false;
		}
		return true;
	}

	@Override
	public int compareTo(Kmer o)
	{
		return this.toString().compareTo(o.toString());
	}

	
}
