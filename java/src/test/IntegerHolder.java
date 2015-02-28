package test;

public class IntegerHolder 
{
	private Integer value;
	public IntegerHolder(int v)
	{
		this.value=v;
	}
	public Integer getInt()
	{
		return value;
	}
	public void setInt(Integer value)
	{
		this.value = value;
	}
	
	public void add(Integer v)
	{
		this.value += v;
	}
	
}
