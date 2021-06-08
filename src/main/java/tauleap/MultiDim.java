
package tauleap;

/**
 * This interface represents a point or array of @f$d@f$ dimensions
 * in @f$\mathbb{R}^d@f$. The value of the @f$j@f$th dimension can be accessed
 * with the method {@link #getCoordinate() getCoordinate(j)}. 
 * 
 * 
 */
public interface MultiDim {

	/**
	 * This method returns the number dimensions of this point.
	 */
	public int dimension();

	/**
	 * Returns the @f$d@f$ coordinates of this point.
	 */
	public double[] getPoint();

	/**
	 * Returns the value of @f$j@f$th coordinate (or dimension).
	 */
	public double getCoordinate(int j);



}
