package tauleap;

/**
 * Implements a #MultDimToOneDimSort for the #PKA class, where the score
 * function is the projection onto PKAc.
 * Using a batch sort with corresponding batch-sizes/batch-components
 * is probably faster!
 * 
 * @author florian
 *
 */
public class PKAbyPKAcSort extends MultDimToOneDimSort {



	public PKAbyPKAcSort() {
		this.dimension = 6;
	}

	

	@Override
	public double scoreFunction(double[] v) {
		return v[5];
	}

	@Override
	public String toString() {
		return "byPKAc-Sort";
	}

}
