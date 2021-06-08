package tauleap;

/**
 * Implements a #MultDimToOneDimSort for the #SchloeglSystem class with
 * the one-step look ahead importance function (OSLAIF), as explained
 * in the reference paper.
 * @author florian
 *
 */

public class SchloeglOSLASort extends MultDimToOneDimSort {

	int coord; //observed coordinate
	double tau; //stepsize
	double[] c; //reaction rates
	double N0; //total number of molecules in the system
	double[][]S; //stoichiometric matrix
	double[]a; //propensity functions
	int K; //numer of reactions
	
	
	public SchloeglOSLASort(int coord, double[] c, double tau, double N0) {
		this.coord = coord;
		this.dimension=2;
		this.c = c;
		this.tau = tau;
		S = new double[][] { { 1, -1, 1, -1 }, { -1, 1, 0, 0 }, { 0, 0, -1, 1 } };
		this.N0 = N0;
		K=c.length;
		a=new double[K];
	}
	
	public void computePropensities(double [] v) {
		double x2 = (N0 - v[0] - v[1]);
		a[0] = 0.5 * c[0] * v[0] * (v[0] - 1.0) * v[1];
		a[1] = c[1] * v[0] * (v[0] - 1.0) * (v[0] - 2.0) / 6.0;
		a[2] = c[2] * x2;
		a[3] = c[3] * v[0];
	}

	
	@Override
	public double scoreFunction(double[] v) {
		computePropensities(v);

//		for (int n = 0; n < N; n++) {
			double score = v[coord];
			for (int k = 0; k < K; k++) {
				score += S[coord][k] * tau * a[k];
			}

//		}
		return score;
	}

	@Override
	public String toString() {
		return "Schloegl-OSLA";
	}

}
