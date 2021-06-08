package tauleap;

/**
 * Implements a #MultDimToOneDimSort for the #PKA class, where the score
 * function is given by the one-step look ahead importance function (OSLA),
 * explained in the reference paper. 
 * 
 * 
 * @author florian
 */

public class PKAOSLASort extends MultDimToOneDimSort  {
	int coord; //observed coordinate
	double tau; //stepsize
	double[] c; //reaction rates
	double[][]S; //stoichiometric matrix
	double[]a; //pronpensity functions
	int K; //number of reactions
	
	public PKAOSLASort(int coord, double[] c, double tau) {
		this.coord = coord;
		this.dimension=6;
		this.c = c;
		this.tau = tau;
		S = new double[][] { { -1, 1, 0, 0, 0, 0 }, { -2, 2, -2, 2, 0, 0 }, { 1, -1, -1, 1, 0, 0 },
			{ 0, 0, 1, -1, -1, 1 }, { 0, 0, 0, 0, 1, -1 }, { 0, 0, 0, 0, 2, -2 } };
		K=c.length;
//		K=2;
		a=new double[K];
	}
	
	public void computePropensities(double [] v) {
		a[0] = 0.5 * c[0] * v[0] * v[1] * (v[1] - 1.0);
		a[1] = c[1] * v[2];
		a[2] = 0.5 * c[2] * v[2] * v[1] * (v[1] - 1.0);
		a[3] = c[3] * v[3];
		a[4] = c[4] * v[3];
		a[5] = 0.5 * c[5] * v[4] * v[5] * (v[5] - 1.0);
	}
	
		@Override
		public double scoreFunction(double[] v) {
			computePropensities(v);

				double score = v[coord];
				for (int k = 0; k < K; k++) {
					score += S[coord][k] * tau * a[k];
				}

//			}
			return score;
		}

		@Override
		public String toString() {
			return "PKA-OSLA";
		}

}
