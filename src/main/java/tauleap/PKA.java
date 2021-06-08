package package tauleap;


import umontreal.ssj.markovchainrqmc.MarkovChainComparable;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.util.sort.MultiDim01;


/**
 * Implements the #ChemicalReactionNetwork for the cyclic AMP activation of PKA,
 * see Koh and Blackwell '12, Strehl and Ilie '15.
 * 
 * @author florian
 *
 */
public class PKA extends ChemicalReactionNetwork implements MultiDim {

	/**
	 * Constructor
	 * @param c reaction rates
	 * @param X0 initial states
	 * @param tau step length
	 * @param T final time
	 */
	public PKA(double[] c, double[] X0, double tau, double T) {
		this.c = c;
		this.X0 = X0;
		this.tau = tau;
		this.T = T;
		S = new double[][] { 
			{ -1, 1, 0, 0, 0, 0 }, 
			{ -2, 2, -2, 2, 0, 0 }, 
			{ 1, -1, -1, 1, 0, 0 },
			{ 0, 0, 1, -1, -1, 1 }, 
			{ 0, 0, 0, 0, 1, -1 }, 
			{ 0, 0, 0, 0, 2, -2 } };
		init();
	}



	
	public String toString() {
		StringBuffer sb = new StringBuffer("----------------------------------------------\n");
		sb.append(" cAMP activation of PKA:\n");
		sb.append("Number of reactions K = " + K + "\n");
		sb.append("Number of species N = " + N + "\n");
		sb.append("X0 =\t" + "{" + X0[0]);
		for (int i = 1; i < X0.length; i++)
			sb.append(", " + X0[i]);
		sb.append("}\n");

		sb.append("c =\t" + "{" + c[0]);
		for (int i = 1; i < c.length; i++)
			sb.append(", " + c[i]);
		sb.append("}\n");
		sb.append("T =\t" + T + "\n");
		sb.append("tau =\t" + tau + "\n");
		sb.append("steps =\t" + numSteps + "\n");
		sb.append("----------------------------------------------\n\n");

		return sb.toString();
	}

	@Override
	public int compareTo(MarkovChainComparable m, int i) {
		if (!(m instanceof PKA)) {
			throw new IllegalArgumentException("Can't compare a PKA with other types of Markov chains.");
		}
		double mx;

		mx = ((PKA) m).X[i];
		return (X[i] > mx ? 1 : (X[i] < mx ? -1 : 0));
	}

	@Override
	public double[] getPoint() {
		return X;
	}

	@Override
	public double getCoordinate(int j) {
		return X[j];

	}

	@Override
	public void computePropensities() {
		a[0] = 0.5 * c[0] * X[0] * X[1] * (X[1] - 1.0);
		a[1] = c[1] * X[2];
		a[2] = 0.5 * c[2] * X[2] * X[1] * (X[1] - 1.0);
		a[3] = c[3] * X[3];
		a[4] = c[4] * X[3];
		a[5] = 0.5 * c[5] * X[4] * X[5] * (X[5] - 1.0);
	}

	@Override
	public double getPerformance() {
		return X[0]; // PKA
//		return X[1]; //cAMP
//		return X[2]; //PKA-cAMP2
//		return X[3];
//		return X[4]; //PKAr
//		return X[5]; //PKAc
	}

}
