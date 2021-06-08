package tauleap;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Scanner;

import umontreal.ssj.markovchainrqmc.MarkovChainComparable;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.Chrono;

import umontreal.ssj.util.sort.MultiDim01;
/**
 * Implements the Schloegl system, as a  #ChemicalReactionNetwork see, e.g., Beentjes and Baker '18, 
 * with the modification that the number of S2 and S3 molecules is not constant. Since the total number 
 * of molecules does not change through any reaction, the state space is, in fact, two-dimensional.
 * @author florian
 *
 */
public class SchloeglSystem2D extends ChemicalReactionNetwork implements MultiDim{

	/**
	 * total number of molecules
	 */
	double N0; // Total number of molecules

	

	public SchloeglSystem2D(double[] c, double[] X0, double tau, double T, double N0) {
		this.c = c;
		this.X0 = X0;
		this.tau = tau;
		this.T = T;
		S = new double[][] { { 1, -1, 1, -1 }, { -1, 1, 0, 0 }, { 0, 0, -1, 1 } };
		init();
		this.N0 = N0;
	}

	
	@Override
	public int compareTo(MarkovChainComparable m, int i) {
		if (!(m instanceof SchloeglSystem2D)) {
			throw new IllegalArgumentException("Can't compare an SchloeglSystem2D with other types of Markov chains.");
		}
		double mx;

		mx = ((SchloeglSystem2D) m).X[i];
		return (X[i] > mx ? 1 : (X[i] < mx ? -1 : 0));
	}

	public String toString() {
		StringBuffer sb = new StringBuffer("----------------------------------------------\n");
		sb.append(" SchloeglSystem2D:\n");
		sb.append("X0 =\t" + "{" + X0[0] + ", " + X0[1] + ", " + (N0 - X0[0] - X0[1]) + "}\n");
		sb.append("c =\t" + "{" + c[0] + ", " + c[1] + ", " + c[2] + "}\n");
		sb.append("T =\t" + T + "\n");
		sb.append("tau =\t" + tau + "\n");
		sb.append("steps =\t" + numSteps + "\n");
		sb.append("----------------------------------------------\n\n");

		return sb.toString();
	}

	@Override
	public double getPerformance() {
		return X[0];
//			return X[1];
	}

	

	@Override
	public double[] getPoint() {
		return X;
	}

	@Override
	public void computePropensities() {
		double x2 = (N0 - X[0] - X[1]);
		a[0] = 0.5 * c[0] * X[0] * (X[0] - 1.0) * X[1];
		a[1] = c[1] * X[0] * (X[0] - 1.0) * (X[0] - 2.0) / 6.0;
		a[2] = c[2] * x2;
		a[3] = c[3] * X[0];
	}

	@Override
	public double getCoordinate(int j) {

return X[j];
		}

	


}
