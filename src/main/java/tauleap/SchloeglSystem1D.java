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
 * with the number of S2 and S3 molecules are constant. Thus, the state space is, in fact, one-dimensional.
 * @author florian
 *
 */
public class SchloeglSystem1D extends ChemicalReactionNetwork implements MultiDim{

	/**
	 * number of S2-molecules (const.)
	 */
	double S2; 
	/**
	 * number of S3-molecules (const.)
	 */
	double S3;

	public SchloeglSystem1D(double[] c, double[] X0, double tau, double T) {
		this.c = c;
		this.X0 = new double [] {X0[0]};
		this.S2 = X0[1];
		this.S3 = X0[2];
		this.tau = tau;
		this.T = T;
		S = new double[][] { { 1, -1, 1, -1 } };
		init();
	}

	
	@Override
	public int compareTo(MarkovChainComparable m, int i) {
		if (!(m instanceof SchloeglSystem1D)) {
			throw new IllegalArgumentException("Can't compare an SchloeglSystem1D with other types of Markov chains.");
		}
		double mx;

		mx = ((SchloeglSystem1D) m).X[i];
		return (X[i] > mx ? 1 : (X[i] < mx ? -1 : 0));
	}

	public String toString() {
		StringBuffer sb = new StringBuffer("----------------------------------------------\n");
		sb.append(" SchloeglSystem1D:\n");
		sb.append("X0 =\t" + "{" + X0[0] + ", " + S2 + ", " + S3 + "}\n");
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
//		if(X[0] > 300)
//			return 1.0;
//		else
//			return 0.0;
	}

	

	@Override
	public double[] getPoint() {
		return X;
	}

	@Override
	public void computePropensities() {
		a[0] = 0.5 * c[0] * X[0] * (X[0] - 1.0) * S2;
		a[1] = c[1] * X[0] * (X[0] - 1.0) * (X[0] - 2.0) / 6.0;
		a[2] = c[2] * S3;
		a[3] = c[3] * X[0];
	}

	@Override
	public double getCoordinate(int j) {

return X[j];
		}



}
