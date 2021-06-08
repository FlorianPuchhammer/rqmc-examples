package tauleap;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Arrays;
import java.util.Scanner;

import umontreal.ssj.markovchainrqmc.MarkovChainComparable;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.sort.MultiDim01;

//Implements the SchloeglSystem2D and maps each state into [0,1) via a logistic transformation.
//The means and variances needed for this transformation are estimated from data files.
//NOTE: this class requires that data files are stored in "../../data/".
//The files are named "MCData_Step_s.csv", where "s" stands for the number of the respective step.

public class SchloeglSystem01 extends SchloeglSystem2D implements MultiDim01 {

	double[][] means;
	double[][] stdDevs;

	public SchloeglSystem01(double[] c, double[] X0, double tau, double T, double N0) throws FileNotFoundException {
		super(c, X0, tau, T, N0);
		readMeanVar("../../data/");

	}

	private void readMeanVar(String path) throws FileNotFoundException {
		int nSample = 524288;
		means = new double[numSteps][N];
		stdDevs = new double[numSteps][N];
		Arrays.fill(stdDevs[0], 1.0);
		Arrays.fill(means[0], 0.0);
		Tally[] tallyStep=new Tally[nSample];
		for (int j = 0; j < N; ++j)
			tallyStep[j] = new Tally("");
		Scanner sc;
		String[] line;
		for (int s = 1; s < numSteps; ++s) {
			for (int j = 0; j < N; ++j)
				tallyStep[j].init();
			sc = new Scanner(new BufferedReader(new FileReader(path + "MCData_Step_" + s + ".csv")));
			for (int i = 0; i < nSample; ++i) {
				// read line
				line = sc.nextLine().trim().split(",");
				// add value to Tally
				for (int j = 0; j < N; ++j)
					tallyStep[j].add(Double.parseDouble(line[j]));
			}
			// write tally-stats into arrays
			for (int j = 0; j < N; ++j) {
				means[s][j] = tallyStep[j].average();
				stdDevs[s][j] = tallyStep[j].standardDeviation();
			}

			sc.close();
		}
	}

	@Override
	public int compareTo(MarkovChainComparable m, int i) {
		if (!(m instanceof SchloeglSystem01)) {
			throw new IllegalArgumentException("Can't compare an SchloeglSystem01 with other types of Markov chains.");
		}
		double mx;

		mx = ((SchloeglSystem01) m).X[i];
		return (X[i] > mx ? 1 : (X[i] < mx ? -1 : 0));
	}

	public String toString() {
		StringBuffer sb = new StringBuffer("----------------------------------------------\n");
		sb.append(" SchloeglSystem01:\n");
		sb.append("X0 =\t" + "{" + X0[0] + ", " + X0[1] + ", " + (N0 - X0[0] - X0[1]) + "}\n");
		sb.append("c =\t" + "{" + c[0] + ", " + c[1] + ", " + c[2] + "}\n");
		sb.append("T =\t" + T + "\n");
		sb.append("tau =\t" + tau + "\n");
		sb.append("steps =\t" + numSteps + "\n");
		sb.append("----------------------------------------------\n\n");

		return sb.toString();
	}

	@Override
	public double[] getPoint() {
		double[] state01 = new double[N];
		for (int i = 0; i < N; i++)
			state01[i] = getCoordinate(i);
		return state01;
	}
	


	@Override
	public double getCoordinate(int j) {
		double xObar=means[step][j] + 2.0 * stdDevs[step][j];
		double xUbar = means[step][j] - 2.0*stdDevs[step][j];

		return 1.0/ (1.0 + Math.exp(- (X[j] - xUbar) / (xObar - xUbar) )  );
	}

}
