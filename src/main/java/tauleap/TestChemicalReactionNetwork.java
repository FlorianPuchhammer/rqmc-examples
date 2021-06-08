package tauleap;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import umontreal.ssj.functionfit.LeastSquares;
import umontreal.ssj.hups.BakerTransformedPointSet;
import umontreal.ssj.hups.CachedPointSet;
import umontreal.ssj.hups.IndependentPointsCached;
import umontreal.ssj.hups.KorobovLattice;
import umontreal.ssj.hups.LMScrambleShift;
import umontreal.ssj.hups.NestedUniformScrambling;
import umontreal.ssj.hups.PointSet;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RQMCPointSet;
import umontreal.ssj.hups.RandomShift;
import umontreal.ssj.hups.Rank1Lattice;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.hups.SortedAndCutPointSet;
import umontreal.ssj.hups.StratifiedUnitCube;
import umontreal.ssj.markovchainrqmc.ArrayOfComparableChains;
import umontreal.ssj.markovchainrqmc.MarkovChainComparable;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.Num;
import umontreal.ssj.util.sort.BatchSort;
import umontreal.ssj.util.sort.HilbertCurveBatchSort;
import umontreal.ssj.util.sort.HilbertCurveSort;
import umontreal.ssj.util.sort.MultiDimSort;
import umontreal.ssj.util.sort.SplitSort;

/**
 * Class that provides a main() method to simulate a #ChemicalReactionNetwork
 * with array-RQMC when one uses one and the same sorting algorithm in every
 * step.
 * 
 * @author florian
 *
 */

public class TestChemicalReactionNetwork {

	public static void main(String[] args) throws IOException {

		ChemicalReactionNetwork model;

		/*
		 * ******************* REVERSIBLE ISO
		 **********************/
// 		double epsInv = 1E2;
// 		double alpha = 1E-4;
// 		double[]c = {1.0,alpha};
// 		double[] x0 = {epsInv};
// 		double N0 = epsInv + epsInv/alpha;
// 		double T = 1.6; //T=25.6; T =102.4; T=819.2;
// 		double tau = T/8.0;  //tau = T/128.0;  tau = T/1024.0;
// 		
// 		 model = new ReversibleIsomerization(c,x0,tau,T,N0);
		
 


		/*
		 * ******************* SCHLOEGL System
		 **********************/
		double[] c = { 3E-7, 1E-4, 1E-3, 3.5 };
		double[] x0 = { 250.0, 1E5};
		double N0 = 2E5 + 1E5 + 250.0;
		double T = 4.0; // T=32.0;
		double tau =T/16.0; tau = T/128.0;

		model = new SchloeglSystem2D(c, x0, tau, T,N0);
//		model = new SchloeglSystem01(c, x0, tau, T,N0); //for Hilbert Curve-related sort. Requires cls-data files (or needs to be adapted).


		/*
		 * ******************* PKA
		 **********************/
// 		double cFac = 2.0 /(11.0 * 6.02214);
// 		double[] c = { 8.696E-5 * cFac, 0.02, 1.154E-4*cFac, 0.02, 0.016, 0.0017*cFac };
// 		double[] x0 = { 33000.0, 33030.0, 1100.0, 1100.0, 1100.0, 1100.0 };
// //		double[] x0 = { 1100.0, 33030.0, 1100.0, 1100.0, 1100.0, 33000.0 };
// 		double  T=0.05;
// 		double tau = T / 256.0;
// 
// 		model = new PKA(c, x0, tau, T);
// //		model = new PKA01(c, x0, tau, T); //for Hilbert Curve-related sort. Requires cls-data files (or needs to be adapted).


		model.init();
		System.out.println(model.toString());
		

		String modelDescription = "schloegl2D";


		ArrayOfComparableChains chain = new ArrayOfComparableChains(model); //class for array-RQMC experiments

		
		int [] N = {8192, 16384,32768, 65536, 131072, 262144, 524288}; //number of points
		int mink = 13; //base-2-logarithm of smallest "N".

		
		int numSets = N.length;
		int m = 100; //number of independent repetitions

		// generating vector for lattice using order dependent weights 0.6^k
		int[][] aa = {
				{ 1, 3455, 1967, 1029, 2117, 3871, 533, 2411, 1277, 2435, 1723, 3803, 1469, 569, 1035, 3977, 721, 797,
						297, 1659 }, // N=2^13

				{ 1, 6915, 3959, 7743, 3087, 5281, 6757, 3369, 7107, 6405, 7753, 1641, 3613, 1819, 5827, 2087, 4417,
						6909, 5623, 4739 }, //  N=2^14

				{ 1, 12031, 14297, 677, 6719, 15787, 10149, 7665, 1017, 2251, 12105, 2149, 16273, 14137, 8179, 6461,
						15051, 6593, 12763, 8497 }, //  N=2^15

				{ 1, 19463, 8279, 14631, 12629, 26571, 30383, 1337, 6431, 3901, 12399, 20871, 5175, 3111, 26857, 15111,
						22307, 30815, 25901, 27415 }, //  N=2^16

				{ 1, 38401, 59817, 33763, 32385, 2887, 45473, 48221, 3193, 63355, 40783, 37741, 54515, 11741, 10889,
						17759, 6115, 18687, 19665, 26557 }, //  N=2^17

				{ 1, 100135, 28235, 46895, 82781, 36145, 36833, 130557, 73161, 2259, 3769, 2379, 80685, 127279, 45979,
						66891, 8969, 56169, 92713, 67743 }, //  N=2^18

				{ 1, 154805, 242105, 171449, 27859, 76855, 183825, 38785, 178577, 18925, 260553, 130473, 258343, 79593,
						96263, 36291, 2035, 198019, 15473, 148703 }, //  N=2^19

				{ 1, 387275, 314993, 50301, 174023, 354905, 303021, 486111, 286797, 463237, 211171, 216757, 29831,
						155061, 315509, 193933, 129563, 276501, 395079, 139111 } //  N=2^20
		};

////		 generating vector for lattice using order dependent weights 0.05^k
//				int[][] aa = {
//						{ 1, 3455, 1899, 2921, 3663, 2823, 3977, 2761, 255, 845, 3029, 3831, 2089, 3691, 1771, 3907, 337, 3735,
//								1373, 1795 }, //  N=2^13
//
//						{ 1, 6915, 4877, 7479, 1203, 3941, 2159, 3225, 5219, 6307, 2643, 633, 7139, 869, 7239, 7019, 8151, 3853,
//								8019, 5731 }, //  N=2^14
//
//						{ 1, 12033, 3801, 5023, 10647, 14127, 12751, 7461, 11901, 1167, 14349, 1951, 2209, 7397, 2505, 5675,
//								12195, 1801, 7707, 13443 }, //  N=2^15
//
//						{ 1, 25015, 11675, 7425, 3289, 17821, 5649, 32161, 10285, 12031, 26337, 13403, 14547, 18661, 7993, 1299,
//								15111, 12735, 13129, 12655 }, //  N=2^16
//
//						{ 1, 38401, 48799, 17301, 59639, 20297, 26805, 53109, 4365, 14055, 5023, 48499, 37937, 5155, 44255,
//								61671, 11409, 38529, 61887, 19183 }, //  N=2^17
//
//						{ 1, 96407, 36479, 31333, 63411, 80945, 24597, 41083, 70179, 42983, 62013, 48035, 80011, 105415, 108151,
//								68869, 104973, 20719, 72257, 59193 }, //  N=2^18
//
//						{ 1, 154805, 243089, 211205, 258913, 18107, 174117, 67287, 3585, 155767, 31401, 154275, 35513, 36509,
//								162377, 51021, 88413, 190981, 145989, 257551 }, //  N=2^19
//
//						{ 1, 387275, 457903, 282967, 117983, 355873, 439959, 109733, 382437, 297385, 267803, 68841, 343399,
//								171303, 420841, 136437, 423733, 355591, 415917, 406205 } //  N=2^20
//				};
		
		
		int numSteps = (int) (T / tau); //number of steps

		ArrayList<MultiDimSort> sortList = new ArrayList<MultiDimSort>(); //which sorts to use. The experiment is carried out for each sort from this list.
		ArrayList<Integer> sortCoordPtsList = new ArrayList<Integer>(); //how many coordinates of the points are used for sorting (one entry per entry in sortList)
		ArrayList<MultiDimSort> sortPts = new ArrayList<MultiDimSort>(); //how are the points sorted. For most standard sorts the same as in "sortList". For Hilbert Curve sort and importance function based sorts usually different. One sort per entry in sortList



    // Some choices of possible importance-function-based sorts
    // Multiple choices are possible. In that case, the simulation will be run for every
    // selected sorting method.
    
    
//		sortList.add(new PKAbyPKAcSort());
//		sortList.add(new PKAOSLASort(4,c,tau));
		sortList.add(new SchloeglOSLASort(0,c,tau,N0));
//	 	sortList.add(new BatchSort(new double[] { 1.0 })); //only first coordinate
// 

	//When using importance-function-based sorts, the point sets will be sorted by their 1st coordinate. 
	//For each sort selected above, add a one dimensional sort to "sortPts" and a "1! to sortCoordPtsList.
	
		sortPts.add(new BatchSort(new double[] { 1.0 })); //I think batch is faster than split
		sortCoordPtsList.add(1);

		
		
		//standard sorts
		
		/* BATCH SORT */

		double[] batchExp = {0.5,0.5};
	
		sortList.add(new BatchSort<MarkovChainComparable>(batchExp));
		sortPts.add(new BatchSort<MarkovChainComparable>(batchExp));		
		sortCoordPtsList.add(model.dimension());
		modelDescription = "schloegl-batch-sort";

		/* SPLIT SORT */
		sortList.add(new SplitSort<MarkovChainComparable>(model.dimension()));
		sortPts.add(new SplitSort<MarkovChainComparable>(model.dimension()));
		sortCoordPtsList.add(model.dimension());
		modelDescription = "schloegl-split-sort";


		//variables for writing results
		StringBuffer sb = new StringBuffer("");
		String str;
		String outFile = modelDescription + ".txt";

		RandomStream stream = new MRG32k3a(); //stream for random numbers
		
		//Variable declaration for RQMC point set construction
		RQMCPointSet[] rqmcPts; //container for RQMC points
		PointSet[] pointSets = new PointSet[numSets]; //container for Point Sets (used to construct rqmcPts)
		PointSetRandomization rand; //point set randomization (used to construct rqmcPts)
		RQMCPointSet prqmc; //temporary variable (used to construct rqmcPts)
		int i, s;

		int nMC = (int) 1E6; // number of samples to estimate MC variance.
		Tally statMC = new Tally(); // statistical container for MC run
		statMC.init();
		str = model.simulRunsFormat(nMC, model.numSteps, stream, statMC); //run with MC
		double varMC = statMC.variance() ;
		sb.append(str);
		System.out.println(str);

		i = 0; // Sorts indexed by i
		for (MultiDimSort sort : sortList) { //for each sort in sortList
			str = "****************************************************\n";
			str += "*\t" + sort.toString() + "\n";
			str += "****************************************************\n\n";
			sb.append(str);
			System.out.println(str);
			ArrayList<RQMCPointSet[]> listP = new ArrayList<RQMCPointSet[]>();
			
			//uncomment the corresponding blocks for the desired point sets. Multiple selections
			//are possible.
			
			// Independent points (Monte Carlo) implemented as RQMC points. Only for debugging purposes.
//			 rqmcPts = new RQMCPointSet[numSets];
//			 for (s = 0; s < numSets; ++s) {
//			 pointSets[s] = new IndependentPointsCached(N[s], model.K + model.N);
//			 rand = new RandomShift(stream);
//			 prqmc = new RQMCPointSet(pointSets[s], rand);
//			 rqmcPts[s] = prqmc;
//			 }
//			 rqmcPts[0].setLabel("Independent points");
//			 listP.add(rqmcPts);
//
//			// Stratification
//				rqmcPts = new RQMCPointSet[numSets];
//				int k;
//				for (s = 0; s < numSets; ++s) {
//					k = (int) Math.round(Math.pow(Num.TWOEXP[s + mink], 1.0 / (double) (sortCoordPtsList.get(i) + model.K)));
//					pointSets[s] = new StratifiedUnitCube(k, sortCoordPtsList.get(i) + model.K);
//					// Here the points must be sorted at each step, always.
//					// In the case of Hilbert map, the points should be 2d and sorted
//					// based on one coordinate,
//					// whereas the states are 2d and sorted by the Hilbert sort.
//					rand = new RandomShift(stream);
//					prqmc = new RQMCPointSet(pointSets[s], rand);
//					rqmcPts[s] = prqmc;
//				}
//				rqmcPts[0].setLabel("Stratification");
//				listP.add(rqmcPts);

// 			 Lattice + Shift 
			rqmcPts = new RQMCPointSet[numSets];
			for (s = 0; s < numSets; ++s) {

				//for array-RQMC define point set as SortedAndCutPointSet
				pointSets[s] = new SortedAndCutPointSet(
						new Rank1Lattice(N[s], aa[s], sortCoordPtsList.get(i) + model.K), sortPts.get(i)); 

				rand = new RandomShift(stream);
				prqmc = new RQMCPointSet(pointSets[s], rand);
				rqmcPts[s] = prqmc;
			}
			rqmcPts[0].setLabel("lattice+shift");
			listP.add(rqmcPts);

			// Rank1Lattice + baker + shift (or shift + baker?)
			rqmcPts = new RQMCPointSet[numSets];
			for (s = 0; s < numSets; ++s) {

				// The points are sorted here, but only once.
				pointSets[s] = new SortedAndCutPointSet(
						new BakerTransformedPointSet(new Rank1Lattice(N[s], aa[s], sortCoordPtsList.get(i) + model.K)),
						sortPts.get(i));
//							 
				rand = new RandomShift(stream);
				prqmc = new RQMCPointSet(pointSets[s], rand);
				rqmcPts[s] = prqmc;
			}
			rqmcPts[0].setLabel("lattice + baker ");
			listP.add(rqmcPts);

//			 Sobol + LMS
			rqmcPts = new RQMCPointSet[numSets];
			for (s = 0; s < numSets; ++s) {

				pointSets[s] = new SortedAndCutPointSet(
						new SobolSequence(s + mink, 31, sortCoordPtsList.get(i) + model.K), sortPts.get(i));

				rand = new LMScrambleShift(stream);
				prqmc = new RQMCPointSet(pointSets[s], rand);
				rqmcPts[s] = prqmc;
			}
			rqmcPts[0].setLabel("Sobol+LMS");
			listP.add(rqmcPts);

			// Sobol+NUS
			rqmcPts = new RQMCPointSet[numSets];
			for (s = 0; s < numSets; ++s) {

				CachedPointSet p = new CachedPointSet(
                    new SobolSequence(s + mink, 31, sortCoordPtsList.get(i) + model.K));
				p.setRandomizeParent(false);
				// The points are sorted here, but only once.
				pointSets[s] = new SortedAndCutPointSet(p, sortPts.get(i));

				rand = new NestedUniformScrambling(stream);
				prqmc = new RQMCPointSet(pointSets[s], rand);
				rqmcPts[s] = prqmc;
			}
			rqmcPts[0].setLabel("Sobol+NUS");
			listP.add(rqmcPts);

			//Run the experiment for every point set.
			for (RQMCPointSet[] ptSeries : listP) {
				String label = ptSeries[0].getLabel();
				str = label;
				str += "\n-----------------------------\n";
				sb.append(str + "\n");
				str += "*****************************\n";
				str += model.toString() + "*****************************\n\n";
				System.out.println(str);
				// sortedCoords gives the number of coordinates by which the point set is sorted in EACH step.
				// most constructions have sortedCoords=0, because they are sorted once in the beginning and the 
				// initial sortedCoordPtsList.get(i) are not randomized again.
				// If we use Stratification, then we need to sort point set in every step
				int sortedCoords = label.startsWith("St") ? sortCoordPtsList.get(i) : 0; 
				str = (chain.testVarianceRateFormat(ptSeries, sort, sortedCoords, model.numSteps, m, varMC,
						modelDescription + "-" + label, label));
				System.out.println(str);
				sb.append(str + "\n");

			}
			i++;

		}
		FileWriter file = new FileWriter(outFile);
		file.write(sb.toString());
		file.close();

	}

}
