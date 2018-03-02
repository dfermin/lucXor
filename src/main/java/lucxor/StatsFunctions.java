/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import java.util.ArrayList;
import java.util.concurrent.ExecutionException;

import gnu.trove.list.array.TIntArrayList;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;

/**
 *
 * @author dfermin
 */
public class StatsFunctions {
	
	/***************
	 * Function returns the number of combinations of (n choose k) 
	 */
	public double combinatorial(double n, double k) {
		double ret = 1.0;
		
		// n = number of characters in alphabet
		// k = desired length of the words to be construted
		
		if(n <= 1.0) ret = 1.0;
		else {
			double diff = n - k;
			double facN = factorial(n);
			double facK = factorial(k);
			double facDiff = factorial(diff);
			
			ret = facN / ( facK * facDiff );
		}
		
		return ret;
	}

	
	// Recursive function to compute the Factorial (X!) of a given number
	public double factorial(double x) {
		double ret = 1.0;
		
		if( x <= 1 ) ret = 1;
		else {
			for(double i = 1; i <= x; i++) ret *= i;
		}
		return ret;
	}
	
	
	
	ArrayList<TIntArrayList> getAllCombinations(TIntArrayList candModSites, int k) {
		ArrayList<TIntArrayList> ret = new ArrayList();
		int L = candModSites.size(); // total number of values
		
		// http://code.google.com/p/combinatoricslib/#3._Simple_combinations

        ArrayList<Integer> vec = new ArrayList();
        for(int i : candModSites.toArray()) vec.add(i);

		ICombinatoricsVector<Integer> initialVector = Factory.createVector( vec );

		Generator<Integer> gen = Factory.createSimpleCombinationGenerator(initialVector, k);
		
		for(ICombinatoricsVector<Integer> combination : gen) {
			TIntArrayList curList = new TIntArrayList();
			for(int i : combination) curList.add(i);
			
			ret.add(curList);
			curList = null;
		}
		
		return ret;
	}

	
	
	double log_gaussianProb(double mu, double sigma2, double x) {
		double ret = 0d;
		ret = -0.5 * Math.pow((x-mu), 2.0) / sigma2 - 0.5 * Math.log( (2.0 * Math.PI * sigma2) );
		return ret;
	}

	
	// Function to compute the False Localization Rate of the PSMs
	void calcFLR() throws InterruptedException, ExecutionException {
		double maxDeltaScore = -1.0;
		FLRClass flr = new FLRClass();
		ArrayList<FLRClass> flrAry = new ArrayList();
		
		System.err.println("\nComputing False Localization Rate (FLR)");
		
		// Identify maxDeltaScore
		for(PSM psm : Globals.PSM_list) {
			if(psm.deltaScore > maxDeltaScore) maxDeltaScore = psm.deltaScore;

			if(psm.deltaScore > Constants.MIN_DELTA_SCORE) {
				if(psm.isDecoy) flr.decoyPSMs.add(psm);
				else flr.realPSMs.add(psm);
			}
		}
		
		flr.maxDeltaScore = maxDeltaScore;
		flr.prepArrays();
		
		flr.initializeTickMarks();
		flr.evalTickMarks(Constants.REAL);
		flr.evalTickMarks(Constants.DECOY);

		flr.calcBothFDRs();
		flr.setMinorMaps();
        flr.performMinorization();
		flr.assignFDRs();
//		flr.debugFLR();
	}
	
	
	
	
	
	
	
}
