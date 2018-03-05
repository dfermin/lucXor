/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import gnu.trove.list.array.TIntArrayList;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;
import org.paukov.combinatorics.util.Util;

/**
 * @author dfermin
 */
class StatsFunctions {

  /**
   * The number of combinations of (n choose k).
   * @param n = number of characters in alphabet
   * @param k = desired length of the words to be construted
   */
  public static long combinatorial(int n, int k) {
    return CombinatoricsUtils.binomialCoefficient(n, k);
  }

  public static ArrayList<TIntArrayList> getAllCombinations(TIntArrayList candModSites, int k) {
    ArrayList<TIntArrayList> ret = new ArrayList<>();

    // http://code.google.com/p/combinatoricslib/#3._Simple_combinations

    ArrayList<Integer> vec = new ArrayList<>();
    for (int i : candModSites.toArray()) {
      vec.add(i);
    }

    ICombinatoricsVector<Integer> initialVector = Factory.createVector(vec);

    Generator<Integer> gen = Factory.createSimpleCombinationGenerator(initialVector, k);

    for (ICombinatoricsVector<Integer> combination : gen) {
      TIntArrayList curList = new TIntArrayList();
      for (int i : combination) {
        curList.add(i);
      }

      ret.add(curList);
    }

    return ret;
  }


  public static double log_gaussianProb(double mu, double sigma2, double x) {
    return -0.5 * Math.pow((x - mu), 2.0) / sigma2 - 0.5 * Math.log((2.0 * Math.PI * sigma2));
  }


  // Function to compute the False Localization Rate of the PSMs
  public static void calcFLR() throws InterruptedException, ExecutionException {
    double maxDeltaScore = -1.0;
    FLRClass flr = new FLRClass();

    System.err.println("\nComputing False Localization Rate (FLR)");

    // Identify maxDeltaScore
    for (PSM psm : Globals.PSM_list) {
      if (psm.deltaScore > maxDeltaScore) {
        maxDeltaScore = psm.deltaScore;
      }

      if (psm.deltaScore > Constants.MIN_DELTA_SCORE) {
        if (psm.isDecoy) {
          flr.decoyPSMs.add(psm);
        } else {
          flr.realPSMs.add(psm);
        }
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
