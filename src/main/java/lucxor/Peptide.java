/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import java.util.ArrayList;
import java.util.Map.Entry;
import org.apache.commons.math3.util.FastMath;

/**
 * @author dfermin
 */
class Peptide {

  String peptide;
  String modPeptide;
  int charge;
  private TIntDoubleHashMap modPosMap; // total mass of each modified amino acid
  THashMap<Integer, String> nonTargetMods; // holds the mod coordinates for non-target mods
  private THashMap<String, Double> b_ions = null;
  private THashMap<String, Double> y_ions = null;

  int pepLen = 0;  // peptide's length
  int numRPS; // number of _R_eported _P_hospho_S_ites
  int numPPS; // number of _P_otential _P_hospho_S_ites

  private long numPermutations;
  private long numDecoyPermutations = 0;
  double score;

  ArrayList<PeakClass> matchedPeaks;

  Peptide() {
    modPosMap = new TIntDoubleHashMap(2);
    nonTargetMods = new THashMap<>(2);
    charge = 0;
    numPermutations = 0;
    score = 0;
  }


  /**
   * This initialization function is called by the PSM class when reading
   * a new entry from a pepXML file.
   */
  void initialize(TIntDoubleHashMap modCoordMap) {

    pepLen = peptide.length();

    peptide = peptide.toUpperCase();

    // remove any modifications from modCoordMap that are fixed modifications
    for (int pos : modCoordMap.keys()) {
      double mass = modCoordMap.get(pos);

      if ((pos == Constants.NTERM_MOD) || (pos == Constants.CTERM_MOD)) {
        modPosMap.put(pos, mass);
      } else {
        String C = Character.toString(peptide.charAt(pos));
        String c = C.toLowerCase();

//        if (Globals.fixedModMap.containsKey(C)) {
//          continue;
//        } else
        if (Globals.targetModMap.containsKey(C)) {
          modPosMap.put(pos, mass);
        } else if (Globals.varModMap.containsKey(c)) {
          nonTargetMods.put(pos, c); // record position of mod.
        }
      }
    }

    // construct the single-character representation of all the modifications
    // in this peptide sequence
    modPeptide = "";
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < pepLen; i++) {
      String AA = Character.toString(peptide.charAt(i));

      if (nonTargetMods.containsKey(i)) {
        AA = nonTargetMods.get(i);
      } else if (modPosMap.containsKey(i)) { // this residue is a target site
        String c = AA;
        if (!Globals.isDecoyResidue(c)) {
          AA = c.toLowerCase();
        }
      }
      sb.append(AA);
    }
    modPeptide = sb.toString();

    // Determine how many potential target-mod sites this peptide has
    numPPS = 0;
    for (int i = 0; i < pepLen; i++) {
      String AA = Character.toString(peptide.charAt(i));
      if (Globals.targetModMap.containsKey(AA)) {
        numPPS++;
      }
    }

    // Determine how many reported target-mod sites this peptide has
    numRPS = 0;
    for (int i = 0; i < pepLen; i++) {
      String aa = Character.toString(modPeptide.charAt(i));
      int score = 0;
      if (Globals.targetModMap.containsKey(aa.toUpperCase())) {
        score++;
      }
      if (Character.isLowerCase(aa.charAt(0))) {
        score++;
      }

      if (score == 2) {
        numRPS++;
      }
    }

    numPermutations = StatsFunctions.combinatorial(numPPS, numRPS);
    calcNumDecoyPermutations();
  }


  // This initialize function is really for declaring a new permutation of a peptide sequence
  // In this function call, the origPepSeq string is the original *MODIFIED* string
  void initialize(String origPep, String curSeq, int z, TIntDoubleHashMap mcp) {

    peptide = origPep;
    pepLen = peptide.length();
    charge = z;
    modPosMap = mcp;

    // deal with the N-terminus and C-terminus values
    if (curSeq.startsWith("[")) {
      modPosMap.put(Constants.NTERM_MOD, Globals.ntermMass);
      modPeptide = curSeq;
    } else {
      modPeptide = curSeq;
    }

    if (curSeq.endsWith("]")) {
      modPosMap.put(Constants.CTERM_MOD, Globals.ctermMass);
    }
  }


  // Function returns the modified peptide assigned to this object in "TPP format"
  public String getModPepTPP() {
    StringBuilder ret = new StringBuilder();
    int pepLen_local = pepLen;

    // Append N-term modification (if any is present)
    if (modPeptide.startsWith("[")) {
      int d = (int) Math.round(Globals.ntermMass);
      ret = new StringBuilder("n[" + String.valueOf(d) + "]");
      pepLen_local += 1; // need this to capture last residue in this function
    }

    // Now label peptide residues
    for (int i = 0; i < pepLen_local; i++) {
      String aa = Character.toString(modPeptide.charAt(i));

      if (aa.equals("[") || aa.equals("]")) {
        continue;
      }

      // this is a modified amino acid
      if (Character.isLowerCase(modPeptide.charAt(i))) {
        ret.append(Globals.getTPPresidue(aa));
      } else {
        ret.append(aa);
      }
    }

    // Append C-term modification (if any is present)
    if (modPeptide.endsWith("]")) {
      int d = (int) Math.round(Globals.ctermMass);
      ret.append("c[").append(String.valueOf(d)).append("]");
    }

    return ret.toString();
  }


  void build_ion_ladders() {
    int ionSeriesLen = charge * pepLen;
    b_ions = new THashMap<>(ionSeriesLen);
    y_ions = new THashMap<>(ionSeriesLen);

    String b, y;
    double bm, ym; // ion mass
    double bmz, ymz; // ion mass divided by charge z
    double ntermM = 0d, ctermM = 0d;

    final int minLen = 2; // minimum number of residues a fragment must contain

    if (modPosMap.containsKey(Constants.NTERM_MOD)) {
      ntermM = Globals.ntermMass;
    }
    if (modPosMap.containsKey(Constants.CTERM_MOD)) {
      ctermM = Globals.ctermMass;
    }

    for (int Z = 1; Z < charge; Z++) {
      for (int i = 1; i < pepLen; i++) { // starting at 1 ensures our shortest fragment ion is 2 char long

        String prefix = modPeptide.substring(0, i);
        String suffix = modPeptide.substring(i);

        String prefix_len = Integer.toString(prefix.length());
        String suffix_len = Integer.toString(suffix.length());

        if (prefix.length() >= minLen) {
          b = "b" + prefix_len + ":" + prefix;

          bm = Globals.getFragmentIonMass(b, Z, 0.0) + ntermM;
          bmz = bm / Z;

          if (Z > 1.0) {
            b += "/+" + Integer.toString(Z);
          }

          if (bmz > Globals.minMZ) {
            b_ions.put(b, bmz);

            // Determine if this ion sequence can under go a neutral loss
            // if it can, record the neutral loss variant
            record_NL_ions(b, Z, bm);
          }

        }

        if ((suffix.length() >= minLen) && !suffix.equalsIgnoreCase(modPeptide)) {
          y = "y" + suffix_len + ":" + suffix;
          ym = Globals.getFragmentIonMass(y, Z, Constants.WATER) + ctermM;
          ymz = ym / Z;

          if (Z > 1.0) {
            y += "/+" + Integer.toString((int) Z);
          }

          if (ymz > Globals.minMZ) {
            y_ions.put(y, ymz);

            // Determine if this ion sequence can under go a neutral loss
            // if it can, record the neutral loss variant
            record_NL_ions(y, Z, ym);
          }
        }
      }
    } // end iteration over charge state
  }


  // Function determines if the passed ion sequence can under go a neutral loss
  // If it can, the NL is recorded
  private void record_NL_ions(String ion, double z, double orig_ion_mass) {

    int start = ion.indexOf(":") + 1;
    int end = ion.length();
    if (ion.contains("/")) {
      end = ion.indexOf("/") - 1;
    }

    String x = ion.substring(start, end);

    if (!Globals.isDecoySeq(x)) {
      for (String nlr : Globals.nlMap.keySet()) {

        String candResidues = nlr.substring(0, nlr.indexOf("-"));
        int numCandRes = 0;

        for (int i = 0; i < x.length(); i++) {
          String residue = Character.toString(x.charAt(i));
          if (candResidues.contains(residue)) {
            numCandRes++;
          }
        }

        if (numCandRes > 0) { // this ion contains residues that can result in a NL
          double nl_mass = Globals.nlMap.get(nlr);
          String NL_tag = nlr.substring(nlr.indexOf("-"));

          String nl_str = ion + NL_tag;
          double mass = orig_ion_mass + nl_mass;
          double mz = Globals.round_dbl((mass / z), 4);

          if (mz > Globals.minMZ) {
            if (ion.startsWith("b")) {
              b_ions.put(nl_str, mz);
            }
            if (ion.startsWith("y")) {
              y_ions.put(nl_str, mz);
            }
          }
        }
      } // end loop Globals.nlMap
    } else {
      // The 'ion' string variable must contain at least 1 decoy residue modification
      // if you got to this point of the code.
      for (Entry<String, Double> e : Globals.decoyNLmap.entrySet()) {
        String NL_tag = e.getKey();
        double nl_mass = e.getValue();

        String nl_str = ion + NL_tag;
        double mass = orig_ion_mass + nl_mass;
        double mz = Globals.round_dbl((mass / z), 4);

        if (mz > Globals.minMZ) {
          if (ion.startsWith("b")) {
            b_ions.put(nl_str, mz);
          }
          if (ion.startsWith("y")) {
            y_ions.put(nl_str, mz);
          }
        }
      }
    }
  }


  private void calcNumDecoyPermutations() {

    int ctr = 0;
    for (int i = 0; i < pepLen; i++) {
      String aa = Character.toString(peptide.charAt(i));
      if (!Globals.targetModMap.containsKey(aa)) {
        ctr++;
      }
    }

    if (ctr >= numPPS) { // you have enough non-target reidues to make a decoy peptide
      numDecoyPermutations = StatsFunctions.combinatorial(ctr, numRPS);
    } else {
      numDecoyPermutations = 0;
    }
  }


  long getNumPerm() {
    return (numPermutations + numDecoyPermutations);
  }


  void printIons() {
    for (String p : b_ions.keySet()) {
      double m = Globals.round_dbl(b_ions.get(p), 4);
      System.out.println(modPeptide + "\t" + p + "\t" + m);
    }

    for (String p : y_ions.keySet()) {
      double m = Globals.round_dbl(y_ions.get(p), 4);
      System.out.println(modPeptide + "\t" + p + "\t" + m);
    }
  }


  /**
   * Returns a map where the key is each possible permutation
   * of the peptide and the value is the score this permutation receives later on.
   * @param permType 0 = generate forward (positive) permutations,
*           any other value = generate decoy (negative) permutations.
   * @return
   */
  THashMap<String, Double> getPermutations(int permType) {
    THashMap<String, Double> ret = new THashMap<>();
    TIntArrayList candModSites = new TIntArrayList();
    ArrayList<TIntArrayList> x;
    int i;

    if (permType == 0) { // generate forward (positive) permutations

      // Get candidate "true" sites to modify
      for (i = 0; i < pepLen; i++) {
        String aa = Character.toString(peptide.charAt(i));
        if (Globals.targetModMap.containsKey(aa)) {
          candModSites.add(i);
        }
      }

      // For the given candidate sites that can undergo modifications,
      // generate all possible permutations involving them.
      x = StatsFunctions.getAllCombinations(candModSites, numRPS);

      for (TIntList a : x) {
        StringBuilder modPep = new StringBuilder();
        if (modPosMap.containsKey(Constants.NTERM_MOD)) {
          modPep = new StringBuilder("[");
        }

        for (i = 0; i < pepLen; i++) {
          String aa = Character.toString(peptide.charAt(i)).toLowerCase();

          if (a.contains(i)) { // site to be modified
            modPep.append(aa); //Character.toString( Character.toLowerCase(peptide.charAt(i)) );
          } else if (nonTargetMods.containsKey(i)) {
            modPep.append(aa);
          } else {
            modPep.append(aa.toUpperCase());
          }
        }

        if (modPosMap.containsKey(Constants.CTERM_MOD)) {
          modPep = new StringBuilder("]");
        }

        ret.put(modPep.toString(), 0d);
      }
    } else { // generate decoy (negative) permutations

      for (i = 0; i < pepLen; i++) {
        String AA = Character.toString(peptide.charAt(i));
        String aa = AA.toLowerCase();
        int score = 0;
        if (!Globals.targetModMap.containsKey(AA)) {
          score++;
        }

        // TODO: XXX: BUG: this.nonTargetMods map never contains keys of type String
        if (!this.nonTargetMods.containsKey(aa)) {
          score++;
        }

        if (score == 2) {
          candModSites.add(i);
        }
      }

      // For the given candidate sites that can undergo modifications,
      // generate all possible permutations involving them.
      x = StatsFunctions.getAllCombinations(candModSites, numRPS);

      for (TIntList a : x) {
        StringBuilder modPep = new StringBuilder();
        if (modPosMap.containsKey(Constants.NTERM_MOD)) {
          modPep = new StringBuilder("[");
        }

        for (i = 0; i < pepLen; i++) {
          String aa = Character.toString(peptide.charAt(i)).toLowerCase();

          if (a.contains(i)) { // site to be modified
            String decoyChar = Globals.getDecoySymbol(peptide.charAt(i));
            modPep.append(decoyChar);
          } else if (nonTargetMods.containsKey(i)) {
            modPep.append(aa);
          } else {
            modPep.append(aa.toUpperCase());
          }
        }

        if (modPosMap.containsKey(Constants.CTERM_MOD)) {
          modPep = new StringBuilder("]");
        }
        ret.put(modPep.toString(), 0d);
      }
    }

    return ret;
  }


  // Function returns the ion ladder for this peptide sequence
  public THashMap<String, Double> getIonLadder() {
    THashMap<String, Double> ret = new THashMap<>();

    ret.putAll(b_ions);
    ret.putAll(y_ions);

    return ret;
  }


  // Function scores the current peptide
  void calcScore_CID() {

    if (matchedPeaks.isEmpty()) {
      score = 0;
    } else {
      ModelData_CID MD = Globals.modelingMap_CID.get(this.charge);

      // Now compute the scores for these peaks
      double intensityM = 0;
      double intensityU;
      double distM = 0;
      double distU;
      double Iscore;
      double Dscore;
      score = 0;

      for (PeakClass pk : matchedPeaks) {

        intensityU = StatsFunctions.log_gaussianProb(MD.mu_int_U, MD.var_int_U, pk.norm_intensity);
        distU = StatsFunctions.log_gaussianProb(MD.mu_dist_U, MD.var_dist_U, pk.dist);

        if (pk.matchedIonStr.startsWith("b")) {
          intensityM = StatsFunctions.log_gaussianProb(MD.mu_int_B, MD.var_int_B, pk.norm_intensity);
          distM = StatsFunctions.log_gaussianProb(MD.mu_dist_B, MD.var_dist_B, pk.dist);
        }

        if (pk.matchedIonStr.startsWith("y")) {
          intensityM = StatsFunctions.log_gaussianProb(MD.mu_int_Y, MD.var_int_Y, pk.norm_intensity);
          distM = StatsFunctions.log_gaussianProb(MD.mu_dist_Y, MD.var_dist_Y, pk.dist);
        }

        Iscore = intensityM - intensityU;
        Dscore = distM - distU;

        double intense_wt = 1.0 / (1.0 + FastMath.exp(-Iscore));
        double x;

        if (Double.isNaN(Dscore) || Double.isInfinite(Dscore)) {
          x = 0;
        } else {
          x = intense_wt * Dscore;
        }

        if (x < 0) {
          x = 0; // this prevents the score from going negative
        }

        pk.score = x;
        pk.distScore = Dscore;
        pk.intensityScore = Iscore;

        score += x;
      }
    }
  }


  // Function returns true if the modPeptide assigned to this object is a decoy
  boolean isDecoyPep() {
    for (int i = 0; i < pepLen; i++) {
      String aa = Character.toString(modPeptide.charAt(i));
      if (Globals.isDecoyResidue(aa)) {
        return true;
      }
    }

    return false;
  }


  /****************
   * Function computes the score for the HCD algorithm
   */
  void calcScore_HCD() {

    if (matchedPeaks.isEmpty()) {
      score = 0;
    } else {
      ModelData_HCD MD = Globals.modelingMap_HCD.get(this.charge);

      // we double the error window size for decoys. Otherwise they may not get matched peaks
      isDecoyPep();

      // Now compute the scores for these peaks
      double intensityM;
      double intensityU;
      double distM;
      double distU;
      double Iscore;
      double Dscore;
      score = 0;

      for (PeakClass pk : matchedPeaks) {
        intensityU = MD.getLogNPdensityInt('n', pk.norm_intensity);
        distU = 0; // log of uniform distribution between -1 and 1 is zero

        char ionType = pk.matchedIonStr.charAt(0);
        intensityM = MD.getLogNPdensityInt(ionType, pk.norm_intensity);
        distM = MD.getLogNPdensityDistPos(pk.dist);

        Iscore = intensityM - intensityU;
        Dscore = distM - distU;

        if (Double.isNaN(Iscore) || Double.isInfinite(Iscore)) {
          Iscore = 0;
        }
        if (Double.isNaN(Dscore) || Double.isInfinite(Dscore)) {
          Dscore = 0;
        }

        double x = Iscore + Dscore;
        if (x < 0) {
          x = 0; // this prevents the score from going negative
        }

        pk.score = x;
        pk.distScore = Dscore;
        pk.intensityScore = Iscore;

        score += x;
      }
    }
  }


  /**********
   * Function identifies all of the peaks in the passed spectrum that can be matched to the
   * theoretical ions for this peptide.
   */
  public void matchPeaks(SpectrumClass obsPeakList) {
//        if( !Globals.modelingMap_CID.containsKey(this.charge) ) {
//            System.err.println("\nError! " + peptide + "/+" + this.charge +
//                    ": a CID Model does not exist for this charge state!\nExiting now.\n");
//            System.exit(0);
//        }

    double matchErr;
    double a, b;

    matchedPeaks = new ArrayList<>();

    int N = b_ions.size() + y_ions.size();

    // The same observed peak might be matchable to multiple theoretical ions for this peptide.
    // This HashMap is used to keep track of the best possible match for an observed peak
    // k = observed peak's m/z value
    // v = the best match for this observed peak based upon the abs(mz_dist)
    TDoubleObjectHashMap<PeakClass> bestMatchMap = new TDoubleObjectHashMap<>(N);

    // Try to match y-ions first, they tend to be of higher intensity which is what we want
    for (String theo_ion : y_ions.keySet()) {
      double theo_mz = y_ions.get(theo_ion);

      if (Globals.ms2tol_units == Constants.PPM_UNITS) {
        double ppmErr = Globals.ms2tol / Constants.PPM;
        matchErr = theo_mz * ppmErr;
      } else {
        matchErr = Globals.ms2tol;
      }

      matchErr *= 0.5; // split in half
      a = Globals.round_dbl((theo_mz - matchErr), 4);
      b = Globals.round_dbl((theo_mz + matchErr), 4);

      // Hold candidate matching peaks here
      ArrayList<PeakClass> cand = new ArrayList<>();
      for (int i = 0; i < obsPeakList.N; i++) {
        double obsMZ = Globals.round_dbl(obsPeakList.mz[i], 4);

        if ((obsMZ >= a) && (obsMZ <= b)) {
          cand.add(obsPeakList.getPeakClassInstance(i));
        }
      }

      if (!cand.isEmpty()) {
        // Take the most intense peak as the best match
        cand.sort(PeakClass.comparator_intensity_hi2low);
        PeakClass pk = cand.get(0);

        pk.matched = true;
        pk.matchedIonStr = theo_ion;
        pk.dist = pk.mz - theo_mz;

        // Check to see if you have already assigned this peak to a theoretical ion
        if (bestMatchMap.containsKey(pk.mz)) {
          PeakClass oldPK = bestMatchMap.get(pk.mz);
          if (Math.abs(oldPK.dist) > (Math.abs(pk.dist))) {
            bestMatchMap.put(pk.mz, pk);
          }
        } else {
          bestMatchMap.put(pk.mz, pk);
        }
      }
      cand.clear();
    }

    // Now try to match b-ions
    for (String theo_ion : b_ions.keySet()) {
      double theo_mz = b_ions.get(theo_ion);

      if (Globals.ms2tol_units == Constants.PPM_UNITS) {
        double ppmErr = Globals.ms2tol / Constants.PPM;
        matchErr = theo_mz * ppmErr;
      } else {
        matchErr = Globals.ms2tol;
      }

      matchErr *= 0.5; // split in half
      a = Globals.round_dbl((theo_mz - matchErr), 4);
      b = Globals.round_dbl((theo_mz + matchErr), 4);

      // Hold candidate matching peaks here
      ArrayList<PeakClass> cand = new ArrayList<>();
      for (int i = 0; i < obsPeakList.N; i++) {
        double obsMZ = Globals.round_dbl(obsPeakList.mz[i], 4);

        if ((obsMZ >= a) && (obsMZ <= b)) {
          cand.add(obsPeakList.getPeakClassInstance(i));
        }
      }

      if (!cand.isEmpty()) {
        // Take the most intense peak as the best match
        cand.sort(PeakClass.comparator_intensity_hi2low);
        PeakClass pk = cand.get(0);

        pk.matched = true;
        pk.matchedIonStr = theo_ion;
        pk.dist = pk.mz - theo_mz;

        // Check to see if you have already assigned this peak to a theoretical ion
        if (bestMatchMap.containsKey(pk.mz)) {
          PeakClass oldPK = bestMatchMap.get(pk.mz);
          if (Math.abs(oldPK.dist) > (Math.abs(pk.dist))) {
            bestMatchMap.put(pk.mz, pk);
          }
        } else {
          bestMatchMap.put(pk.mz, pk);
        }
      }
      cand.clear();
    }

    // MM holds the best match for each observed peak.
    for (double mz : bestMatchMap.keys()) {
      matchedPeaks.add(bestMatchMap.get(mz));
    }
    bestMatchMap.clear();
  }
}
