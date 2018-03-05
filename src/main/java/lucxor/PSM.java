/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author dfermin
 */
class PSM {

  String specId; // TPP-based name
  String srcFile; // name of file from which spectrum is derived
  int scanNum;
  int charge;
  double PSMscore; // score assigned to this PSM by original search engine
  double deltaScore;
  double localFDR;  // local false localization rate
  double globalFDR; // global false localization rate
  boolean isKeeper; // true means the PSM can be scored
  boolean useForModel; // true means the PSM can be used for modeling
  boolean isDecoy; // true means the top match for this PSM is a decoy
  private boolean isUnambiguous; // true means the number potential PTM sites is equal to the number of reported PTM sites

  TIntDoubleHashMap modCoordMap; // holds modified amino acid positions

  private THashMap<String, Double> posPermutationScoreMap = null;
  private THashMap<String, Double> negPermutationScoreMap = null;

  // matched and unmatched peaks
  ArrayList<PeakClass> posPeaks = null;
  ArrayList<PeakClass> negPeaks = null;

  private SpectrumClass PeakList;


  final Peptide origPep;
  private Peptide score1pep; // top scoring peptide permutation
  private Peptide score2pep; // 2nd best scoring peptide permutation


  // default constructor
  PSM() {
    origPep = new Peptide();
    modCoordMap = new TIntDoubleHashMap(3);
    isKeeper = false;
    useForModel = false;
    isDecoy = false;
    isUnambiguous = false;
    PeakList = new SpectrumClass();
  }

  void init() {


    // TODO: Continue here, check what all these initializations do
    origPep.initialize(modCoordMap);
    origPep.charge = this.charge;


    if (origPep.numPPS == origPep.numRPS) {
      isUnambiguous = true;
    }

    // Determine if the PSM is to be kept and if it should be used for modeling
    if (origPep.numPPS > 0 && origPep.numRPS > 0 &&
        PSMscore >= Globals.scoreTH && charge <= Globals.maxChargeState &&
        origPep.pepLen <= Globals.maxPepLen)
    {
      isKeeper = true;
    }

    if (isKeeper && (PSMscore >= Globals.modelTH)) {
      useForModel = true;
    }

    if (isKeeper) {
      // We will reconstruct the specId here to remove left-padded zeros
      String suffix = "";

      if (Globals.inputType == Constants.PEPXML) {
        if (Globals.spectrumSuffix.equalsIgnoreCase("mzXML")) {
          suffix = "mzXML";
        }
        if (Globals.spectrumSuffix.equalsIgnoreCase("mzML")) {
          suffix = "mzML";
        }
        if (Globals.spectrumSuffix.equalsIgnoreCase("mgf")) {
          suffix = "mgf";
        }

        /*
         ** Some scans might have a name formatted like this:
         ** expt_file_name.c.14033.14033.2
         ** we need to be able to handle these cases
         */
        Pattern r = Pattern.compile("(.+)\\.\\d+\\.\\d+\\.\\d$"); // pattern of a spectrum ID
        Matcher m = r.matcher(specId);

        if (m.find()) {
          srcFile = m.group(1) + "." + suffix;
          specId = m.group(1) + "." +
              Integer.toString(scanNum) + "." +
              Integer.toString(scanNum) + "." +
              Integer.toString(charge);
        } else {
          System.err.println("\nERROR: PSM.java::Unable to correctly format specId variable.\n");
          System.exit(0);
        }
      } else {
        // We need to construct a valid TPP-specId variable for
        // this PSM instance because it was read in from a tab-delimited file
        Pattern r = Pattern.compile("(.+)\\.(mgf|mzXML|mzML)$"); // pattern of a spectrum ID
        Matcher m = r.matcher(srcFile);

        if (m.find()) {
          specId = m.group(1) + "." +
              Integer.toString(scanNum) + "." +
              Integer.toString(scanNum) + "." +
              Integer.toString(charge);
        } else {
          System.err.println("\nERROR: PSM.java::Unable to correctly format specId variable.\n");
          System.exit(0);
        }
      }

      origPep.build_ion_ladders();
    }

  }


  void recordSpectra(SpectrumClass S) {

    PeakList = S;

    // reduce the impact of the neutral loss peak (if any)
    if (Globals.reduceNL == 1) {
      reduceNLpeak();
    }

    // normalize the peak intensities to the median
    PeakList.medianNormalizeSpectra();
  }


  // If the user asked for it, reduce the intensity of the precursor neutral loss peak
  private void reduceNLpeak() {
    double pepMHplus = Globals.getFragmentIonMass(origPep.modPeptide, 1.0, Constants.WATER);
    double NLmass = MathHelper.roundDouble((pepMHplus + Globals.precursorNLmass),
        3); // the Globals.precursorNLmass is a negative number
    double NLmz = MathHelper.roundDouble((NLmass / (double) charge), 3);

    // Do not perform this NL reduction if the peptide doesn't contain
    // a serine or theronine
    if (!origPep.peptide.contains("S") && !origPep.peptide.contains("T")) {
      return;
    }

    // Find the most intense peak in the spectrum
    double maxPk_mz = PeakList.mz[PeakList.maxI_index];

    double mzDelta = MathHelper.roundDouble(Math.abs((maxPk_mz - NLmz)), 5);
    double mzErr = Constants.PROTON;
    if (mzDelta <= mzErr) { // the max peak is a neutral loss peak

      double orig_maxI = PeakList.maxI; // record the original intensity of the peak
      int orig_maxI_index = PeakList.maxI_index;

      PeakList.raw_intensity[PeakList.maxI_index] = 0; // bring the peak down to zero

      // identify the second-most intense peak in the spectrum
      double new_maxI = 0d;
      for (int i = 0; i < PeakList.N; i++) {
        if (PeakList.raw_intensity[i] > new_maxI) {
          new_maxI = PeakList.raw_intensity[i];
        }
      }
      PeakList.maxI = new_maxI;

      PeakList.calcRelativeIntensity(); // recompute the relative peak intensities
      PeakList.raw_intensity[orig_maxI_index] = orig_maxI;
      PeakList.rel_intensity[orig_maxI_index] = 50.0; // set the peak intensity to be 50%

//          PeakList.maxI *= 0.5; // reduce the peak intensity by 50%
    }

  }


  // Function returns a map of the modified residues in the given modPeptide
  private TIntDoubleHashMap getModMap(String modPep) {
    TIntDoubleHashMap ret = new TIntDoubleHashMap();

    // record n-terminal modification (if any)
    if (modCoordMap.containsKey(Constants.NTERM_MOD)) {
      ret.put(Constants.NTERM_MOD, Globals.ntermMass);
    }

    // record c-terminal modification (if any)
    if (modCoordMap.containsKey(Constants.CTERM_MOD)) {
      ret.put(Constants.CTERM_MOD, Globals.ntermMass);
    }

    // record all variable modifications
    int startAt = 0;

    for (int i = startAt; i < modPep.length(); i++) {
      double mass;
      char c = modPep.charAt(i);

      // record decoy-modified residues
      if (Globals.isDecoyResidue(Character.toString(c))) {
        mass = Globals.AAmassMap.get(Character.toString(c));
        ret.put(i, mass);
      }
      // record normal-modified residues
      else if (Character.isLowerCase(c)) {
        mass = Globals.AAmassMap.get(Character.toString(c));
        ret.put(i, mass);
      }
    }

    // Search the original sequence assigned to this PSM to identify
    // any residues that are already modified for this peptide. These
    // are residues that cannot be used for building decoys
    for (int p : origPep.nonTargetMods.keySet()) {
      if ((p == Constants.NTERM_MOD) || (p == Constants.CTERM_MOD)) {
        continue;
      }

      String aa = origPep.nonTargetMods.get(p);
      double mass = Globals.AAmassMap.get(aa);
      ret.put(p, mass);
    }

    return ret;
  }


  // This function identifies all of the peaks in PeakList that can be matched to any fragment ion of any permutation
  // of the peptide sequence associated with this PSM
  void matchAllPeaks() {

    posPeaks = new ArrayList<>(PeakList.N);

    TIntDoubleHashMap mcp = new TIntDoubleHashMap(); // use this to record modifications for each permutation

    // This variable will hold all of the fragment ions that are possible for every
    // permutation of the peptide sequence assigned to this peptide
    THashMap<String, Double> all_ions = new THashMap<>();

    for (String pep : posPermutationScoreMap.keySet()) {
      mcp.clear();
      mcp = getModMap(pep); // get the modified residues for this sequence permutation
      Peptide curPep = new Peptide();
      curPep.initialize(origPep.peptide, pep, charge, mcp);
      curPep.build_ion_ladders();

      // get the ion ladder for this permutation and match it against the
      // spectrum assigned to this PSM
      all_ions.putAll(curPep.getIonLadder());
    }

    // The same observed peak could be assigned to multiple theoretical peaks.
    // We use this data object to collect all of these cases for later processing.
    // k = observed Peak MZ value
    // v = array list of all the theoretical peaks that match to it
    TDoubleObjectHashMap<ArrayList<PeakClass>> candMatchedPks = new TDoubleObjectHashMap<>(
        PeakList.N);

    // Iterate over each theoretical ion
    for (String theo_ion : all_ions.keySet()) {

      double theo_mz = all_ions.get(theo_ion);

      PeakClass matchedPk = getMatchedPeak(theo_ion, theo_mz);
      if (null == matchedPk) {
        continue; // no match was found for this theoretical peak
      }

      if (candMatchedPks.containsKey(matchedPk.mz)) {
        ArrayList<PeakClass> ary = candMatchedPks.get(matchedPk.mz);
        PeakClass newEntry = new PeakClass(matchedPk);
        ary.add(newEntry);
        candMatchedPks.put(matchedPk.mz, ary);
      } else {
        ArrayList<PeakClass> ary = new ArrayList<>();
        ary.add(matchedPk);
        candMatchedPks.put(matchedPk.mz, ary);
      }
    }

    // For each matched peak in 'candMatchedPks' now select the match that is the closest by m/z value
    for (double mz : candMatchedPks.keys()) {
      ArrayList<PeakClass> ary = candMatchedPks.get(mz);
      ary.sort(PeakClass.comparator_mz_abs_dist);
      posPeaks.add(ary.get(0));
    }
    candMatchedPks.clear();

    // This can happen if there are no matched peaks
    // identified for this spectrum due to the exclusion of neutral loss ions
    // or due to a very small fragment ion tolerance window.
    // When it happens we just exit the function.
    if (posPeaks.isEmpty()) {
      useForModel = false;
      posPeaks.clear();
      return;
    }

    // Now that you have matched all the peaks you can for this spectrum
    // All the remaining peaks are noise, we put them in the negPeaks arraylist
    negPeaks = new ArrayList<>(PeakList.N);

    for (int i = 0; i < PeakList.N; i++) {

      PeakClass negPk = PeakList.getPeakClassInstance(i);

      if (!posPeaks.contains(negPk)) {

        double finalDist;

        // record the distance of this unmatched peak to the nearest peak in 'posPeaks'
        ArrayList<Double> dist = new ArrayList<>(posPeaks.size());
        ArrayList<Double> absDist = new ArrayList<>(posPeaks.size());

        int j = 0;
        for (PeakClass matchedPk : posPeaks) {
          double d = matchedPk.mz - negPk.mz;
          dist.add(d);
          absDist.add(Math.abs(d));
          j++;
        }

        Collections.shuffle(dist);
        Collections.sort(absDist); // sort the m/z distances from low to high
        finalDist = absDist.get(0);
        for (double d : dist) {
          if (Math.abs(d) == finalDist) {
            negPk.matched = false;
            negPk.dist = d;
            negPeaks.add(negPk);
            break; // leave loop
          }
        }
      }
    }
  }


  /*******
   * This function returns a PeakClass object matched to the given theoretical mz value.
   * The function returns a null object if no match was found.
   */
  private PeakClass getMatchedPeak(String theo_ion, double theo_mz) {

    PeakClass ret = null;

    double matchErr;
    double a, b;

    // Compute the fragment error tolerance that will be used
    if (Globals.ms2tol_units == Constants.PPM_UNITS) {
      double ppmErr = Globals.ms2tol / Constants.PPM;
      matchErr = theo_mz * ppmErr;
    } else {
      matchErr = Globals.ms2tol;
    }

    matchErr *= 0.5; // split in half

    a = theo_mz - matchErr;
    b = theo_mz + matchErr;

    // Iterate over the peaks in PeakList
    ArrayList<PeakClass> candMatches = new ArrayList<>();
    for (int i = 0; i < PeakList.N; i++) {

      if ((PeakList.mz[i] >= a) && (PeakList.mz[i] <= b)) {
        PeakClass pk = PeakList.getPeakClassInstance(i);
        candMatches.add(pk);
      }
    }

    // At least one match was found for this theoretical mz value
    if (!candMatches.isEmpty()) {

      // identify the most intense peak in candMatches
      candMatches.sort(PeakClass.comparator_intensity_hi2low);

      ret = candMatches.get(0);
      ret.matched = true;
      ret.dist = ret.mz - theo_mz; // obs - expected
      ret.matchedIonStr = theo_ion;
      ret.matchedIonMZ = theo_mz;
    }

    return ret;
  }


  /**
   * Creates all of the permutations for the sequence assigned to this PSM.
   * @param RN
   */
  void generatePermutations(int RN) {
    posPermutationScoreMap = origPep.getPermutations(0);

    if (!isUnambiguous) {
      if (RN == 0) { // RN = runMode is selected to create decoys
        negPermutationScoreMap = origPep.getPermutations(1);
      } else {
        negPermutationScoreMap = new THashMap<>(0);
      }
    }
  }

  /**
   * Function to fill in required variables when a thread instance of this class is killed due to
   * either running over time or when the executor thread pool is shutdown.
   */
  private void killThreadResults() {
    deltaScore = -1;
    score1pep = origPep;
    score2pep = origPep;
  }


  // Function scores every candidate sequence associated with this PSM
  void scorePermutations() throws IOException {

    if (Thread.interrupted()) { // this code just prevents excessive run times
      killThreadResults();
      return;
    }

    File debugF = null;
    if (Globals.debugMode == Constants.WRITE_PERM_SCORES) {
      debugF = new File("all_scores.debug");
    } else if (Globals.debugMode == Constants.WRITE_ALL_MATCHED_PK_SCORES) {
      debugF = new File("all_matched_pks.debug");
    }

    BufferedWriter bw = null;
    TIntDoubleHashMap mcp;
    try {
      if (debugF != null) {
        boolean existed = debugF.exists();
        bw = new BufferedWriter(new FileWriter(debugF, existed));
        if (!existed) {
          if (Globals.debugMode == Constants.WRITE_PERM_SCORES) {
            bw.write("specId\torigModPep\tcurPermutation\tisDecoy\tscore\n");
          } else if (Globals.debugMode == Constants.WRITE_ALL_MATCHED_PK_SCORES) {
            bw.write(
                "specId\tcurPermutation\tisDecoy\tionSeq\tmz\tintensity\tmzDist\tIscore\tDscore\tscore\n");
          }
        }
      }

      // Scores for forward sequences
      String line;
      for (String curSeq : posPermutationScoreMap.keySet()) {

        mcp = this.getModMap(curSeq);

        Peptide curPep = new Peptide();
        curPep.initialize(origPep.peptide, curSeq, charge, mcp);
        curPep.build_ion_ladders();
        curPep.numPPS = origPep.numPPS;
        curPep.numRPS = origPep.numRPS;

        curPep.matchPeaks(PeakList); // match all the peaks you can for this peptide permutation

        if (Globals.scoringAlgorithm == Constants.CID) {
          curPep.calcScore_CID();
        }
        if (Globals.scoringAlgorithm == Constants.HCD) {
          curPep.calcScore_HCD();
        }

        if (Globals.debugMode == Constants.WRITE_PERM_SCORES && bw != null) {
          line = specId + "\t" + origPep.modPeptide + "\t" +
              curPep.modPeptide + "\t0\t" + curPep.score + "\n";
          bw.write(line);
        }

        writeAllMatchedPkScores(bw, curSeq, curPep);

        posPermutationScoreMap.put(curSeq, curPep.score);
      }

      if (!isUnambiguous) { // at least 2 possible permutations exist for this PSM

        // Scores for decoy sequences
        for (String curSeq : negPermutationScoreMap.keySet()) {
          Peptide curPep = new Peptide();
          mcp = this.getModMap(curSeq);
          curPep.initialize(origPep.modPeptide, curSeq, charge, mcp);
          curPep.build_ion_ladders();
          curPep.numPPS = origPep.numPPS;
          curPep.numRPS = origPep.numRPS;

          curPep.matchPeaks(PeakList);

          if (Globals.scoringAlgorithm == Constants.CID) {
            curPep.calcScore_CID();
          }
          if (Globals.scoringAlgorithm == Constants.HCD) {
            curPep.calcScore_HCD();
          }

          if (Globals.debugMode == Constants.WRITE_PERM_SCORES && bw != null) {
            line = specId + "\t" + origPep.modPeptide + "\t" +
                curPep.modPeptide + "\t1\t" + curPep.score + "\n";
            bw.write(line);
          }

          writeAllMatchedPkScores(bw, curSeq, curPep);

          negPermutationScoreMap.put(curSeq, curPep.score);
        }
      }

    } finally {
      if (bw != null) {
        bw.close();
      }
    }
    // *************  debug writer closed here

    // collect all of the scores into a single arraylist and select the top
    // two matches to compute the delta score
    ArrayList<Double> allScores = new ArrayList<>();
    for (String curSeq : posPermutationScoreMap.keySet()) {
      double d = posPermutationScoreMap.get(curSeq);
      allScores.add(d);
    }

    if (!isUnambiguous) {
      for (String curSeq : negPermutationScoreMap.keySet()) {
        double d = negPermutationScoreMap.get(curSeq);
        allScores.add(d);
      }
    }

    Collections.sort(allScores); // low to high
    Collections.reverse(allScores);// high to low

    double score1 = MathHelper.roundDouble(allScores.get(0), 6);
    double score2 = 0;
    if (!isUnambiguous) {
      score2 = MathHelper.roundDouble(allScores.get(1), 6);
    }

    String pep1 = "";
    String pep2 = "";
    int numAssigned = 0;

    if (!isUnambiguous) {
      // Find the permutations that have the values stored in score1 and score2
      // We'll first try to find the top scores among the non-decoy
      // permutations
      for (Entry<String, Double> e : posPermutationScoreMap.entrySet()) {
        String curSeq = e.getKey();
        double x = e.getValue();
        double d = MathHelper.roundDouble(x, 6);

        if ((d == score1) && (pep1.isEmpty())) {
          pep1 = curSeq;
          numAssigned++;
        } else if ((d == score2) && (pep2.isEmpty())) {
          pep2 = curSeq;
          numAssigned++;
        }

        if (numAssigned == 2) {
          break;
        }
      }

      // if this is true, then you need to search among the decoys
      if (numAssigned != 2) {
        for (Entry<String, Double> e : negPermutationScoreMap.entrySet()) {
          String curSeq = e.getKey();
          double x = e.getValue();
          double d = MathHelper.roundDouble(x, 6);

          if ((d == score1) && (pep1.isEmpty())) {
            pep1 = curSeq;
            numAssigned++;
          } else if ((d == score2) && (pep2.isEmpty())) {
            pep2 = curSeq;
            numAssigned++;
          }

          if (numAssigned == 2) {
            break;
          }
        }
      }
    } else { // special case for unambiguous PSMs
      for (Entry<String, Double> e : posPermutationScoreMap.entrySet()) {
        pep1 = e.getKey();
      }
    }

    // Record the best match into score1
    mcp = this.getModMap(pep1);
    score1pep = new Peptide();
    score1pep.initialize(origPep.peptide, pep1, charge, mcp);
    score1pep.score = score1;
    isDecoy = score1pep.isDecoyPep(); // determine if this PSM is a decoy hit

    if (isUnambiguous) {
      // This is an unambiguous case.
      // Make a duplicate of score1pep and assign it to score2pep.
      //mcp = new HashMap();
      mcp = this.getModMap(pep1);
      score2pep = new Peptide();
      score2pep.initialize(origPep.peptide, pep1, charge, mcp);
      score2pep.score = 0;
      deltaScore = score1;
    } else {
      // Record second best match
      //mcp = new HashMap();
      mcp = this.getModMap(pep2);
      score2pep = new Peptide();
      score2pep.initialize(origPep.peptide, pep2, charge, mcp);
      score2pep.score = score2;
      deltaScore = score1 - score2;
    }
  }

  private void writeAllMatchedPkScores(BufferedWriter bw, String curSeq, Peptide curPep)
      throws IOException {
    String line;
    if (Globals.debugMode == Constants.WRITE_ALL_MATCHED_PK_SCORES) {
      for (PeakClass pk : curPep.matchedPeaks) {
        line = specId + "\t" + curSeq + "\t" + Globals.isDecoySeq(curSeq) + "\t" +
            pk.matchedIonStr + "\t" + pk.mz + "\t" +
            pk.rel_intensity + "\t" + pk.dist + "\t" +
            pk.intensityScore + "\t" + pk.distScore + "\t" + pk.score + "\n";
        bw.write(line);
      }
      bw.write("\n");
    }
  }


  // Function returns results as a string
  String getResults() {

    String ret = specId + "\t";

    if (Globals.peptideRepresentation == Constants.SINGLE_CHAR) {
      ret += origPep.peptide + "\t";
      ret += score1pep.modPeptide + "\t";
      ret += score2pep.modPeptide + "\t";

    } else {
      ret += origPep.getModPepTPP() + "\t";
      ret += score1pep.getModPepTPP() + "\t";
      ret += score2pep.getModPepTPP() + "\t";
    }

    ret += Integer.toString(origPep.numPPS) + "\t" + Integer.toString(origPep.numRPS) + "\t";

    ret += Double.toString(PSMscore) + "\t";

    if (Globals.runMode == Constants.REPORT_DECOYS) {
      ret += (score1pep.isDecoyPep() ? 1 : 0) + "\t" + (score2pep.isDecoyPep() ? 1 : 0) + "\t";
    }

    ret += Double.toString(deltaScore) + "\t";

    ret += Double.toString(score1pep.score) + "\t"
        + Double.toString(score2pep.score) + "\t";

    ret += Double.toString(globalFDR) + "\t" + Double.toString(localFDR);

    ret += "\n";

    return ret;
  }


  // This function will write the scored peaks for the top permutation assigned
  // to this PSM. The data is written in a TAB-delimited format.
  void debug_writeScoredPeaks() throws IOException {

    String tsvFileName;
    if (Globals.inputType == Constants.PEPXML) {
      tsvFileName = specId + ".tsv";
    } else {
      int i = srcFile.lastIndexOf(Globals.spectrumSuffix);
      String k = String.format("%05d", scanNum) + Integer.toString(charge);
      tsvFileName = srcFile.substring(0, i) + k + ".tsv";
    }

    // With the next 4 lines you create output directory if necessary
    String outDirName = "debug_scored_peaks." + Globals.dateStamp;
    String outFileStr = outDirName + "/" + tsvFileName;

    File outF = new File(outFileStr);
    outF.getParentFile().mkdirs();

    FileWriter fw = new FileWriter(outF.getAbsoluteFile());
    BufferedWriter bw = new BufferedWriter(fw);

    bw.write("num\tmodPeptide\tfragmentIon\tobs m/z\tmatchDist\trelIntensity\t" +
        "normIntensity\tDscore\tIscore\tfinalScore\n");

    score1pep.build_ion_ladders();
    score1pep.matchPeaks(PeakList);

    if (Globals.scoringAlgorithm == Constants.CID) {
      score1pep.calcScore_CID();
    }
    if (Globals.scoringAlgorithm == Constants.HCD) {
      score1pep.calcScore_HCD();
    }

    ArrayList<PeakClass> mPK = score1pep.matchedPeaks;
    mPK.sort(PeakClass.comparator_mz);
    for (PeakClass pk : mPK) {
      if (!pk.matched) {
        continue;
      }

      String mz = Double.toString(MathHelper.roundDouble(pk.mz, 4));
      String relI = Double.toString(MathHelper.roundDouble(pk.rel_intensity, 4));
      String normI = Double.toString(MathHelper.roundDouble(pk.norm_intensity, 4));
      String Iscore = Double.toString(MathHelper.roundDouble(pk.intensityScore, 4));
      String Dscore = Double.toString(MathHelper.roundDouble(pk.distScore, 4));
      String score = Double.toString(MathHelper.roundDouble(pk.score, 4));
      String dist = Double.toString(MathHelper.roundDouble(pk.dist, 4));

      bw.write(
          "0\t" + score1pep.modPeptide + "\t" + pk.matchedIonStr + "\t" + mz + "\t" + dist + "\t" +
              relI + "\t" + normI + "\t" + Dscore + "\t" + Iscore + "\t" +
              score + "\n");
    }

    score2pep.build_ion_ladders();
    score2pep.matchPeaks(PeakList);
    if (Globals.scoringAlgorithm == Constants.CID) {
      score2pep.calcScore_CID();
    }
    if (Globals.scoringAlgorithm == Constants.HCD) {
      score2pep.calcScore_HCD();
    }

    mPK.clear();
    mPK = score2pep.matchedPeaks;
    mPK.sort(PeakClass.comparator_mz);
    for (PeakClass pk : mPK) {
      if (!pk.matched) {
        continue;
      }

      String mz = Double.toString(MathHelper.roundDouble(pk.mz, 4));
      String relI = Double.toString(MathHelper.roundDouble(pk.rel_intensity, 4));
      String normI = Double.toString(MathHelper.roundDouble(pk.norm_intensity, 4));
      String Iscore = Double.toString(MathHelper.roundDouble(pk.intensityScore, 4));
      String Dscore = Double.toString(MathHelper.roundDouble(pk.distScore, 4));
      String score = Double.toString(MathHelper.roundDouble(pk.score, 4));
      String dist = Double.toString(MathHelper.roundDouble(pk.dist, 4));

      bw.write(
          "1\t" + score2pep.modPeptide + "\t" + pk.matchedIonStr + "\t" + mz + "\t" + dist + "\t" +
              relI + "\t" + normI + "\t" + Dscore + "\t" + Iscore + "\t" +
              score + "\n");
    }

    bw.close();
  }


  // This function returns the matched peaks for the top 2 permutations assigned
  // to this PSM. The data is returned in a TAB-delimited format.
  String writeMatchedPks() {

    StringBuilder ret = new StringBuilder();

    score1pep.build_ion_ladders();
    score1pep.matchPeaks(PeakList);

    if (Globals.scoringAlgorithm == Constants.CID) {
      score1pep.calcScore_CID();
    }
    if (Globals.scoringAlgorithm == Constants.HCD) {
      score1pep.calcScore_HCD();
    }

    ArrayList<PeakClass> mPK = score1pep.matchedPeaks;
    mPK.sort(PeakClass.comparator_mz);
    for (PeakClass pk : mPK) {
      if (!pk.matched) {
        continue;
      }

      String mz = Double.toString(MathHelper.roundDouble(pk.mz, 4));
      String relI = Double.toString(MathHelper.roundDouble(pk.rel_intensity, 4));
      String Iscore = Double.toString(MathHelper.roundDouble(pk.intensityScore, 4));
      String Dscore = Double.toString(MathHelper.roundDouble(pk.distScore, 4));
      String score = Double.toString(MathHelper.roundDouble(pk.score, 4));

      ret.append(specId).append("\t1\t").append(score1pep.modPeptide).append("\t")
          .append(pk.matchedIonStr).append("\t").append(mz).append("\t").append(relI).append("\t")
          .append(Dscore).append("\t").append(Iscore).append("\t").append(score).append("\n");
    }

    score2pep.build_ion_ladders();
    score2pep.matchPeaks(PeakList);
    if (Globals.scoringAlgorithm == Constants.CID) {
      score2pep.calcScore_CID();
    }
    if (Globals.scoringAlgorithm == Constants.HCD) {
      score2pep.calcScore_HCD();
    }

    mPK.clear();
    mPK = score2pep.matchedPeaks;
    mPK.sort(PeakClass.comparator_mz);
    for (PeakClass pk : mPK) {
      if (!pk.matched) {
        continue;
      }

      String mz = Double.toString(MathHelper.roundDouble(pk.mz, 4));
      String relI = Double.toString(MathHelper.roundDouble(pk.rel_intensity, 4));
      String Iscore = Double.toString(MathHelper.roundDouble(pk.intensityScore, 4));
      String Dscore = Double.toString(MathHelper.roundDouble(pk.distScore, 4));
      String score = Double.toString(MathHelper.roundDouble(pk.score, 4));

      ret.append(specId).append("\t2\t").append(score2pep.modPeptide).append("\t")
          .append(pk.matchedIonStr).append("\t").append(mz).append("\t").append(relI).append("\t")
          .append(Dscore).append("\t").append(Iscore).append("\t").append(score).append("\n");
    }

    return (ret.toString());
  }


  // This function prepares the PSM to be scored again *AFTER* the FLR has been estimated.
  public void clearScores() {

    deltaScore = 0;
    globalFDR = Double.NaN;
    localFDR = Double.NaN;
    isDecoy = false;
    score1pep = null;
    score2pep = null;
    posPermutationScoreMap = null;
    negPermutationScoreMap = null;

  }

}
