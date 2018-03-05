/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;
import umich.ms.fileio.exceptions.FileParsingException;

/**
 * @author dfermin
 */
class Globals {

  private static File spectrumPath = null;
  static String spectrumSuffix = null;
  static String matchedPkFile = null;
  static File inputFile = null;
  static String outputFile = null;
  static String timeStamp = null;
  static String dateStamp = null;
  static int ms2tol_units;
  static int inputType;
  static int debugMode;
  static int reduceNL;
  static int peptideRepresentation;
  static int scoringMethod;
  static int scoringAlgorithm;
  static int maxChargeState;
  static int minNumPSMsForModeling;
  static int numThreads;
  static int runMode;
  static int tsvHdr;
  static int maxPepLen;
  static double ms2tol;
  static double modelTH;
  static double scoreTH;
  private static double decoyMass;
  static double minMZ;
  static double max_num_permutations;
  static double precursorNLmass;
  static double ntermMassProt;
  static double ntermMass;
  static double ctermMassProt;
  static double ctermMass;
  static boolean writeMatchedPeaks;

  static ArrayList<PSM> PSM_list = null;

  /** Mods user wants to search for. */
  static THashMap<String, Double> targetModMap = null;
  /** Fixed mods observed in data. */
  static THashMap<String, Double> fixedModMap = null;
  /** Variable mods observed in data. */
  static THashMap<String, Double> varModMap = null;
  /** Holds all neutral loss masses, k= list of amino acids v= NL mass. */
  static THashMap<String, Double> nlMap = null;
  static THashMap<String, Double> decoyNLmap = null;
  private static THashMap<Double, double[]> FLRestimateMap = null; // ary[0] = globalFLR, ary[1] = localFLR

  // TODO: ACHTUNG: XXX: Delete this monstrosity (replace with Char->Stirng or an array?)
  static THashMap<String, Double> AAmassMap = null;

  // TODO: ACHTUNG: XXX: Delete this monstrosity (replace with Char->Stirng or an array?)
  private static THashMap<Character, String> decoyAAMap = null;

  static THashMap<Integer, ModelData_CID> modelingMap_CID = null;
  static THashMap<Integer, ModelData_HCD> modelingMap_HCD = null;


  static void initialize() {
    PSM_list = new ArrayList<>();

    decoyAAMap = new THashMap<>();
    AAmassMap = new THashMap<>();
    targetModMap = new THashMap<>();
    fixedModMap = new THashMap<>();
    varModMap = new THashMap<>();
    nlMap = new THashMap<>();
    decoyNLmap = new THashMap<>();

    ntermMass = 0d;
    ntermMassProt = 0d;
    ctermMass = 0d;
    ctermMassProt = 0d;

    // in case we need it, initialize the timestamp signature variable
    java.util.Date date = new java.util.Date();
    SimpleDateFormat sdf = new SimpleDateFormat("yyyyMMMdd-hh_mm_ss");
    timeStamp = sdf.format(date);
    sdf = new SimpleDateFormat("yyyMMMdd");
    dateStamp = sdf.format(date);

    AAmassMap.put("A", 71.03711);
    AAmassMap.put("R", 156.10111);
    AAmassMap.put("N", 114.04293);
    AAmassMap.put("D", 115.02694);
    AAmassMap.put("C", 103.00919);
    AAmassMap.put("E", 129.04259);
    AAmassMap.put("Q", 128.05858);
    AAmassMap.put("G", 57.02146);
    AAmassMap.put("H", 137.05891);
    AAmassMap.put("I", 113.08406);
    AAmassMap.put("L", 113.08406);
    AAmassMap.put("K", 128.09496);
    AAmassMap.put("M", 131.04049);
    AAmassMap.put("F", 147.06841);
    AAmassMap.put("P", 97.05276);
    AAmassMap.put("S", 87.03203);
    AAmassMap.put("T", 101.04768);
    AAmassMap.put("W", 186.07931);
    AAmassMap.put("Y", 163.06333);
    AAmassMap.put("V", 99.06841);

    decoyAAMap.put('2', "A");
    decoyAAMap.put('3', "R");
    decoyAAMap.put('4', "N");
    decoyAAMap.put('5', "D");
    decoyAAMap.put('6', "C");
    decoyAAMap.put('7', "E");
    decoyAAMap.put('8', "Q");
    decoyAAMap.put('9', "G");
    decoyAAMap.put('0', "H");
    decoyAAMap.put('@', "I");
    decoyAAMap.put('#', "L");
    decoyAAMap.put('$', "K");
    decoyAAMap.put('%', "M");
    decoyAAMap.put('&', "F");
    decoyAAMap.put(';', "P");
    decoyAAMap.put('?', "W");
    decoyAAMap.put('~', "V");
    decoyAAMap.put('^', "S");
    decoyAAMap.put('*', "T");
    decoyAAMap.put('=', "Y");
  }


  // Function to incorporate information read in from input file into the
  // AAmassMap object
  public static void loadUserMods() {

    // Assign the information you got from the input file to the
    // relevant map. If none were given, then this should just proceed
    // without changing anything.
    for (String c : targetModMap.keySet()) {
      double mass = AAmassMap.get(c) + targetModMap.get(c);
      String symbol = c.toLowerCase();
      AAmassMap.put(symbol, mass);
    }

    for (String c : fixedModMap.keySet()) {
      double mass = AAmassMap.get(c) + fixedModMap.get(c);
      String symbol = c.toUpperCase();
      AAmassMap.put(symbol, mass);
    }

    for (String c : varModMap.keySet()) { // this will be a lower case key

      if (c.equalsIgnoreCase("[")) {
        Globals.ntermMass = varModMap.get(c);
      } else if (c.equalsIgnoreCase("]")) {
        Globals.ctermMass = varModMap.get(c);
      } else {
        double mass = AAmassMap.get(c.toUpperCase()) + varModMap.get(c);
        AAmassMap.put(c, mass);
      }
    }

    // Now add the decoy masses to the AAmassMap
    for (Character c : decoyAAMap.keySet()) {
      String trueAA = decoyAAMap.get(c);

      if (varModMap.containsKey(trueAA)) {
        continue;
      }
      if (targetModMap.containsKey(trueAA)) {
        continue;
      }

      double mass = AAmassMap.get(trueAA) + Globals.decoyMass;
      AAmassMap.put(c.toString(), mass);
    }

  }


  /**
   * This function is only called if we are getting our modifications from a pepXML file.
   */
  static void recordModsFromPepXML() {

    String alphabet = PepXML.AA_CHARS_UPPER;

    for (String c : fixedModMap.keySet()) {

      if (!alphabet.contains(c)) {
        continue; // skip non-amino acid characters
      }

      double mass = AAmassMap.get(c) + fixedModMap.get(c);
      String symbol = c.toUpperCase();
      AAmassMap.put(symbol, mass);
    }

    // First filter the data in varModMap.
    // We want to remove non-standard amino acid characters and remove
    // the amino acids that are in our 'targetModMap' variable.
    HashMap<String, Double> tmp = new HashMap<>(varModMap);
    varModMap.clear();

    for (String c : tmp.keySet()) {
      String C = c.toUpperCase();
      double mass = tmp.get(c);
      if (!alphabet.contains(C)) {
        continue; // skip non-amino acid characters
      }
      if (targetModMap.containsKey(C)) {
        continue; // skip target residues
      }
      varModMap.put(c, mass);
    }

    for (String c : varModMap.keySet()) {
      String C = c.toUpperCase();
      double mass = AAmassMap.get(C) + varModMap.get(c);
      AAmassMap.put(c, mass);
    }
  }


  static void read_in_spectra()
      throws IOException, IllegalStateException, FileParsingException {

    System.err.println("\nReading spectra from " + Globals.spectrumPath.getCanonicalPath() + "  ("
        + Globals.spectrumSuffix.toUpperCase() + " format)");
    System.err.println("This can take a while so please be patient.");

    final List<String> rawFileNames = PSM_list.stream()
        .map(psm -> psm.srcFile).distinct().sorted()
        .collect(Collectors.toList());
    System.err.printf("Input PSMs originated from %d spectral data files:%n", rawFileNames.size());
    rawFileNames.forEach(file -> System.err.printf("\t%s%n", file));
    // test file paths
    final Set<String> rawPathNonExistent = rawFileNames.stream()
        .filter(file -> !Files.exists(Paths.get(Globals.spectrumPath.toString(), file)))
        .collect(Collectors.toSet());
    if (!rawPathNonExistent.isEmpty()) {
      System.err.printf("Files don't exist, PSMs from them will be skipped:%n");
      rawPathNonExistent.forEach(path -> System.err.printf("\t%s%n", path));
    }

    // map whatever paths are in PSMs list to absolute paths
    final Map<String, Path> mapPaths = rawFileNames.stream()
        .filter(file -> !rawPathNonExistent.contains(file))
        .collect(Collectors.toMap(f -> f, f -> Paths.get(Globals.spectrumPath.toString(), f).toAbsolutePath()));

    TreeMap<String, TreeMap<Integer, PSM>> mapFilesToPsms = new TreeMap<>();
    final ArrayList<PSM> psmsKept = new ArrayList<>(PSM_list.size());
    PSM_list.stream()
        .filter(psm -> !rawPathNonExistent.contains(psm.srcFile))
        .forEach(psm -> {
          psmsKept.add(psm);
          TreeMap<Integer, PSM> map = mapFilesToPsms.get(psm.srcFile);
          if (map == null) map = new TreeMap<>();
          map.put(psm.scanNum, psm);
          mapFilesToPsms.put(mapPaths.get(psm.srcFile).toString(), map);
        });
    PSM_list = psmsKept;

    if (Globals.spectrumSuffix.equalsIgnoreCase("mgf")) {

      for (String specFile : mapFilesToPsms.keySet()) {
        TIntObjectHashMap<SpectrumClass> curSpectra = ReaderHelper.read_mgf(specFile);
        String fileName = new File(specFile).getName(); // get just the file name of specFile
        int assignedSpectraCtr = 0;

        // Assign the spectra to their respective PSMs
        for (PSM p : PSM_list) {
          if (p.srcFile.equalsIgnoreCase(fileName)) {
            if (curSpectra.containsKey(p.scanNum)) {
              p.recordSpectra(curSpectra.get(p.scanNum));
              assignedSpectraCtr++;
            }
          }
        }
        System.err.println(fileName + ": " + assignedSpectraCtr + " spectra read in.");
      }
    } else if (Globals.spectrumSuffix.equalsIgnoreCase("mzXML")) {
      ReaderHelper.read_mzXML(mapFilesToPsms, numThreads);
    } else if (Globals.spectrumSuffix.equalsIgnoreCase("mzML")) {
      ReaderHelper.read_mzML(mapFilesToPsms, numThreads);
    }
  }


  static double getFragmentIonMass(String x, double z, double addl_mass) {
    double ret = 1.00728 * z;

    int start = x.indexOf(":") + 1;
    int stop = x.length();

    if (x.contains("-")) {
      stop = x.indexOf("-");
    }

    for (int i = start; i < stop; i++) {
      String c = Character.toString(x.charAt(i));

      if (AAmassMap.containsKey(c)) {
        double mass = AAmassMap.get(c);
        ret += mass;
      }
    }

    ret += addl_mass; // y-ions have a water molecule extra

    return ret;
  }


  // Function returns the decoy 'key' from the decoy map for the given residue
  static String getDecoySymbol(char c) {
    String srcChar = Character.toString(c);

    // TODO: ACHTUNG: Continue here, change the underlying map and this whole method.

    for (Character k : decoyAAMap.keySet()) {
      String v = decoyAAMap.get(k);

      if (v.equalsIgnoreCase(srcChar)) {
        return k.toString();
        //break;
      }
    }

    return "";
  }


  // Function returns the TPP-formatted representation of the given single
  // character modification
  static String getTPPresidue(String c) {
    String ret = "";
    String orig;

    if (c.equalsIgnoreCase("[")) {
      int d = (int) Math.round(Globals.ntermMass) + 1; // adds a proton
      ret = "n[" + String.valueOf(d) + "]";
    } else if (c.equals("]")) {
      int d = (int) Math.round(Globals.ctermMass);
      ret += "c[" + String.valueOf(d) + "]";
    } else {
      int i = (int) Math.round(AAmassMap.get(c));

      if (isDecoyResidue(c)) {
        orig = decoyAAMap.get(c);
      } else {
        orig = c.toUpperCase();
      }

      ret = orig + "[" + String.valueOf(i) + "]";
    }

    return ret;
  }


  // Function returns true if the given sequence contains decoy characters
  static boolean isDecoySeq(String seq) {

    for (int i = 0; i < seq.length(); i++) {
      char c = seq.charAt(i);
      if (Character.compare('[', c) == 0) {
        continue; // n-term mod symbol
      }
      if (Character.compare(']', c) == 0) {
        continue; // c-term mod symbol
      }

      if (isDecoyResidue(Character.toString(c))) {
        return true;
      }
    }

    return false;
  }


  // Function returns true if the given residue is from the decoy list
  static boolean isDecoyResidue(String AA) {
    return decoyAAMap.containsKey(AA);
  }


  // Record the global and local FLR values estimated for all of the delta scores
  public static void recordFLRestimates() {
    FLRestimateMap = new THashMap<>();

    for (PSM p : Globals.PSM_list) {
      if (p.isDecoy) {
        continue; // skip FLR data from decoys
      }
      double[] d = new double[2];
      d[0] = p.globalFDR;
      d[1] = p.localFDR;
      FLRestimateMap.put(p.deltaScore, d);
    }
  }


  // Function assigns the global and local FLR for the current PSM from the FLRestimateMap
  public static void assignFLR() {

    ArrayList<Double> obsDeltaScores = new ArrayList<>(FLRestimateMap.keySet());
    Collections.sort(obsDeltaScores);  // sort them from low to high
    int N = obsDeltaScores.size();
    boolean assigned;
    for (PSM p : Globals.PSM_list) {
      double obs_ds = p.deltaScore;
      assigned = false;

      // iterate over the delta scores until you find the value closest to this one
      for (int i = 1; i < N; i++) {
        double curDS = obsDeltaScores.get(i);
        if (curDS > obs_ds) { // hit the limit, get the *previous* delta score
          double[] d = FLRestimateMap.get(obsDeltaScores.get((i - 1)));
          p.globalFDR = d[0];
          p.localFDR = d[1];
          assigned = true;
          break;
        }
      }

      if (!assigned) { // very high scoring PSM
        double[] d = FLRestimateMap.get(obsDeltaScores.get((N - 1)));
        p.globalFDR = d[0];
        p.localFDR = d[1];
      }
    }
  }


  // This function prepares each PSM for the second iteration (the one after the FLR has been estimated)
  public static void clearPSMs() {
    for (PSM p : PSM_list) {
      p.clearScores();
    }
  }
}
