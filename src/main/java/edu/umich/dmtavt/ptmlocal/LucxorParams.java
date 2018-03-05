package edu.umich.dmtavt.ptmlocal;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.TreeMap;
import lucxor.Constants;

/**
 * Extracted user-defined params from Globals.java.
 *
 * @author Dmitry Avtonomov
 */
public class LucxorParams {
  public static final String DEFAULT_PARAMS_FILE_NAME = "luciphor2_input_template.txt";

  private Path spectralFilesPath = null;
  private File spectrumPath = null;
  private String spectrumSuffix = null;
  private String matchedPkFile = null;
  private File inputFile = null;
  private String outputFile = null;
  private String timeStamp = null;
  private String dateStamp = null;
  private int ms2tol_units = Integer.MIN_VALUE;
  private int inputType = Integer.MIN_VALUE;
  private int debugMode = Integer.MIN_VALUE;
  private int reduceNL = Integer.MIN_VALUE;
  private int peptideRepresentation = Integer.MIN_VALUE;
  private int scoringMethod = Integer.MIN_VALUE;
  private int scoringAlgorithm = Integer.MIN_VALUE;
  private int maxChargeState = Integer.MIN_VALUE;
  private int minNumPSMsForModeling = Integer.MIN_VALUE;
  private int numThreads = Integer.MIN_VALUE;
  private int runMode = Integer.MIN_VALUE;
  private int tsvHdr = Integer.MIN_VALUE;
  private int maxPepLen = Integer.MIN_VALUE;
  private double precursorNLmass = Double.NaN;
  private double ms2tol = Double.NaN;
  private double modelTH = Double.NaN;
  private double scoreTH = Double.NaN;
  private double decoyMass = Double.NaN;
  private double minMZ = Double.NaN;
  private double max_num_permutations = Double.NaN;
  private boolean writeMatchedPeaks = false;

  /** Holds all neutral loss masses, k= list of amino acids v= NL mass. */
  private TreeMap<String, Double> nlMap = new TreeMap<>();
  private TreeMap<String, Double> nlMapDecoy = new TreeMap<>();
  /** Mods user wants to search for. */
  private TreeMap<String, Double> targetModMap = new TreeMap<>();
  /** Fixed mods observed in data. This will only be filled from user params if input is TSV,
   * if input is PEPXML, this stays empty.*/
  private TreeMap<String, Double> fixedModMap = new TreeMap<>();
  /** Variable mods observed in data. This will only be filled from user params if input is TSV,
   * if input is PEPXML, this stays empty.*/
  private TreeMap<String, Double> varModMap = new TreeMap<>();

  public static LucxorParams parseParameterFile(String str) throws IOException {

    final LucxorParams p = new LucxorParams();

    File fileIn = new File(str);
    if (!fileIn.exists()) {
      System.err.print("\nERROR! Unable to open " + str + "\n\n");
      System.exit(0);
    }

    p.debugMode = 0; // 0 means no debugging output
    p.runMode = 0; // calculate FLR and rescore PSMs without decoys (2 iterations)
    p.minNumPSMsForModeling = 50;
    p.maxPepLen = 40;
    p.reduceNL = 0;
    p.numThreads = Runtime.getRuntime().availableProcessors();
    p.tsvHdr = 0; // default is no header
    p.writeMatchedPeaks = false; // true means you will generate the matched peaks

    try (BufferedReader br = new BufferedReader(new FileReader(fileIn))) {
      String line;
      while ((line = br.readLine()) != null) {
        if (line.startsWith("#")) {
          continue;
        }
        if (line.length() < 2) {
          continue;
        }

        if (line.startsWith("SPECTRUM_PATH")) {
          String s = parseInputLine(line);
          p.spectrumPath = new File(s).getCanonicalFile();
        }

        if (line.startsWith("SPECTRUM_SUFFIX")) {
          String s = parseInputLine(line);
          p.spectrumSuffix = s.toLowerCase();
        }

        if (line.startsWith("INPUT_DATA")) {
          String s = parseInputLine(line);
          p.inputFile = new File(s);
        }

        if (line.startsWith("OUTPUT_FILE")) {
          p.outputFile = parseInputLine(line);
        }

        if (line.startsWith("INPUT_TYPE")) {
          String s = parseInputLine(line);
          p.inputType = Integer.valueOf(s);
        }

        if (line.startsWith("MAX_PEP_LEN")) {
          String s = parseInputLine(line);
          p.maxPepLen = Integer.valueOf(s);
        }

        if (line.startsWith("MAX_NUM_PERM")) {
          String s = parseInputLine(line);
          p.max_num_permutations = Double.valueOf(s);
        }

        if (line.startsWith("MIN_NUM_PSMS_MODEL")) {
          String s = parseInputLine(line);
          p.minNumPSMsForModeling = Integer.valueOf(s);
        }

        if (line.startsWith("MS2_TOL") && !line.contains("_UNITS")) {
          String s = parseInputLine(line);
          p.ms2tol = Double.valueOf(s);
        }

        if (line.startsWith("MS2_TOL_UNITS")) {
          String s = parseInputLine(line);
          p.ms2tol_units = Integer.valueOf(s);
        }

        if (line.startsWith("ALGORITHM")) {
          String s = parseInputLine(line);
          p.scoringAlgorithm = Integer.valueOf(s);
        }

        if (line.startsWith("TSV_HEADER")) {
          String s = parseInputLine(line);
          p.tsvHdr = Integer.valueOf(s);
        }

        if (line.startsWith("REDUCE_PRECURSOR_NL")) {
          String s = parseInputLine(line);
          p.reduceNL = Integer.valueOf(s);
        }

        if (line.startsWith("PRECURSOR_NL_MASS_DIFF")) {
          String s = parseInputLine(line);
          p.precursorNLmass = Double.valueOf(s);
        }

        if (line.startsWith("SELECTION_METHOD")) {
          String s = parseInputLine(line);
          p.scoringMethod = Integer.valueOf(s);
        }

        if (line.startsWith("MODELING_SCORE_THRESHOLD")) {
          String s = parseInputLine(line);
          p.modelTH = Double.valueOf(s);
        }

        if (line.startsWith("MAX_CHARGE_STATE")) {
          String s = parseInputLine(line);
          p.maxChargeState = Integer.valueOf(s);
        }

        if (line.startsWith("NUM_THREADS")) {
          String s = parseInputLine(line);
          int x = Integer.valueOf(s);
          // We do minus 1 because 1 thread already goes to
          // running the whole program
          if (x < 0) {
            p.numThreads = 1;
          } else if (x > 1) {
            p.numThreads = (x - 1);
          } else if (x == 0) {
            p.numThreads = Runtime.getRuntime().availableProcessors();
          } else {
            p.numThreads = x;
          }
        }

        if (line.startsWith("DEBUG_MODE")) {
          String s = parseInputLine(line);
          p.debugMode = Integer.valueOf(s);
        }

        if (line.startsWith("WRITE_MATCHED_PEAKS_FILE")) {
          String s = parseInputLine(line);
          if (s.equals("1")) {
            p.writeMatchedPeaks = true;
          }
        }

        if (line.startsWith("RUN_MODE")) {
          String s = parseInputLine(line);
          p.runMode = Integer.valueOf(s);
          if (p.runMode > 1) {
            p.runMode = 0;
          }
        }

        if (line.startsWith("SCORING_THRESHOLD")) {
          String s = parseInputLine(line);
          p.scoreTH = Double.valueOf(s);
        }

        if (line.startsWith("DECOY_MASS")) {
          String s = parseInputLine(line);
          p.decoyMass = Double.valueOf(s);
        }

        if (line.startsWith("DECOY_NL")) {
          String[] ary = parseNeutralLossLine(line);
          String k = ary[0].substring(1);
          double m = Double.valueOf(ary[1]);
          p.nlMapDecoy.put(k, m);
        }

        if (line.startsWith("NL")) {
          String[] ary = parseNeutralLossLine(line);
          double m = Double.valueOf(ary[1]);
          p.nlMap.put(ary[0], m);
        }

        if (line.startsWith("MIN_MZ")) {
          String s = parseInputLine(line);
          p.minMZ = Double.valueOf(s);
        }

        if (line.startsWith("MOD_PEP_REP")) {
          String s = parseInputLine(line);
          p.peptideRepresentation = Integer.valueOf(s);
        }

        if (line.startsWith("TARGET_MOD")) {
          String[] ary = parseInputModLine(line);
          double m = Double.valueOf(ary[1]);
          p.targetModMap.put(ary[0].toUpperCase(), m);
        }

        // You only need to extract the VAR_MOD and FIXED_MOD values
        // from the input file if you are using TSV files
        if (p.inputType == Constants.TSV) {
          if (line.startsWith("VAR_MOD")) {
            String[] ary = parseInputModLine(line);
            double m = Double.valueOf(ary[1]);
            p.varModMap.put(ary[0].toLowerCase(), m); // variable mods are lower case
          }

          if (line.startsWith("FIXED_MOD")) {
            String[] ary = parseInputModLine(line);
            double m = Double.valueOf(ary[1]);
            p.fixedModMap.put(ary[0].toUpperCase(), m);
          }
        }
      }
      br.close();
    }

    if ((null == p.outputFile) || (p.outputFile.isEmpty())) {
      p.outputFile = "luciphor_results.tsv";
    }

    String classStr = "";
    switch (p.scoringMethod) {
      case Constants.PEPPROPHET:
        classStr = "Peptide Prophet Prob.";
        break;
      case Constants.MASCOTIONSCORE:
        classStr = "Mascot Ion Score";
        break;
      case Constants.NEGLOGEXPECT:
        classStr = "-log(Expect Value) (X!Tandem or Comet)";
        break;
      case Constants.XTDHYPERSCORE:
        classStr = "X!Tandem Hyperscore";
        break;
      case Constants.XCORR:
        classStr = "Sequest XCorr";
        break;
      default:
        System.err.print("\nERROR! Unknown scoring method: " + p.scoringMethod + "\n\n");
        System.exit(0);
        break;
    }

    int NCPU = p.numThreads;
    if ((p.numThreads > 1) && (p.numThreads < Runtime.getRuntime().availableProcessors())) {
      NCPU = (p.numThreads + 1);
    }

    System.err.println("Spectrum Path:           " + p.spectrumPath.getAbsolutePath());
    System.err.println("Spectrum Suffix:         " + p.spectrumSuffix);
    System.err.println("Input file:              " + p.inputFile);
    System.err
        .println("Input type:              " + (p.inputType == Constants.PEPXML ? "pepXML" : "tsv"));
    System.err.println(
        "MS2 tolerance:           " + p.ms2tol + (p.ms2tol_units == Constants.DALTONS ? " Da"
            : " ppm"));
    System.err
        .println("Luciphor Algorithm:      " + (p.scoringAlgorithm == Constants.CID ? "CID" : "HCD"));
    System.err.println("Classifying on:          " + classStr);
    System.err.println("Run Mode:                " + (p.runMode == 0 ? "Default" : "Report Decoys"));
    System.err.println("Num of Threads:          " + NCPU);
    System.err.println("Modeling Threshold:      " + p.modelTH);
    System.err.println("Scoring Threshold:       " + p.scoreTH);
    System.err.println("Permutation Limit:       " + p.max_num_permutations);
    System.err.println("Max peptide length:      " + p.maxPepLen);
    System.err.println("Min num PSMs for model:  " + p.minNumPSMsForModeling);
    System.err.println("Decoy Mass Adduct:       " + p.decoyMass);
    System.err.println("Max Charge State:        " + p.maxChargeState);
    System.err.println("Reduce NL:               " + (p.reduceNL == 0 ? "no" : "yes"));
    System.err.println("Output File:             " + p.outputFile);
    System.err.println("Write matched Peaks:     " + (p.writeMatchedPeaks ? "yes" : "no"));
    System.err.print("\n");

    if (p.debugMode != 0) {
      System.err.println("Debug mode:              " + p.debugMode + "  (Limiting to 1 CPU)\n");
      p.numThreads = 1;
    }

    System.err.println("Mods to score:");
    for (String s : p.targetModMap.keySet()) {
      System.err.println(s + "\t" + p.targetModMap.get(s));
    }

    if (!p.nlMap.isEmpty()) {
      System.err.println("\nAllowed Neutral Losses:");
      for (String s : p.nlMap.keySet()) {
        System.err.println(s + "\t" + p.nlMap.get(s));
      }
      for (String s : p.nlMapDecoy.keySet()) {
        System.err.println("<X>" + s + "\t" + p.nlMapDecoy.get(s) + "  (Decoy NL)");
      }
    }

    return p;
  }

  private static String parseInputLine(String line) {
    StringBuilder sb = new StringBuilder();
    int N = line.length();

    int b = line.indexOf("=") + 1;

    for (int i = b; i < N; i++) {
      char c = line.charAt(i);
      if (c == '#') {
        break;
      }
      if (c == ' ') {
        continue;
      }
      sb.append(c);
    }

    return sb.toString();
  }

  private static String[] parseInputModLine(String line) {
    String[] ret = new String[2];

    char aa = 0;
    StringBuilder sb = new StringBuilder();
    double mass;
    int N = line.length();

    int b = line.indexOf("=") + 1;

    for (int i = b; i < N; i++) {
      char c = line.charAt(i);

      if (c == '#') {
        break;
      }
      if (c == ' ') {
        continue;
      }

      if (Character.isAlphabetic(c)) {
        aa = c;
      } else if ((c == '[') || (c == ']')) {
        aa = c;
      } else {
        sb.append(c);
      }
    }

    mass = Double.valueOf(sb.toString());

    ret[0] = Character.toString(aa);
    ret[1] = String.valueOf(mass);

    return ret;
  }

  private static String[] parseNeutralLossLine(String line) {
    String[] ret = new String[2];

    line = line.replaceAll("#", "");

    String[] tmp = line.split("\\s+");

    ret[0] = tmp[2] + tmp[3];
    ret[1] = tmp[4];

    return ret;
  }


  /**
   * Write a sample input file for LucXor to disk.
   * @throws IOException
   */
  public static Path writeTemplateInputFile() throws IOException {

    final Path path = Paths.get(DEFAULT_PARAMS_FILE_NAME).toAbsolutePath();
    if (java.nio.file.Files.exists(path))
      throw new IllegalStateException(String.format("File already exists, won't overwrite: %s", path));

    try (BufferedWriter bw = new BufferedWriter(new FileWriter(path.toFile()))) {

      bw.write(
          "## Input file for Luciphor2 (aka: LucXor).\n## Anything after a hash '#' is ignored\n");
      bw.write(
          "## By default, these initial parameters are for performing a phosphorylation search\n\n");
      bw.write("SPECTRUM_PATH = <fill_in> ## specify the path to the spectra here\n");
      bw.write("SPECTRUM_SUFFIX = MGF     ## available options are MGF, mzML, or mzXML\n\n");

      bw.write("## Specify your input PSM results format\n" +
          "INPUT_DATA = <fill_in> ## specify the path to the pepXML or tab-delimited file here\n" +
          "INPUT_TYPE = 0   ## 0 = TPP pepXML, 1 = tab-delimited file\n");
      bw.write("ALGORITHM = 0    ## 0 = CID method, 1 = HCD method\n\n");

      bw.write("TSV_HEADER = 0 ## This pertains ONLY to tab-delimited (TSV) input data\n" +
          "               ## 0 = the input file does NOT contain column names as the first row\n" +
          "               ## 1 = the first row of the input file is the column names\n\n");

      bw.write("MS2_TOL = 0.5      ## MS/MS fragment ion tolerance\n" +
          "MS2_TOL_UNITS = 0  ## 0 = Daltons, 1 = PPM\n\n");

      bw.write("MIN_MZ = 150.0 ## do not consider peaks below this value " +
          "for matching fragment ions\n\n");

      bw.write("OUTPUT_FILE =  ## Specify the path to your desired output filename here\n" +
          "               ## A default value will be used if nothing is specified here.\n\n");

      bw.write(
          "WRITE_MATCHED_PEAKS_FILE = 0 ## Generate a tab-delimited file of all the matched peaks\n"
              +
              "                             ## for the top 2 predictions of each spectra\n" +
              "                             ## Useful for plotting spectra\n" +
              "                             ## 0 = no, 1 = yes\n\n");

      bw.write("## Place here any FIXED modifications that were used in your search\n" +
          "## This field is ONLY used for tab-delimited input\n" +
          "## Syntax: FIXED_MOD = <RESIDUE> <MODIFICATION_MASS>\n" +
          "## For N-terminal modifications use '[' for the residue character\n" +
          "## For C-terminal modifications use ']' for the residue character\n" +
          "FIXED_MOD = C 57.021464\n\n");

      bw.write(
          "## Place here any VARIABLE modifications that might be encountered that you don't\n" +
              "## want luciphor to score.\n" +
              "## This field is ONLY used for tab-delimited input\n" +
              "## Syntax: VAR_MOD = <RESIDUE> <MODIFICATION_MASS>\n" +
              "## For N-terminal modifications use '[' for the residue character\n" +
              "## For C-terminal modifications use ']' for the residue character\n" +
              "VAR_MOD = M 15.994915\n\n"
      );

      bw.write("## List the amino acids to be searched for and their mass modifications\n" +
          "## Syntax: TARGET_MOD = <RESIDUE> <MODIFICATION_MASS>\n" +
          "TARGET_MOD = S 79.966331\n" +
          "TARGET_MOD = T 79.966331\n" +
          "TARGET_MOD = Y 79.966331\n\n");

      bw.write("## List the types of neutral losses that you want to consider\n" +
          "## The residue field is case sensitive. For example: lower case 'sty' implies\n" +
          "## that the neutral loss can only occur if the specified modification is present\n" +
          "## Syntax: NL = <RESDIUES> -<NEUTRAL_LOSS_MOLECULAR_FORMULA> <MASS_LOST>\n" +
          "#NL = STE -H2O -18.01056    ## a correctly formatted example, (not actually recommended for phospho-searches)\n"
          +
          "#NL = RKQN -NH3 -17.026548  ## another correctly formatted example, (again not recommended for phospho-searches)\n"
          +
          "NL = sty -H3PO4 -97.97690\n\n");

      bw.write(
          "DECOY_MASS = 79.966331  ## how much to add to an amino acid to make it a decoy\n\n");
      bw.write("## For handling the neutral loss from a decoy sequence.\n" +
          "## The syntax for this is identical to that of the normal neutral losses given\n" +
          "## above except that the residue is always 'X'\n" +
          "## Syntax: DECOY_NL = X -<NEUTRAL_LOSS_MOLECULAR_FORMULA> <MASS_LOST>\n" +
          "DECOY_NL = X -H3PO4 -97.97690\n\n");

      bw.write(
          "MAX_CHARGE_STATE = 5 ## do not consider PSMs with a charge state above this value\n");
      bw.write(
          "MAX_PEP_LEN = 40 ## restrict scoring to peptides with a length shorter than this value\n");
      bw.write(
          "MAX_NUM_PERM = 16384 ## the maximum number of permutations a sequence can have\n\n");

      bw.write("SELECTION_METHOD = 0   ## 0 = Peptide Prophet probability (default)\n" +
          "                       ## 1 = Mascot Ion Score\n" +
          "                       ## 2 = -log(E-value)\n" +
          "                       ## 3 = X!Tandem Hyperscore\n" +
          "                       ## 4 = Sequest Xcorr\n\n");

      bw.write(
          "MODELING_SCORE_THRESHOLD = 0.95 ## minimum score a PSM needs to be considered for modeling\n");
      bw.write("SCORING_THRESHOLD = 0    ## PSMs below this value will be discarded\n");
      bw.write(
          "MIN_NUM_PSMS_MODEL = 50  ## The minimum number of PSMs you need for any charge state in order to build a model for it\n\n");

      bw.write(
          "MOD_PEP_REP = 0 ## 0 = show single character modifications, 1 = show TPP-formatted modifications\n\n");

      bw.write("NUM_THREADS = 0 ## For multi-threading, zero = use all CPU found by JAVA\n\n");

      bw.write("RUN_MODE = 0 ## Determines how Luciphor will run.\n" +
          "             ## 0 = Default: calculate FLR then rerun scoring without decoys (two iterations)\n"
          +
          "             ## 1 = Report Decoys: calculate FLR but don't rescore PSMs, all decoy hits will be reported\n\n");

      bw.write(
          "## This option can be used to help diagnose problems with Luciphor. Multi-threading is disabled in debug mode.\n"
              +
              "DEBUG_MODE = 0 ## 0 = default: turn off debugging\n" +
              "               ## 1 = write peaks selected for modeling to disk\n" +
              "               ## 2 = write the scores of all permutations for each PSM to disk\n" +
              "               ## 3 = write the matched peaks for the top scoring permutation to disk\n"
              +
              "               ## 4 = write HCD non-parametric models to disk (HCD-mode only option)\n\n");

    }

    return path;
  }

  public Path getSpectralFilesPath() {
    return spectralFilesPath;
  }

  public File getSpectrumPath() {
    return spectrumPath;
  }

  public String getSpectrumSuffix() {
    return spectrumSuffix;
  }

  public String getMatchedPkFile() {
    return matchedPkFile;
  }

  public File getInputFile() {
    return inputFile;
  }

  public String getOutputFile() {
    return outputFile;
  }

  public String getTimeStamp() {
    return timeStamp;
  }

  public String getDateStamp() {
    return dateStamp;
  }

  public int getMs2tol_units() {
    return ms2tol_units;
  }

  public int getInputType() {
    return inputType;
  }

  public int getDebugMode() {
    return debugMode;
  }

  public int getReduceNL() {
    return reduceNL;
  }

  public int getPeptideRepresentation() {
    return peptideRepresentation;
  }

  public int getScoringMethod() {
    return scoringMethod;
  }

  public int getScoringAlgorithm() {
    return scoringAlgorithm;
  }

  public int getMaxChargeState() {
    return maxChargeState;
  }

  public int getMinNumPSMsForModeling() {
    return minNumPSMsForModeling;
  }

  public int getNumThreads() {
    return numThreads;
  }

  public int getRunMode() {
    return runMode;
  }

  public int getTsvHdr() {
    return tsvHdr;
  }

  public int getMaxPepLen() {
    return maxPepLen;
  }

  public double getPrecursorNLmass() {
    return precursorNLmass;
  }

  public double getMs2tol() {
    return ms2tol;
  }

  public double getModelTH() {
    return modelTH;
  }

  public double getScoreTH() {
    return scoreTH;
  }

  public double getDecoyMass() {
    return decoyMass;
  }

  public double getMinMZ() {
    return minMZ;
  }

  public double getMax_num_permutations() {
    return max_num_permutations;
  }

  public boolean isWriteMatchedPeaks() {
    return writeMatchedPeaks;
  }

  public TreeMap<String, Double> getNlMap() {
    return nlMap;
  }

  public TreeMap<String, Double> getNlMapDecoy() {
    return nlMapDecoy;
  }

  public TreeMap<String, Double> getTargetModMap() {
    return targetModMap;
  }

  public TreeMap<String, Double> getFixedModMap() {
    return fixedModMap;
  }

  public TreeMap<String, Double> getVarModMap() {
    return varModMap;
  }
}
