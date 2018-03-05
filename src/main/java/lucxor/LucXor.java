/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import edu.umich.dmtavt.ptmlocal.LucxorParams;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.hash.THashSet;
import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.pepxml.PepXmlParser;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.MsmsPipelineAnalysis;

/**
 * @author dfermin
 */
public class LucXor {
  public static final String releaseVersion = "1.2014Oct10";

  private static long startTime;
  private static long endTime;

  public static void main(String[] args)
      throws IOException, IllegalStateException, InterruptedException, ExecutionException, FileParsingException {

    // Get the release version of the program from the ANT build.xml file
    //String releaseVersion = LucXor.class.getPackage().getImplementationVersion();

//    System.out.println("Program started, press enter to continue...");
//    Scanner scanner = new Scanner(System.in);
//    scanner.nextLine();

    System.err.print("\nluciphor2 (JAVA-based version of Luciphor)\n" +
        "Version: " + releaseVersion + "\n" +
        "Original C++ version available at: http://luciphor.sf.net\n\n");

    if (args.length < 1) {
      System.err.print("USAGE: java -jar luciphor2.jar <input_file>\n\n");
      System.err.print("\tGenerate a luciphor2 input file with: java -jar luciphor2.jar -t\n");
      System.err
          .print("\tModify the input file to suit your needs and submit it to the program.\n");
      System.err.print("\tExample: java -jar luciphor2.jar input_file_you_edited\n\n");
      System.exit(0);
    }

    if (args[0].equalsIgnoreCase("-t")) {
      LucxorParams.writeTemplateInputFile();
    }

    // Start recording the running of the program
    startTime = System.nanoTime();

    Globals.initialize(); // Initialize global variables

    final LucxorParams params = LucxorParams.parseParameterFile(args[0]);

    Globals.loadUserMods(params);

    // read in the score results
    if (Globals.inputType == Constants.PEPXML) {
      parse_pepXML();
    } else {
      parse_TSV_src();
    }

    Globals.read_in_spectra(); // Read in spectra for these PSMs

    switch (Globals.scoringAlgorithm) {
      case Constants.CID:
        runCIDcode();
        break;
      case Constants.HCD:
        runHCDcode();
        break;
    }

    writeResults();

    endTime = System.nanoTime();
    reportElapsedTime();
  }


  // Function computes the elapsed run time of the program and reports it
  private static void reportElapsedTime() {
    long elapsedTime = endTime - startTime;

    long milliseconds = TimeUnit.NANOSECONDS.toMillis(elapsedTime);
    int ss = (int) (milliseconds / 1000) % 60;
    int mm = (int) ((milliseconds / (1000 * 60)) % 60);
    int hh = (int) ((milliseconds / (1000 * 60 * 60)) % 24);

    System.err.println("\nTotal run time (HH:MM:SS) = " +
        String.format("%02d", hh) + ":" +
        String.format("%02d", mm) + ":" +
        String.format("%02d", ss) + "\n");
  }


  // Function parses a pepXML file (this is the default input for luciphor)
  private static void parse_pepXML() throws FileParsingException {
    System.err.println("\nReading PSMs from pepXML file: " + Globals.inputFile.getAbsolutePath());


    final MsmsPipelineAnalysis pepxml = PepXmlParser.parse(Paths.get(Globals.inputFile.getAbsolutePath()));
    PepXML p = new PepXML();
    p.load(pepxml);
    //p.load(Globals.inputFile);

    System.err.println(Globals.PSM_list.size() + " Candidate PSMs read in.");
  }


  // Function parses a tab-delimited file of PSMs
  private static void parse_TSV_src() throws IOException {
    System.err.println("\nReading PSM from TSV file: " + Globals.inputFile.getAbsolutePath());
    PSM curPSM;

    if (!Globals.inputFile.exists()) {
      System.err.println("ERROR: Unable to find " + Globals.inputFile.getAbsolutePath());
      System.exit(0);
    }

    BufferedReader br = new BufferedReader(new FileReader(Globals.inputFile));
    String line;
    boolean passedHdr = false; // used only if the input file has a header line

    while ((line = br.readLine()) != null) {
      if (line.startsWith("#")) {
        continue;
      }

      // Skip over the header line for the TSV file (if it has one)
      if ((Globals.tsvHdr == 1) && (!passedHdr)) {
        passedHdr = true;
        continue;
      }

      String[] vars = line.split("\t");

      curPSM = new PSM();
      curPSM.srcFile = vars[0];
      curPSM.scanNum = Integer.valueOf(vars[1]);
      curPSM.charge = Integer.valueOf(vars[2]);
      curPSM.PSMscore = Double.valueOf(vars[3]);
      curPSM.origPep.peptide = vars[4].toUpperCase();

      if (vars.length < 6) { // no modifications so skip this PSM
        continue;
      }

      // modification string syntax: <pos>=<mass_of_modified_AA>
      // multiple modifications are comma separated (no spaces)
      // N-term pos = -100, C-term pos = 100
      // Fixed modifications are not reported
      for (String modSites : vars[5].split(",")) {
        String[] ary = modSites.split("=");
        int pos = Integer.valueOf(ary[0]); // assume 0-based coordinates
        double mass = Double.valueOf(ary[1]);

        curPSM.modCoordMap.put(pos, mass);
      }

      // Determine if this PSM has any non-standard AA characters.
      // If so, discard it
      int numBadChars = 0;
      for (int i = 0; i < vars[4].length(); i++) {
        String c = Character.toString(vars[4].charAt(i));
        if (!"ACDEFGHIKLMNPQRSTVWY".contains(c)) {
          numBadChars++;
        }
      }

      // Skip PSMs that exceed the number of candidate permutations
      if (curPSM.origPep.getNumPerm() > Globals.max_num_permutations) {
        numBadChars = 100;
      }

      if (numBadChars == 0) {
        curPSM.init();
        if (curPSM.isKeeper) {
          Globals.PSM_list.add(curPSM);
        }
      }
    }
    br.close();

    System.err.println("Read in " + Globals.PSM_list.size() + " PSMs");
  }


  // Function executes all the commands that are specific to the CID-mode
  // scoring algorithm
  private static void runCIDcode()
      throws IOException, InterruptedException, ExecutionException {

    System.err.println("\nRunning in CID mode.\n");

    int NCPU = Globals.numThreads;
    if (NCPU > 1) {
      if (NCPU < Runtime.getRuntime().availableProcessors()) {
        NCPU = Globals.numThreads + 1;
      }
    } else {
      NCPU = 1;
    }

    Globals.modelingMap_CID = new THashMap<>();
    int numPSM = 0; // track number of PSMs used for modeling

    // First make sure you have enough PSMs for each charge state to create
    // an accurate model.
    TIntIntHashMap chargeMap = new TIntIntHashMap();
    for (int z = 2; z <= Globals.maxChargeState; z++) {
      chargeMap.put(z, 0);
    }
    for (PSM p : Globals.PSM_list) {
      if (p.useForModel) {
        int old = 1;
        if (chargeMap.containsKey(p.charge)) {
          old = chargeMap.get(p.charge) + 1;
        }
        chargeMap.put(p.charge, old);
        numPSM++;
      }
    }

    // remove charge states you can't model
    TIntHashSet badZ = new TIntHashSet();
    System.err.print("PSMs for modeling:\n------------------\n");
    for (int z = 2; z <= Globals.maxChargeState; z++) {
      int n = chargeMap.get(z);
      System.err.println("+" + z + ": " + n + " PSMs");

      if (n < Globals.minNumPSMsForModeling) {
        badZ.add(z);
      }
    }
    System.err.print("\n");

    // Check to see if you have enough hi-scoring PSMs to construct an accurate model.
    // If not quit now
    if ((numPSM < Globals.minNumPSMsForModeling) || (badZ.size() == chargeMap.size())) {
      System.err.println("You do not have enough PSMs with a score > " + Globals.modelTH +
          " to accurately model the data. (Minimum number of PSMs required per charge state: "
          + Globals.minNumPSMsForModeling + ")\n" +
          "Exiting now.\n"
      );
      System.exit(0);
    }

    // iterate over the charge states collecting data for all modeling PSMs
    System.err.println("(CID) Building Parametric Models from high-scoring PSMs...");
    for (int z = 2; z <= Globals.maxChargeState; z++) {
      ArrayList<PeakClass> modelingPks = new ArrayList<>();

      if (badZ.contains(z)) {
        continue; // skip this charge state, you don't have enough data to model it
      }

      numPSM = 0; // track number of PSMs used for modeling

      RunModelParameterWorkers(z);

      // At this point, each PSM (whether you used threading or not) has all of it's
      // peaks classified as 'pos' or 'neg'. Record them into the modelingPks object
      numPSM = RecordClassifiedPeaks(numPSM, z, modelingPks);

      if (modelingPks.isEmpty()) {
        continue;
      }

      ModelData_CID M = new ModelData_CID(z, modelingPks);

      if (Globals.debugMode == Constants.WRITE_MODEL_PKS) {
        M.writeModelPks();
      }

      M.numPSM = numPSM;
      Globals.modelingMap_CID.put(z, M);
      modelingPks.clear();
    } // end loop over charge state for modeling

    if (Globals.modelingMap_CID.size() < 1) {
      System.err.println("\nInsufficient data to construct model.\nExiting now.\n");
      System.exit(0);
    }

    // Compute the statistics for the collected charge states
    TIntArrayList missedChargeStates = new TIntArrayList();
    int maxObsZ = 0;  // highest modeled charge state.
    for (int z = 2; z <= Globals.maxChargeState; z++) {
      if (!Globals.modelingMap_CID.containsKey(z)) {
        missedChargeStates.add(z);
        continue;
      } else {
        ModelData_CID m = Globals.modelingMap_CID.get(z);
        m.calcMean();
        m.calcVar();
        m.printSummaryStats();
        m.clearArrays();
        Globals.modelingMap_CID.put(z, m);
        if (z > maxObsZ) {
          maxObsZ = z;
        }
      }
    }

    // Back fill any missing charge states using the data for charge state 'maxObsZ'
    ModelData_CID m = Globals.modelingMap_CID.get(maxObsZ);
    for (int missedZ : missedChargeStates.toArray()) {
      Globals.modelingMap_CID.put(missedZ, m);
    }
    missedChargeStates.clear();

    System.gc(); // try and reclaim some memory

    //***********************************************************************
    // Score the PSMs
    //************************************************************************/

    for (int RN = 0; RN < 2; RN++) { // RN = run number, 0 = calc. FDR 1 = assign FDR estimates

      printRunInfo(NCPU, RN);

      // for single threaded analysis
      if ((Globals.numThreads == 1) || (Globals.debugMode != 0)) {
        scorePermutationsAndWriteScoredPeaks(RN);
      } else { // multi-threaded scoring

        // A pool of 'N' worker threads to do the work
        final int NTHREADS = Globals.numThreads;
        ExecutorService executor = Executors.newFixedThreadPool(NTHREADS);
        int ctr = 1;
        for (PSM p : Globals.PSM_list) {
          Runnable worker = new ScoringWorkerThread(p, RN, ctr++);
          executor.execute(worker);
        }
        executor.shutdown();

        // wait for all jobs to finish
        if( !executor.awaitTermination((long) Constants.FUNCTION_TIME_LIMIT, TimeUnit.SECONDS) ) {
          System.err.println("\nThread queue failed to terminate nicely\n");
          executor.shutdownNow();
        }

      }
      System.err.print("\n");

      // first iteration is done, compute the FLR from the PSMs
      if (computeFlrFromPsms(RN)) {
        break;
      }

      if (RN == 1) { // now assign FLRs to each PSM
        Globals.assignFLR();
      }
    } // end loop over RN
    //*******************************************************/
  }

  private static boolean computeFlrFromPsms(int RN)
      throws InterruptedException, ExecutionException {
    if (RN == 0) {
      StatsFunctions.calcFLR();

      if (Globals.runMode == Constants.DEFAULT_RUN_MODE) {
        Globals.recordFLRestimates();
        Globals.clearPSMs();
      } else { // you are done, exit loop over RN variable
        return true;
      }
    }
    return false;
  }

  private static void scorePermutationsAndWriteScoredPeaks(int RN) throws IOException {

    int ctr = 0;
    for (PSM p : Globals.PSM_list) {

      p.generatePermutations(RN);
      p.scorePermutations();

      if (Globals.debugMode == Constants.WRITE_SCORED_PKS) {
        p.debug_writeScoredPeaks();
      }


      ctr++;
      if ((ctr % 100) == 0) {
        System.err.print(ctr + " ");
      }
      if ((ctr % 1000) == 0) {
        System.err.print("\n");
      }
    }
  }

  private static void printRunInfo(int NCPU, int RN) {
    if (Globals.runMode == Constants.DEFAULT_RUN_MODE) {

      if (RN == 0) {
        System.err
            .println("\n[ " + RN + " ] Estimating FLR with decoys (" + NCPU
                + " threads)...");
      }
      if (RN == 1) {
        System.err.println(
            "\n[ " + RN + " ] Scoring " + Globals.PSM_list.size() + " PSMs (" + NCPU
                + " threads)...");
      }
    } else {
      System.err
          .println("\nScoring " + Globals.PSM_list.size() + " PSMs (" + NCPU + " threads)...");
    }
  }

  private static int RecordClassifiedPeaks(int numPSM, int z, ArrayList<PeakClass> modelingPks) {
    for (PSM p : Globals.PSM_list) {
      if ((p.charge == z) && p.useForModel) {
        if ((null != p.posPeaks) && (!p.posPeaks.isEmpty())) {
          modelingPks.addAll(p.posPeaks);
          p.posPeaks.clear();
        }
        if ((null != p.negPeaks) && (!p.negPeaks.isEmpty())) {
          modelingPks.addAll(p.negPeaks);
          p.negPeaks.clear();
        }
        numPSM++;
      }
    }
    return numPSM;
  }

  private static void RunModelParameterWorkers(int z) throws InterruptedException {
    if ((Globals.numThreads == 1) || (Globals.debugMode != 0)) {

      for (PSM p : Globals.PSM_list) {
        if ((p.charge == z) && p.useForModel) {
          p.generatePermutations(0); // generate both real permutations
          p.matchAllPeaks();
        }
      }
    } else { // multi-threaded modeling

      // A pool of 'N' worker threads
      final int NTHREADS = Globals.numThreads;
      ExecutorService executor = Executors.newFixedThreadPool(NTHREADS);
      int ctr = 1;
      for (PSM p : Globals.PSM_list) {
        if ((p.charge == z) && p.useForModel) {
          Runnable worker = new ModelParameterWorkerThread(p, ctr++);
          executor.execute(worker);
        }
      }
      executor.shutdown();

      // wait for all jobs to finish
      if( !executor.awaitTermination((long) Constants.FUNCTION_TIME_LIMIT, TimeUnit.SECONDS) ) {
        System.err.println("\nThread queue failed to terminate nicely\n");
        executor.shutdownNow();
      }

    } // end of multi-threaded approach
  }


  // Function executes all the commands that are specific to the HCD-mode
  // scoring algorithm
  private static void runHCDcode() throws IOException, InterruptedException, ExecutionException {
    System.err.println("\nRunning in HCD mode.\n");

    int numPSM = 0; // track number of PSMs used for modeling

    int NCPU = Globals.numThreads;
    if (NCPU < Runtime.getRuntime().availableProcessors()) {
      NCPU = Globals.numThreads + 1;
    }

    Globals.modelingMap_HCD = new THashMap<>();

    // First make sure you have enough PSMs for each charge state to create
    // an accurate model.
    THashMap<Integer, Integer> chargeMap = new THashMap<>();
    for (int z = 2; z <= Globals.maxChargeState; z++) {
      chargeMap.put(z, 0);
    }
    for (PSM p : Globals.PSM_list) {
      if (p.useForModel) {
        int old = 1;
        if (chargeMap.containsKey(p.charge)) {
          old = chargeMap.get(p.charge) + 1;
        }
        chargeMap.put(p.charge, old);
        numPSM++;
      }
    }

    // remove charge states you can't model
    THashSet<Integer> badZ = new THashSet<>();
    System.err.print("PSMs for modeling:\n------------------\n");
    for (int z = 2; z <= Globals.maxChargeState; z++) {
      int n = chargeMap.get(z);
      System.err.println("+" + z + ": " + n + " PSMs");

      if (n < Globals.minNumPSMsForModeling) {
        badZ.add(z);
      }
    }
    System.err.print("\n");

    // Check to see if you have enough hi-scoring PSMs to construct an accurate model.
    // If not quit now
    if ((numPSM < Globals.minNumPSMsForModeling) || (badZ.size() == chargeMap.size())) {
      System.err.println("You do not have enough PSMs with a score > " + Globals.modelTH +
          " to accurately model the data. (Minimum number of PSMs required per charge state: "
          + Globals.minNumPSMsForModeling + ")\n" +
          "Exiting now.\n"
      );
      System.exit(0);
    }

    // iterate over the charge states collecting data for all modeling PSMs
    System.err.println("(HCD) Acquiring Non-Parametric Model features from high-scoring PSMs...");

    for (int z = 2; z <= Globals.maxChargeState; z++) {
      ArrayList<PeakClass> modelingPks = new ArrayList<>();

      if (badZ.contains(z)) {
        continue; // skip this charge state you don't have enough data to model it
      }

      numPSM = 0; // track number of PSMs used for modeling

      RunModelParameterWorkers(z);

      // At this point, each PSM (whether you used threading or not) has all of it's
      // peaks classified as 'pos' or 'neg'. Record them into the modelingPks object
      numPSM = RecordClassifiedPeaks(numPSM, z, modelingPks);

      if (modelingPks.isEmpty()) {
        continue;
      }

      ModelData_HCD M = new ModelData_HCD(z, modelingPks);
      M.numPSM = numPSM;
      Globals.modelingMap_HCD.put(z, M);

      if (Globals.debugMode == Constants.WRITE_MODEL_PKS) {
        M.writeModelPks();
      }

    } // end for loop over charge state for modeling

    if (Globals.modelingMap_HCD.size() < 1) {
      System.err.println("\nInsufficient data to construct model.\nExiting now.\n");
      System.exit(0);
    }

    // Compute the statistics for the collected charge states
    ArrayList<Integer> missedChargeStates = new ArrayList<>();
    int maxObsZ = 0;  // highest modeled charge state.
    for (int z = 2; z <= Globals.maxChargeState; z++) {

      if (!Globals.modelingMap_HCD.containsKey(z)) {
        missedChargeStates.add(z);
        continue;
      } else {
        ModelData_HCD m = Globals.modelingMap_HCD.get(z);

        m.calcMean();
        m.calcVar();

        m.printStats();

        // Compute intensity parameters
        m.estimateNP_intensity('b');
        m.estimateNP_intensity('y');
        m.estimateNP_intensity('n');

        // Compute distance parameters
        m.estimateNP_posDist();
        System.err.print("\n");  // makes for prettier output

        Globals.modelingMap_HCD.put(z, m);
        if (z > maxObsZ) {
          maxObsZ = z;
        }
      }
    }

    // Back fill model data for missing charge states
    ModelData_HCD m = Globals.modelingMap_HCD.get(maxObsZ);
    for (int missedZ : missedChargeStates) {
      Globals.modelingMap_HCD.put(missedZ, m);
    }
    missedChargeStates.clear();

    // for debugging
    if (Globals.debugMode == Constants.WRITE_HCD_NONPARAM) {
      for (ModelData_HCD hcd : Globals.modelingMap_HCD.values()) {
        hcd.write_density_data(1); // intensity
        hcd.write_density_data(2); // distance
      }
    }

    // Clean up before you start scoring
    System.gc();

    //*******************************************************************
    // Score the PSMs
    //********************************************************************/
    for (int RN = 0; RN < 2; RN++) { // RN = run number

      printRunInfo(NCPU, RN);

      // for single threaded analysis
      if ((Globals.numThreads == 1) || (Globals.debugMode != 0)) {
        System.err.println("Scoring single-threaded");
        scorePermutationsAndWriteScoredPeaks(RN);

      } else { // multi-threaded scoring
        System.err.println("Scoring multi-threaded");
        final int NTHREADS = Globals.numThreads;
        ExecutorService executor = Executors.newFixedThreadPool(NTHREADS);
        int ctr = 1;
        for (PSM p : Globals.PSM_list) {
          p.generatePermutations(RN);
          Runnable worker = new ScoringWorkerThread(p, RN, ctr++);
          executor.execute(worker);
        }
        executor.shutdown();

        // wait for all jobs to finish
        if( !executor.awaitTermination((long) Constants.FUNCTION_TIME_LIMIT, TimeUnit.SECONDS) ) {
            System.err.println("\nThread queue failed to terminate nicely\n");
            executor.shutdownNow();
        }

      }
      System.err.print("\n");

      // first iteration is done, compute the FLR from the PSMs
      if (computeFlrFromPsms(RN)) {
        break;
      }

      if (RN == 1) { // now assign FLRs to each PSM
        Globals.assignFLR();
      }
    } // end loop over RN
    //*******************************************************/
  }


  // Function writes the collected results to disk
  private static void writeResults() throws IOException {

    // set the output file's name
    if (Globals.outputFile.equalsIgnoreCase("luciphor_results.tsv")) {
      Globals.outputFile = "luciphor_results." + Globals.timeStamp + ".tsv";
    }

    File fileOut = new File(Globals.outputFile);
    System.err.println("\nResults written to '" + fileOut.getAbsoluteFile() + "'");

    // You will only need these variables when the user wants the matched
    // peaks to be written to disk

    try (BufferedWriter bw = new BufferedWriter(new FileWriter(fileOut.getAbsoluteFile()))) {
      // Print column headers in output file
      String hdr = "specId\t"
          + "peptide\t"
          + "predictedPep1\t"
          + "predictedPep2\t"
          + "numPPS\t"
          + "numRPS\t";

      switch (Globals.scoringMethod) {
        case Constants.MASCOTIONSCORE:
          hdr += "MascotIonScore\t";
          break;
        case Constants.NEGLOGEXPECT:
          hdr += "negLogExpect\t";
          break;
        case Constants.PEPPROPHET:
          hdr += "pepProphet\t";
          break;
      }

      if (Globals.runMode == Constants.REPORT_DECOYS) {
        hdr += "isDecoy1\tisDecoy2\t";
      }

      hdr += "deltaScore\t"
          + "pep1score\tpep2score\t"
          + "globalFLR\tlocalFLR\n";

      bw.write(hdr);

      for (PSM psm : Globals.PSM_list) {
        bw.write(psm.getResults());
      }
    }

    if (Globals.writeMatchedPeaks) {
      String suffixDir = fileOut.getParent();
      Globals.matchedPkFile = suffixDir + "/luciphor_matchedPks." + Globals.timeStamp + ".tsv";
      System.err.print("Matched peaks will be written to '" + Globals.matchedPkFile + "'\n");
      File matchedPksF = new File(Globals.matchedPkFile);
      try (BufferedWriter bwPks = new BufferedWriter(
          new FileWriter(matchedPksF.getAbsoluteFile()))) {
        // Column header for the peak file
        String pkHdr = "specId\tpepNum\tpredictPep\tfragmentIon\t" +
            "m/z\trelIntensity\tDscore\tIscore\tscore\n";
        bwPks.write(pkHdr);
        for (PSM psm : Globals.PSM_list) {
          bwPks.write(psm.writeMatchedPks());
        }
      }
    }

  }

}
