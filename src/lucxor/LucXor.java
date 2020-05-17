/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;

import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.hash.THashSet;
import gnu.trove.set.hash.TIntHashSet;
import lombok.extern.slf4j.Slf4j;
import org.xml.sax.SAXException;
import umich.ms.fileio.exceptions.FileParsingException;

/**
 *
 * @author dfermin
 * @author ypriverol
 */
@Slf4j
public class LucXor {

    private static long startTime;
    private static long endTime;

    public static void main(String[] args) throws ParserConfigurationException,
            SAXException, IOException, IllegalStateException, DataFormatException,
            InterruptedException, ExecutionException, FileParsingException {
		
		String releaseVersion = "1.2014Oct10";
		log.info("\nluciphor2 (JAVA-based version of Luciphor)\n" +
						 "Version: " + releaseVersion + "\n" +
						 "Original C++ version available at: http://luciphor.sf.net\n\n");
				
		
		if(args.length < 1) {
			log.info("USAGE: java -jar luciphor2.jar <input_file>\n\n");
			log.info("\tGenerate a luciphor2 input file with: java -jar luciphor2.jar -t\n");
            log.info("\tModify the input file to suit your needs and submit it to the program.\n");
            log.info("\tExample: java -jar luciphor2.jar input_file_you_edited\n\n");
			System.exit(0);
		}
		
		if(args[0].equalsIgnoreCase("-t")) {
		    File outF = Globals.writeTemplateInputFile();
            log.info("\nPlease edit the input file: " + outF.getPath() +
                    " with your favorite text editor\n\n");
            System.exit(0);
        }

        // Start recording the running of the program
        startTime = System.nanoTime();

		// Read Configuration file
        Globals.parseInputFile(args[0]);

        // Load User Modifications
		Globals.loadUserMods();
		
		
		// read in the score results
		if(Globals.inputType == Constants.PEPXML)
			parsePepXML();
		else
			parseTSVSrc();

        // Read in spectra for these PSMs
		Globals.readInSpectra();

		// Run CID Algorithm
        if(Globals.scoringAlgorithm == Constants.CID)
            runCIDcode();

        // Run HCD Algorithm
    	if(Globals.scoringAlgorithm == Constants.HCD)
    	    runHCDcode();

    	//Write final results
		writeResults();

        endTime = System.nanoTime();
        reportElapsedTime();
	}


    /**
     * Function computes the elapsed run time of the program and reports it
     */
    private static void reportElapsedTime() {
        long elapsedTime = endTime - startTime;

        long milliseconds = TimeUnit.NANOSECONDS.toMillis(elapsedTime);
        int ss = (int) (milliseconds / 1000) % 60 ;
        int mm = (int) ((milliseconds / (1000*60)) % 60);
        int hh   = (int) ((milliseconds / (1000*60*60)) % 24);

         log.info("\nTotal run time (HH:MM:SS) = " +
                String.format("%02d", hh) + ":" +
                String.format("%02d", mm) + ":" +
                String.format("%02d", ss) + "\n");
    }


    /**
     * Function parses a pepXML file (this is the default input for luciphor)
     * @throws ParserConfigurationException
     * @throws SAXException
     * @throws IOException
     */
	private static void parsePepXML() throws ParserConfigurationException,
            SAXException, IOException {

	    log.info("\nReading PSMs from pepXML file: " +
                Globals.inputFile.getAbsolutePath());
		PepXML p = new PepXML(Globals.inputFile);
		log.info(Globals.psmList.size() + " Candidate PSMs read in.");
	}

	
	// Function parses a tab-delimited file of PSMs
	private static void parseTSVSrc() throws IOException {

	    log.info("\nReading PSM from TSV file: " + Globals.inputFile.getAbsolutePath());
		PSM curPSM = null;

		if( !Globals.inputFile.exists() ) {
			log.info("ERROR: Unable to find " + Globals.inputFile.getAbsolutePath());
			System.exit(0);
		}
		
		BufferedReader br = new BufferedReader(new FileReader(Globals.inputFile));
		String line;
        boolean passedHdr = false; // used only if the input file has a header line

		while( (line = br.readLine()) != null) {
			if(line.startsWith("#")) continue;

            // Skip over the header line for the TSV file (if it has one)
            if( (Globals.tsvHdr == 1) && (!passedHdr) ) {
                passedHdr = true;
                continue;
            }

            String[] vars = line.split("\t");
			
			curPSM = new PSM();
			curPSM.setSrcFile(vars[0]);
			curPSM.setScanNum(Integer.valueOf( vars[1] ));
			curPSM.setCharge(Integer.valueOf( vars[2] ));
			double obsScore = Double.valueOf( vars[3] );
			curPSM.setPeptideSequence(vars[4].toUpperCase());

			if(Globals.scoringMethod == Constants.NEGLOGEXPECT) {
			    curPSM.setPSMscore(-1.0 * Math.log(obsScore));
            } else {
			    curPSM.setPSMscore(obsScore);
            }


			if(vars.length < 6) { // no modifications so skip this PSM
				curPSM = null;
				continue;
			}
			
			// modification string syntax: <pos>=<mass_of_modified_AA>
			// multiple modifications are comma separated (no spaces)
			// N-term pos = -100, C-term pos = 100
			// Fixed modifications are not reported
			for(String modSites : vars[5].split(",")) {
				String[] ary = modSites.split("=");
				int pos = Integer.valueOf(ary[0]); // assume 0-based coordinates
				double mass = Double.valueOf(ary[1]);
				
				curPSM.getModCoordMap().put(pos, mass);
			}
			
			// Determine if this PSM has any non-standard AA characters. 
			// If so, discard it
			int numBadChars = 0;
			for(int i = 0; i < vars[4].length(); i++) {
				String c = Character.toString( vars[4].charAt(i) );
				if( !"ACDEFGHIKLMNPQRSTVWY".contains(c) ) numBadChars++;
			}

			// Skip PSMs that exceed the number of candidate permutations
			if(curPSM.getOrigPep().getNumPerm() > Globals.max_num_permutations) numBadChars = 100;

			if(numBadChars == 0) {
				curPSM.process();
				if(curPSM.isKeeper()) Globals.psmList.add(curPSM);
			}
			curPSM = null;
		}
		br.close();

		log.info("Read in " + Globals.psmList.size() + " PSMs");
	}


    /**
     * Function executes all the commands that are specific to the CID-mode
     * scoring algorithm
     * @throws IOException
     * @throws InterruptedException
     * @throws ExecutionException
     */
	private static void runCIDcode() throws IOException,
            InterruptedException, ExecutionException {
		
		log.info("\nRunning in CID mode.\n");

        int NCPU = Globals.numThreads;
        if(NCPU > 1) {
            if (NCPU < Runtime.getRuntime().availableProcessors()) NCPU = Globals.numThreads + 1;
        }
        else NCPU = 1;


		Globals.modelingMapCID = new THashMap<>();
		int numPSM = 0; // track number of PSMs used for modeling

		// First make sure you have enough PSMs for each charge state to create
		// an accurate model.
        TIntIntHashMap chargeMap = new TIntIntHashMap();
		for(int z = 2; z <= Globals.maxChargeState; z++) chargeMap.put(z,0);
		for(PSM p : Globals.psmList) {
			if(p.isUseForModel()) {
                int old = 1;
                if(chargeMap.containsKey(p.getCharge())) old = chargeMap.get(p.getCharge()) + 1;
				chargeMap.put(p.getCharge(), old);
                numPSM++;
			}
        }
		
		// remove charge states you can't model
        TIntHashSet badZ = new TIntHashSet();
        log.info("PSMs for modeling:\n------------------\n");
		for(int z = 2; z <= Globals.maxChargeState; z++) {
			int n = chargeMap.get(z);
            log.info("+" + z + ": " + n + " PSMs");

            if(n < Globals.minNumPSMsForModeling) badZ.add(z);
		}
        log.info("\n");

        // Check to see if you have enough hi-scoring PSMs to construct an accurate model.
        // If not quit now
        if( (numPSM < Globals.minNumPSMsForModeling) || (badZ.size() == chargeMap.size()) ) {
            log.info("You do not have enough PSMs with a score > " + Globals.modelTH +
                " to accurately model the data. (Minimum number of PSMs required per charge state: " + Globals.minNumPSMsForModeling + ")\n" +
                "Exiting now.\n"
            );
            System.exit(0);
        }

		
		// iterate over the charge states collecting data for all modeling PSMs
		log.info("(CID) Building Parametric Models from high-scoring PSMs...");
		for(int z = 2; z <= Globals.maxChargeState; z++) {
			ArrayList<Peak> modelingPks = new ArrayList<>();
			
			if(badZ.contains(z)) continue; // skip this charge state, you don't have enough data to model it
			
			numPSM = 0; // track number of PSMs used for modeling

            if((Globals.numThreads == 1) || (Globals.debugMode != 0)) {

                for (PSM p : Globals.psmList) {
                    if ((p.getCharge() == z) && p.isUseForModel()) {
                        p.generatePermutations(0); // generate both real permutations
                        p.matchAllPeaks();
                    }
                }
            }
            else { // multi-threaded modeling

                // A pool of 'N' worker threads
                final int NTHREADS = Globals.numThreads;
                ExecutorService executor = Executors.newFixedThreadPool(NTHREADS);
                int ctr = 1;
                for (PSM p : Globals.psmList) {
                    if ((p.getCharge() == z) && p.isUseForModel()) {
                        Runnable worker = new ModelParameterWorkerThread(p, ctr++);
                        executor.execute(worker);
                    }
                }
                executor.shutdown();

                // wait for all jobs to finish
                executor.awaitTermination((long) Constants.FUNCTION_TIME_LIMIT, TimeUnit.SECONDS);

                while (!executor.isShutdown()) {}
            } // end of multi-threaded approach


            // At this point, each PSM (whether you used threading or not) has all of it's
            // peaks classified as 'pos' or 'neg'. Record them into the modelingPks object
            for(PSM p : Globals.psmList) {
                if ((p.getCharge() == z) && p.isUseForModel()) {
                    if((null != p.getPosPeaks()) && (!p.getPosPeaks().isEmpty())) {
                        modelingPks.addAll(p.getPosPeaks());
                        p.getPosPeaks().clear();
                    }
                    if((null != p.getNegPeaks()) && (!p.getNegPeaks().isEmpty())) {
                        modelingPks.addAll(p.getNegPeaks());
                        p.getNegPeaks().clear();
                    }
                    numPSM++;
                }
            }
			
			if(modelingPks.isEmpty()) continue; 
			
			ModelDataCID M = new ModelDataCID(z, modelingPks);

            if(Globals.debugMode == Constants.WRITE_MODEL_PKS) M.writeModelPks();

            M.numPSM = numPSM;
			Globals.modelingMapCID.put(z, M);
			M = null;
			modelingPks.clear();
            modelingPks = null;
		} // end loop over charge state for modeling


        if(Globals.modelingMapCID.size() < 1) {
            log.info("\nInsufficient data to construct model.\nExiting now.\n");
            System.exit(0);
        }

		// Compute the statistics for the collected charge states
        TIntArrayList missedChargeStates = new TIntArrayList();
        int maxObsZ = 0;  // highest modeled charge state.
		for(int z = 2; z <= Globals.maxChargeState; z++) {
			if(!Globals.modelingMapCID.containsKey(z)) {
                missedChargeStates.add(z);
            }
			else {
                ModelDataCID m = Globals.modelingMapCID.get(z);
                m.calcMean();
                m.calcVar();
                m.printSummaryStats();
                m.clearArrays();
                Globals.modelingMapCID.put(z, m);
                m = null;
                if(z > maxObsZ) maxObsZ = z;
            }
		}

        // Back fill any missing charge states using the data for charge state 'maxObsZ'
        ModelDataCID m = Globals.modelingMapCID.get(maxObsZ);
        for(int missedZ : missedChargeStates.toArray()) {
            Globals.modelingMapCID.put(missedZ, m);
        }
		m = null;
		missedChargeStates.clear();

        System.gc(); // try and reclaim some memory



		// Score the PSMs
        for(int RN = 0; RN < 2; RN++) { // RN = run number, 0 = calc. FDR 1 = assign FDR estimates

            if(Globals.runMode == Constants.DEFAULT_RUN_MODE) {
                if(RN == 0)
                    log.info("\n[ " + RN + " ] Estimating FLR with decoys (" + NCPU + " threads)...");
                if(RN == 1)
                    log.info("\n[ " + RN + " ] Scoring " + Globals.psmList.size() + " PSMs (" + NCPU + " threads)...");
            }
            else {
                log.info("\nScoring " + Globals.psmList.size() + " PSMs (" + NCPU + " threads)...");
            }

            AtomicInteger ctr = new AtomicInteger(1);

            int currentRN = RN;
            Globals.psmList.parallelStream().forEach(p -> {
                p.generatePermutations(currentRN);
                try {
                    p.scorePermutations();
                    if(Globals.debugMode == Constants.WRITE_SCORED_PKS) p.debugWriteScoredPeaks();
                    ctr.getAndIncrement();
                    if( (ctr.get() % 100) == 0 ) log.info(ctr + " ");
                    if( (ctr.get() % 1000) == 0 ) log.info("\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
            });

            log.info("\n");

            // first iteration is done, compute the FLR from the PSMs
            if(RN == 0) {
               Globals.calcFLR();
                if(Globals.runMode == Constants.DEFAULT_RUN_MODE) {
                    Globals.recordFLRestimates();
                    Globals.clearPSMs();
                }
                else { // you are done, exit loop over RN variable
                    RN = 2;
                    break;
                }
            }
            if(RN == 1) { // now assign FLRs to each PSM
                Globals.assignFLR();
            }
        }
	}

    /**
     * Function executes all the commands that are specific to the HCD-mode
     * scoring algorithm
     * @throws IOException Reading exception
     * @throws InterruptedException Thread exception
     * @throws ExecutionException Execution exception
     */
	private static void runHCDcode() throws IOException,
            InterruptedException, ExecutionException {
		log.info("\nRunning in HCD mode.\n");

        int NCPU = Globals.numThreads;
        if(NCPU < Runtime.getRuntime().availableProcessors())
            NCPU = Globals.numThreads + 1;

		Globals.modelingMapHCD = new THashMap<>();
		
		// remove charge states you can't model
        Map<Integer, Integer> chargeMap = Globals.psmList.getChargeCount();
        AtomicInteger numPSM = new AtomicInteger(chargeMap.values().stream().reduce(0, Integer::sum));

        THashSet<Integer> badZ = new THashSet<>();
        log.info("PSMs for modeling:\n------------------\n");
        for(int z = 2; z <= Globals.maxChargeState; z++) {
			int n = chargeMap.get(z);
            log.info("+" + z + ": " + n + " PSMs");
            if(n < Globals.minNumPSMsForModeling)
                badZ.add(z);
		}
        log.info("\n");

        // Check to see if you have enough hi-scoring PSMs to construct an accurate model.
        // If not quit now
        if( (numPSM.get() < Globals.minNumPSMsForModeling) || (badZ.size() == chargeMap.size()) ) {
            log.info("You do not have enough PSMs with a score > " + Globals.modelTH +
                    " to accurately model the data. (Minimum number of PSMs required per charge state: " + Globals.minNumPSMsForModeling + ")\n" +
                    "Exiting now.\n"
            );
            System.exit(0);
        }

		
		// iterate over the charge states collecting data for all modeling PSMs
		log.info("(HCD) Acquiring Non-Parametric Model features from high-scoring PSMs...");

		for(int z = 2; z <= Globals.maxChargeState; z++) {

		    //We need a SynCronize List to be modified in the ParallelStream
		    final List<Peak> modelingPks = Collections
                    .synchronizedList(new ArrayList<>(1));

            // Skip this charge state you don't have enough data to model it
			if(!badZ.contains(z)) {
                numPSM.set(0); // track number of PSMs used for modeling

                int charge = z;
                Globals.psmList.parallelStream().forEach(p -> {
                    if ((p.getCharge() == charge) && p.isUseForModel()) {
                        p.generatePermutations(0); // generate real permutations
                        p.matchAllPeaks();
                    }
                });

                Globals.psmList.parallelStream().forEach(p-> {
                    if((p.getCharge() == charge) && p.isUseForModel()) {
                        if( (null != p.getPosPeaks()) && (!p.getPosPeaks().isEmpty()) ) {
                            modelingPks.addAll(p.getPosPeaks());
                            p.getPosPeaks().clear();
                        }
                        if( (null != p.getNegPeaks()) && (!p.getNegPeaks().isEmpty()) ) {
                            modelingPks.addAll(p.getNegPeaks());
                            p.getNegPeaks().clear();
                        }
                        numPSM.getAndIncrement();
                    }
                });

                if(!modelingPks.isEmpty()){
                    ModelDataHCD M = new ModelDataHCD(z, modelingPks);
                    M.numPSM = numPSM.get();
                    Globals.modelingMapHCD.put(z, M);

                    if(Globals.debugMode == Constants.WRITE_MODEL_PKS)
                        M.writeModelPks();
                }
            }
			

		} // end for loop over charge state for modeling


        if(Globals.modelingMapHCD.size() < 1) {
            log.info("\nInsufficient data to construct model.\nExiting now.\n");
            System.exit(0);
        }

		// Compute the statistics for the collected charge states
        ArrayList<Integer> missedChargeStates = new ArrayList<>();
        int maxObsZ = 0;  // highest modeled charge state.
		for(int z = 2; z <= Globals.maxChargeState; z++) {

			if(!Globals.modelingMapHCD.containsKey(z)) {
                missedChargeStates.add(z);
            }
			else {
                ModelDataHCD m = Globals.modelingMapHCD.get(z);

                m.calcMean();
                m.calcVar();

                m.printStats();

                // Compute intensity parameters
                m.estimateNP_intensity('b');
                m.estimateNP_intensity('y');
                m.estimateNP_intensity('n');

                // Compute distance parameters
                m.estimateNP_posDist();
                log.info("\n");  // makes for prettier output

                Globals.modelingMapHCD.put(z, m);
                m = null;
                if(z > maxObsZ) maxObsZ = z;
            }
		}

		
		// Back fill model data for missing charge states
        ModelDataHCD m = Globals.modelingMapHCD.get(maxObsZ);
        for(int missedZ : missedChargeStates) {
            Globals.modelingMapHCD.put(missedZ, m);
        }
        m = null;
        missedChargeStates.clear();

		// for debugging
		if(Globals.debugMode == Constants.WRITE_HCD_NONPARAM) {
			for(ModelDataHCD hcd : Globals.modelingMapHCD.values()) {
				hcd.write_density_data(1); // intensity
				hcd.write_density_data(2); // distance
			}
		}

        // Clean up before you start scoring
        System.gc();


        // Score the PSMs
        for(int RN = 0; RN < 2; RN++) { // RN = run number

            // If the number of threads is not defined use the default number of Processors, if
            // debug mode is use, then multithreading is disable due the writing needs of the peaks

            int numThreads = Globals.numThreads;
            if(numThreads <= 1 || Globals.debugMode != 0)
                numThreads =  1;
            if(numThreads > Runtime.getRuntime().availableProcessors())
                numThreads = Runtime.getRuntime().availableProcessors();

            if(Globals.runMode == Constants.DEFAULT_RUN_MODE)
                if(RN == 0)
                    log.info("\n[ " + RN + " ] Estimating FLR with decoys (" + numThreads + " threads)...");
                if(RN == 1)
                    log.info("\n[ " + RN + " ] Scoring " + Globals.psmList.size() +
                            " PSMs (" + numThreads + " threads)...");
             else
                log.info("\nScoring " + Globals.psmList.size() + " PSMs (" + numThreads + " threads)...");


            parallelFLP(numThreads, RN);
            System.err.println("\n");

            // first iteration is done, compute the FLR from the PSMs
            if(RN == 0) {
                Globals.calcFLR();

                if(Globals.runMode == Constants.DEFAULT_RUN_MODE) {
                    Globals.recordFLRestimates();
                    Globals.clearPSMs();
                }
                else { // you are done, exit loop over RN variable
                    RN = 2; // this will guarantee that you'll exit the loop
                    break;
                }
            }

            if(RN == 1) { // now assign FLRs to each PSM
                Globals.assignFLR();
            }
        } // end loop over RN
	}

    /**
     * This function compute the FLP using parallel stream and an a custom {@link ForkJoinPool}
     * @param numThreads Number of threads to be use in the {@link ForkJoinPool}
     * @param RN RN iteration
     * @throws ExecutionException
     * @throws InterruptedException
     */
	private static void parallelFLP(int numThreads, int RN) throws ExecutionException, InterruptedException {

        // Create a new Async ForJoinPool
	    ForkJoinPool forkJoinPool = new ForkJoinPool(numThreads,
                ForkJoinPool.defaultForkJoinWorkerThreadFactory, null, true);
	    forkJoinPool.submit(() -> {
            AtomicInteger ctr = new AtomicInteger(1);
            // When the number of threads is 1 do not use parallel Stream
            if(numThreads == 1){
                Globals.psmList.forEach(p-> {
                    p.generatePermutations(RN);
                    try {
                        p.scorePermutations();
                        if(Globals.debugMode == Constants.WRITE_SCORED_PKS) p.debugWriteScoredPeaks();
                        ctr.getAndIncrement();

                        if( (ctr.get() % 100) == 0 )  System.err.print(ctr + " ");
                        if( (ctr.get() % 1000) == 0 ) System.err.print("\n");
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                });
            }else{
                Globals.psmList.parallelStream().forEach(p-> {
                    p.generatePermutations(RN);
                    try {
                        p.scorePermutations();
                        if(Globals.debugMode == Constants.WRITE_SCORED_PKS) p.debugWriteScoredPeaks();
                        ctr.getAndIncrement();

                        if( (ctr.get() % 100) == 0 )  System.err.print(ctr + " ");
                        if( (ctr.get() % 1000) == 0 ) System.err.print("\n");
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                });
            }
        }).get();
    }
	
	
	/*
	 * Function writes the collected results to disk
	 */
	private static void writeResults() throws IOException {
		
		// set the output file's name
		if(Globals.outputFile.equalsIgnoreCase("luciphor_results.tsv")) {
			Globals.outputFile = "luciphor_results." + Globals.timeStamp + ".tsv";
		}

        File outF = new File(Globals.outputFile);
        log.info("\nResults written to '" + outF.getAbsoluteFile() + "'");

        // You will only need these variables when the user wants the matched
        // peaks to be written to disk
        File matchedPksF = null;
        FileWriter fwPks = null;
        BufferedWriter bwPks = null;

        if(Globals.writeMatchedPeaks) {
            String suffixDir = outF.getParent();
            Globals.matchedPkFile = suffixDir + "/luciphor_matchedPks." + Globals.timeStamp +".tsv";
            log.info("Matched peaks will be written to '" + Globals.matchedPkFile + "'\n");
            matchedPksF = new File(Globals.matchedPkFile);
            fwPks  = new FileWriter(matchedPksF.getAbsoluteFile());
            bwPks = new BufferedWriter(fwPks);

            // Column header for the peak file
            String pkHdr = "specId\tpepNum\tpredictPep\tfragmentIon\t" +
                           "m/z\trelIntensity\tDscore\tIscore\tscore\n";
            bwPks.write(pkHdr);
        }

        FileWriter fw = new FileWriter(outF.getAbsoluteFile());
		BufferedWriter bw = new BufferedWriter(fw);
		
		// Print column headers in output file
		String hdr = "specId\t"
				   + "peptide\t"
				   + "predictedPep1\t"
				   + "predictedPep2\t"
				   + "numPPS\t"
				   + "numRPS\t";
		
		switch(Globals.scoringMethod) {
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

        if(Globals.runMode == Constants.REPORT_DECOYS) hdr += "isDecoy1\tisDecoy2\t";

		hdr += "deltaScore\t"
			+ "pep1score\tpep2score\t"
			+ "globalFLR\tlocalFLR\n";
		
		bw.write(hdr);

		for(PSM psm : Globals.psmList) {
			bw.write( psm.getResults() );
            if(Globals.writeMatchedPeaks) bwPks.write( psm.writeMatchedPks() );
		}
		
		bw.close();
        if(Globals.writeMatchedPeaks) bwPks.close();
		
	}
	
	
}
