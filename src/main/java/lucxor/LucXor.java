/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import javax.xml.parsers.ParserConfigurationException;

import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.THashSet;
import gnu.trove.set.hash.TIntHashSet;
import lombok.extern.slf4j.Slf4j;
import lucxor.algorithm.FLR;
import lucxor.algorithm.ModelDataCID;
import lucxor.algorithm.ModelDataHCD;
import lucxor.algorithm.ModelParameterWorkerThread;
import lucxor.common.*;
import lucxor.pride.PrideJsonRandomAccess;
import lucxor.utils.Constants;
import lucxor.utils.MathFunctions;
import org.apache.commons.cli.*;
import org.xml.sax.SAXException;
import uk.ac.ebi.pride.archive.dataprovider.common.Tuple;
import uk.ac.ebi.pride.archive.dataprovider.data.ptm.IdentifiedModification;
import uk.ac.ebi.pride.archive.dataprovider.data.spectra.ArchiveSpectrum;
import uk.ac.ebi.pride.archive.dataprovider.param.CvParamProvider;
import umich.ms.datatypes.index.Index;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.LCMSDataSource;
import umich.ms.fileio.filetypes.mzml.MZMLFile;
import umich.ms.fileio.filetypes.mzxml.MZXMLFile;

import static lucxor.utils.Constants.*;

/**
 *
 * @author dfermin
 * @author ypriverol
 */

@Slf4j
public class LucXor {

    private static long startTime;
    private static long endTime;
    private static String timeStamp;

    private static PSMList psmList;
    private static TMap<Double, double[]> FLRestimateMap = new THashMap<>(); // ary[0] = globalFLR, ary[1] = localFLR

    private static TMap<Integer, ModelDataCID> modelingMapCID = null;
    private static TMap<Integer, ModelDataHCD> modelingMapHCD = null;

    private static PrideJsonRandomAccess prideJsonReader;

    public static void main(String[] args) throws
            IOException, IllegalStateException,
            InterruptedException, ExecutionException, FileParsingException {

        String prideInput = null;
        String prideOutput = null;

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


        CommandLine cmd;
        CommandLineParser parser = new BasicParser();
        HelpFormatter helper = new HelpFormatter();

        Options options = generateToolOptions();

        try {
            cmd = parser.parse(options, args);
            if(cmd.hasOption("-t")) {
                File outF = LucXorConfiguration.writeTemplateInputFile();
                log.info("\nPlease edit the input file: " + outF.getPath() +
                        " with your favorite text editor\n\n");
                System.exit(0);
            }

            if (cmd.hasOption("c")) {
                String opt_config = cmd.getOptionValue("config");
                // Read Configuration file
                LucXorConfiguration.parseConfigurationFile(opt_config);
                // Load User Modifications
                LucXorConfiguration.loadUserMods();
            }

            if (cmd.hasOption("a")){
                String opt_algorithm = cmd.getOptionValue("algorithm");
                LucXorConfiguration.setAlgorithm(Integer.parseInt(opt_algorithm));
            }

            if (cmd.hasOption("p")){
                prideInput = cmd.getOptionValue("prideFile");
                LucXorConfiguration.setSpectraFileType(PRIDE_TYPE);
                LucXorConfiguration.setScoringMethod(Constants.NEGLOGEXPECT);
            }

            if(cmd.hasOption("o")){
                prideOutput = cmd.getOptionValue("prideOutputFile");
            }

            if(prideInput != null && prideOutput == null)
                throw new ParseException("If pride json file is provided and output json file should be also provided");

        } catch (ParseException e) {
            System.out.println(e.getMessage());
            helper.printHelp("Usage:", options);
            System.exit(0);
        }

        // Start recording the running of the program
        startTime = System.nanoTime();

		// Time Stamp
        SimpleDateFormat sdf = new SimpleDateFormat("yyyyMMMdd-hh_mm_ss");
        timeStamp = sdf.format(new java.util.Date());

		// read in the score results
        psmList = new PSMList();
		if (prideInput != null)
            parsePrideJsonFile(prideInput);
        else if (LucXorConfiguration.getInputType() == Constants.PEPXML)
			parsePepXML();
		else
			parseTSVSrc();

        // Read in spectra for these PSMs
		if(prideInput != null)
            readInPrideSpectra();
        else
            readInSpectra();

		// Run CID Algorithm
        if(LucXorConfiguration.getScoringAlgorithm() == Constants.CID)
            runCIDcode();

        // Run HCD Algorithm
    	if(LucXorConfiguration.getScoringAlgorithm() == Constants.HCD)
    	    runHCDcode();

    	//Write final results
		writeResults();

        endTime = System.nanoTime();
        reportElapsedTime();
	}

    private static void readInPrideSpectra() throws IOException {

        log.info("\nReading PRIDE spectra ...");

        int assignedSpectraCtr = 0;

        // Assign the spectra to their respective PSMs
        for(PSM p : psmList) {
            ArchiveSpectrum spectrum = prideJsonReader.readArchiveSpectrum(p.getSpecId());
            Spectrum spec = new Spectrum(Stream.of(spectrum.getMasses()).mapToDouble(Double::doubleValue).toArray(), Stream.of(spectrum.getIntensities()).mapToDouble(Double::doubleValue).toArray());
            p.recordSpectra(spec);
            assignedSpectraCtr++;
        }
        log.info("Spectra readed in total: " + assignedSpectraCtr );

    }

    private static void parsePrideJsonFile(String prideInput) throws IOException {

        prideJsonReader = new PrideJsonRandomAccess(prideInput);
        prideJsonReader.parseIndex();
        for(String key: prideJsonReader.getKeys()){

            //mzspec:PRD000909:tio2_250ug:scan:92678:SQEPVTLDFLDAELENDIK/2
            Integer scanNumber = Integer.parseInt(key.split(":")[4]);
            String spectrumFile = key.split(":")[2];

            ArchiveSpectrum spectrum = prideJsonReader.readArchiveSpectrum(key);
            PSM psm = PSM.builder()
                    .charge(spectrum.getPrecursorCharge())
                    .isDecoy(spectrum.isDecoy())
                    .srcFile(spectrumFile)
                    .scanNum(scanNumber)
                    .specId(spectrumFile + "." + scanNumber + "." + scanNumber + "." + spectrum.getPrecursorCharge())
                    .build();

            double obsScore = Double.parseDouble( spectrum.getBestSearchEngineScore().getValue());
            psm.setPeptideSequence(spectrum.getPeptideSequence());

            psm.setPSMscore(-1.0 * Math.log(obsScore));

            if(spectrum.getModifications().size() == 0) { // no modifications so skip this PSM
                continue;
            }

            // modification string syntax: <pos>=<mass_of_modified_AA>
            // multiple modifications are comma separated (no spaces)
            // N-term pos = -100, C-term pos = 100
            // Fixed modifications are not reported
            for(IdentifiedModification modSites : spectrum.getModifications()) {
                for( Tuple<Integer, Set<? extends CvParamProvider >> entry: modSites.getPositionMap()){
                    int pos = entry.getKey() - 1;
                    if (pos == -1)
                        pos = -100; // N-Term modification
                    if (pos > spectrum.getPeptideSequence().length())
                        pos = 100;  // C-term modification
                    double mass = Double.parseDouble(modSites.getModification().getValue());
                    psm.getModCoordMap().put(pos, mass);
                }
            }

            // Determine if this PSM has any non-standard AA characters.
            // If so, discard it
            int numBadChars = 0;
            for(int i = 0; i < psm.getPeptideSequence().length(); i++) {
                String c = Character.toString( psm.getPeptideSequence().charAt(i) );
                if( !"ACDEFGHIKLMNPQRSTVWY".contains(c) ) numBadChars++;
            }

            // Skip PSMs that exceed the number of candidate permutations
            if(psm.getOrigPep().getNumPerm() > LucXorConfiguration.getMaxNumPermutations()) numBadChars = 100;

            if(numBadChars == 0) {
                psm.process();
                psm.setSpecId(key);
                if(psm.isKeeper()) psmList.add(psm);
            }
        }
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

        System.err.println("\nTotal run time (HH:MM:SS) = " +
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
	private static void parsePepXML() throws IOException, FileParsingException {

	    log.info("\nReading PSMs from pepXML file: " +
                LucXorConfiguration.getInputFile().getAbsolutePath());

	    if(LucXorConfiguration.getIdentificationSuffix() != null && LucXorConfiguration.getInputFile().isDirectory()){
             Stream<Path> walk = Files.walk(Paths.get(LucXorConfiguration.getInputFile().getAbsolutePath()));
             List<String> result = walk.map(Path::toString)
                        .filter(f -> f.toLowerCase().endsWith(LucXorConfiguration.getIdentificationSuffix().toLowerCase()))
                     .collect(Collectors.toList());

             result.parallelStream().forEach( fileName -> {
                 try {
                     PepXML.readStreamFile(new File(fileName), psmList);
                     log.info("The following file has been read it -- " + fileName);
                 } catch ( Exception e) {
                     log.error("Error reading input file -- " + fileName + " " + e.getMessage());
                 }
             });

        }else if(!LucXorConfiguration.getInputFile().isDirectory()){
            PepXML.readStreamFile(LucXorConfiguration.getInputFile(), psmList);
        }else{
	        throw new IOException("The provided identification file path or file is not valid");
        }
		log.info(psmList + " Candidate PSMs read in.");
	}

	
	// Function parses a tab-delimited file of PSMs
	private static void parseTSVSrc() throws IOException {
	    log.info("\nReading PSM from TSV file: " + LucXorConfiguration.getInputFile().getAbsolutePath());
		PSM curPSM;

		if( !LucXorConfiguration.getInputFile().exists() ) {
			log.info("ERROR: Unable to find " + LucXorConfiguration.getInputFile().getAbsolutePath());
			System.exit(0);
		}
		
		BufferedReader br = new BufferedReader(new FileReader(LucXorConfiguration.getInputFile()));
		String line;
        boolean passedHdr = false; // used only if the input file has a header line

		while( (line = br.readLine()) != null) {
			if(line.startsWith("#")) continue;

            // Skip over the header line for the TSV file (if it has one)
            if( (LucXorConfiguration.getTsvHdr() == 1) && (!passedHdr) ) {
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

			if(LucXorConfiguration.getScoringMethod() == Constants.NEGLOGEXPECT) {
			    curPSM.setPSMscore(-1.0 * Math.log(obsScore));
            } else {
			    curPSM.setPSMscore(obsScore);
            }


			if(vars.length < 6) { // no modifications so skip this PSM
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
			if(curPSM.getOrigPep().getNumPerm() > LucXorConfiguration.getMaxNumPermutations()) numBadChars = 100;

			if(numBadChars == 0) {
				curPSM.process();
				if(curPSM.isKeeper()) psmList.add(curPSM);
			}
        }
		br.close();

		log.info("Read in " + psmList.size() + " PSMs");
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

        int NCPU = LucXorConfiguration.getNumThreads();
        if(NCPU > 1) {
            if (NCPU < Runtime.getRuntime().availableProcessors()) NCPU = LucXorConfiguration.getNumThreads() + 1;
        }
        else NCPU = 1;


		modelingMapCID = new THashMap<>();
		int numPSM = 0; // track number of PSMs used for modeling

		// First make sure you have enough PSMs for each charge state to create
		// an accurate model.
        TIntIntHashMap chargeMap = new TIntIntHashMap();
		for(int z = 2; z <= LucXorConfiguration.getMaxChargeState(); z++) chargeMap.put(z,0);
		for(PSM p : psmList) {
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
		for(int z = 2; z <= LucXorConfiguration.getMaxChargeState(); z++) {
			int n = chargeMap.get(z);
            log.info("+" + z + ": " + n + " PSMs");

            if(n < LucXorConfiguration.getMinNumPsmModeling()) badZ.add(z);
		}
        log.info("\n");

        // Check to see if you have enough hi-scoring PSMs to construct an accurate model.
        // If not quit now
        if( (numPSM < LucXorConfiguration.getMinNumPsmModeling()) || (badZ.size() == chargeMap.size()) ) {
            log.info("You do not have enough PSMs with a score > " + LucXorConfiguration.getMODELTH() +
                " to accurately model the data. (Minimum number of PSMs required per charge state: " + LucXorConfiguration.getMinNumPsmModeling() + ")\n" +
                "Exiting now.\n"
            );
            System.exit(0);
        }

		
		// iterate over the charge states collecting data for all modeling PSMs
		log.info("(CID) Building Parametric Models from high-scoring PSMs...");
		for(int z = 2; z <= LucXorConfiguration.getMaxChargeState(); z++) {
			ArrayList<Peak> modelingPks = new ArrayList<>(100);
			
			if(badZ.contains(z)) continue; // skip this charge state, you don't have enough data to model it
			
			numPSM = 0; // track number of PSMs used for modeling

            if((LucXorConfiguration.getNumThreads() == 1) || (LucXorConfiguration.getDebugMode() != 0)) {

                for (PSM p : psmList) {
                    if ((p.getCharge() == z) && p.isUseForModel()) {
                        p.generatePermutations(0); // generate both real permutations
                        p.matchAllPeaks();
                    }
                }
            }
            else { // multi-threaded modeling

                // A pool of 'N' worker threads
                final int NTHREADS = LucXorConfiguration.getNumThreads();
                ExecutorService executor = Executors.newFixedThreadPool(NTHREADS);
                for (PSM p : psmList) {
                    if ((p.getCharge() == z) && p.isUseForModel()) {
                        Runnable worker = new ModelParameterWorkerThread(p);
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
            for(PSM p : psmList) {
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

            if(LucXorConfiguration.getDebugMode() == Constants.WRITE_MODEL_PKS) M.writeModelPks();

            M.setNumPSM(numPSM);
			modelingMapCID.put(z, M);
            modelingPks.clear();
        } // end loop over charge state for modeling


        if(modelingMapCID.size() < 1) {
            log.info("\nInsufficient data to construct model.\nExiting now.\n");
            System.exit(0);
        }

		// Compute the statistics for the collected charge states
        TIntArrayList missedChargeStates = new TIntArrayList();
        int maxObsZ = 0;  // highest modeled charge state.
		for(int z = 2; z <= LucXorConfiguration.getMaxChargeState(); z++) {
			if(!modelingMapCID.containsKey(z)) {
                missedChargeStates.add(z);
            }
			else {
                ModelDataCID m = modelingMapCID.get(z);
                m.calcMean();
                m.calcVar();
                m.printSummaryStats();
                m.clearArrays();
                modelingMapCID.put(z, m);
                if(z > maxObsZ) maxObsZ = z;
            }
		}

        // Back fill any missing charge states using the data for charge state 'maxObsZ'
        ModelDataCID m = modelingMapCID.get(maxObsZ);
        for(int missedZ : missedChargeStates.toArray()) {
            modelingMapCID.put(missedZ, m);
        }
        missedChargeStates.clear();

        System.gc(); // try and reclaim some memory

		// Score the PSMs
        for(int RN = 0; RN < 2; RN++) { // RN = run number, 0 = calc. FDR 1 = assign FDR estimates

            if(LucXorConfiguration.getRunMode() == Constants.DEFAULT_RUN_MODE) {
                if(RN == 0)
                    log.info("\n[ " + RN + " ] Estimating FLR with decoys (" + NCPU + " threads)...");
                if(RN == 1)
                    log.info("\n[ " + RN + " ] Scoring " + psmList.size() + " PSMs (" + NCPU + " threads)...");
            }
            else {
                log.info("\nScoring " + psmList.size() + " PSMs (" + NCPU + " threads)...");
            }

            AtomicInteger ctr = new AtomicInteger(1);

            int currentRN = RN;
            psmList.parallelStream().forEach(p -> {
                p.generatePermutations(currentRN);
                try {
                    p.scorePermutations(modelingMapHCD, modelingMapCID);
                    if(LucXorConfiguration.getDebugMode() == Constants.WRITE_SCORED_PKS) p.debugWriteScoredPeaks(timeStamp, modelingMapHCD, modelingMapCID);
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
               calcFLR();
                if(LucXorConfiguration.getRunMode() == Constants.DEFAULT_RUN_MODE) {
                    recordFLRestimates();
                    psmList.clearScores();
                }
                else { // you are done, exit loop over RN variable
                    RN = 2;
                    break;
                }
            }
            if(RN == 1) { // now assign FLRs to each PSM
                assignFLR();
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
	private static void runHCDcode() throws IOException, InterruptedException, ExecutionException {

	    log.info("\nRunning in HCD mode.\n");


		modelingMapHCD = new THashMap<>();
		
		// remove charge states you can't model
        Map<Integer, Integer> chargeMap = psmList.getChargeCount();
        AtomicInteger numPSM = new AtomicInteger(chargeMap.values().stream().reduce(0, Integer::sum));

        THashSet<Integer> badZ = new THashSet<>();
        log.info("PSMs for modeling:\n------------------\n");
        for(int z = 2; z <= LucXorConfiguration.getMaxChargeState(); z++) {
            if(!chargeMap.containsKey(z))
                badZ.add(z);
            else {
                int n = chargeMap.get(z);
                log.info("+" + z + ": " + n + " PSMs");
                if(n < LucXorConfiguration.getMinNumPsmModeling())
                    badZ.add(z);
            }
		}
        log.info("\n");

        // Check to see if you have enough hi-scoring PSMs to construct an accurate model.
        // If not quit now
        if( (numPSM.get() < LucXorConfiguration.getMinNumPsmModeling()) || (badZ.size() == chargeMap.size()) ) {
            log.info("You do not have enough PSMs with a score > " + LucXorConfiguration.getMODELTH() +
                    " to accurately model the data. (Minimum number of PSMs required per charge state: " + LucXorConfiguration.getMinNumPsmModeling() + ")\n" +
                    "Exiting now.\n"
            );
            System.exit(0);
        }

		
		// iterate over the charge states collecting data for all modeling PSMs
		log.info("(HCD) Acquiring Non-Parametric Model features from high-scoring PSMs...");
        System.err.println("(HCD) Acquiring Non-Parametric Model features from high-scoring PSMs...");

		for(int z = 2; z <= LucXorConfiguration.getMaxChargeState(); z++) {

		    //We need a Synchronized List to be modified in the ParallelStream
		    final List<Peak> modelingPks = Collections
                    .synchronizedList(new ArrayList<>(1));

            // Skip this charge state you don't have enough data to model it
			if(!badZ.contains(z)) {
                numPSM.set(0); // track number of PSMs used for modeling

                int charge = z;
                psmList.parallelStream().forEach(p -> {
                    if ((p.getCharge() == charge) && p.isUseForModel()) {
                        p.generatePermutations(0); // generate real permutations
                        p.matchAllPeaks();
                    }
                });

                psmList.parallelStream().forEach(p-> {
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
                    M.setNumPSM(numPSM.get());
                    modelingMapHCD.put(z, M);

                    if(LucXorConfiguration.getDebugMode() == Constants.WRITE_MODEL_PKS)
                        M.writeModelPks();
                }
            }
			

		} // end for loop over charge state for modeling


        if(modelingMapHCD.size() < 1) {
            log.info("\nInsufficient data to construct model.\nExiting now.\n");
            System.exit(0);
        }

		// Compute the statistics for the collected charge states
        ArrayList<Integer> missedChargeStates = new ArrayList<>(5);
        int maxObsZ = 0;  // highest modeled charge state.
		for(int z = 2; z <= LucXorConfiguration.getMaxChargeState(); z++) {

			if(!modelingMapHCD.containsKey(z)) {
                missedChargeStates.add(z);
            }
			else {
                ModelDataHCD m = modelingMapHCD.get(z);

                m.calcMean();
                m.calcVar();

                m.printStats();

                // Compute intensity parameters
                m.estimateNPIntensity('b');
                m.estimateNPIntensity('y');
                m.estimateNPIntensity('n');

                // Compute distance parameters
                m.estimateNPPosDist();
                log.info("\n");  // makes for prettier output

                modelingMapHCD.put(z, m);
                if(z > maxObsZ) maxObsZ = z;
            }
		}

		
		// Back fill model data for missing charge states
        ModelDataHCD m = modelingMapHCD.get(maxObsZ);
        for(int missedZ : missedChargeStates) {
            modelingMapHCD.put(missedZ, m);
        }
        missedChargeStates.clear();

		// for debugging
		if(LucXorConfiguration.getDebugMode() == Constants.WRITE_HCD_NONPARAM) {
			for(ModelDataHCD hcd : modelingMapHCD.values()) {
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

            int numThreads = LucXorConfiguration.getNumThreads();
            if(numThreads <= 1 || LucXorConfiguration.getDebugMode() != 0)
                numThreads =  1;
            if(numThreads > Runtime.getRuntime().availableProcessors())
                numThreads = Runtime.getRuntime().availableProcessors();

            if(LucXorConfiguration.getRunMode() == Constants.DEFAULT_RUN_MODE)
                if(RN == 0){
                    log.info("\n[ " + RN + " ] Estimating FLR with decoys (" + numThreads + " threads)...");
                    System.err.println("\n[ " + RN + " ] Estimating FLR with decoys (" + numThreads + " threads)...");
                }
                if(RN == 1){
                    log.info("\n[ " + RN + " ] Scoring " + psmList.size() +
                            " PSMs (" + numThreads + " threads)...");
                    System.err.println("\n[ " + RN + " ] Scoring " + psmList.size() +
                            " PSMs (" + numThreads + " threads)...");
                }
             else{
                 log.info("\nScoring " + psmList.size() + " PSMs (" + numThreads + " threads)...");
                 System.err.println("\nScoring " + psmList.size() + " PSMs (" + numThreads + " threads)...");
             }

            parallelFLP(numThreads, RN);
            System.err.println("\n");

            // first iteration is done, compute the FLR from the PSMs
            if(RN == 0) {
                calcFLR();

                if(LucXorConfiguration.getRunMode() == Constants.DEFAULT_RUN_MODE) {
                    recordFLRestimates();
                    psmList.clearScores();
                }
                else { // you are done, exit loop over RN variable
                    RN = 2; // this will guarantee that you'll exit the loop
                    break;
                }
            }

            if(RN == 1) { // now assign FLRs to each PSM
                assignFLR();
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
                psmList.forEach(p-> {
                    p.generatePermutations(RN);
                    try {
                        p.scorePermutations(modelingMapHCD, modelingMapCID);
                        if(LucXorConfiguration.getDebugMode() == Constants.WRITE_SCORED_PKS)
                            p.debugWriteScoredPeaks(timeStamp, modelingMapHCD, modelingMapCID);
                        ctr.getAndIncrement();

                        if( (ctr.get() % 100) == 0 ){
                            System.err.print(ctr + " ");
                            log.info(ctr + " ");
                        }if( (ctr.get() % 1000) == 0 ){
                            System.err.print("\n");
                        }
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                });
            }else{
                psmList.parallelStream().forEach(p-> {
                    p.generatePermutations(RN);
                    try {
                        p.scorePermutations(modelingMapHCD, modelingMapCID);
                        if(LucXorConfiguration.getDebugMode() == Constants.WRITE_SCORED_PKS)
                            p.debugWriteScoredPeaks(timeStamp, modelingMapHCD, modelingMapCID);
                        ctr.getAndIncrement();

                        if( (ctr.get() % 100) == 0 ) {
                            System.err.print(ctr + " ");
                            log.info(ctr + " ");
                        }
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

	    String outputFile = LucXorConfiguration.getOutputFile();

	    // set the output file's name
		if(LucXorConfiguration.getOutputFile().equalsIgnoreCase("luciphor_results.tsv")) {
			outputFile = "luciphor_results." + timeStamp + ".tsv";
		}

        File outF = new File(outputFile);
        log.info("\nResults written to '" + outF.getAbsoluteFile() + "'");

        // You will only need these variables when the user wants the matched
        // peaks to be written to disk
        File matchedPksF;
        FileWriter fwPks;
        BufferedWriter bwPks = null;

        if(LucXorConfiguration.isWriteMatchedPeaks()) {
            String suffixDir = outF.getParent();
            String matchedFileName = suffixDir + "/luciphor_matchedPks." + timeStamp +".tsv";
            log.info("Matched peaks will be written to '" + matchedFileName + "'\n");
            matchedPksF = new File(matchedFileName);
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
		
		switch(LucXorConfiguration.getScoringMethod()) {
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

        if(LucXorConfiguration.getRunMode() == Constants.REPORT_DECOYS) hdr += "isDecoy1\tisDecoy2\t";

		hdr += "deltaScore\t"
			+ "pep1score\tpep2score\t"
			+ "globalFLR\tlocalFLR\n";
		
		bw.write(hdr);

		for(PSM psm : psmList) {
			bw.write( psm.getResults() );
            if(LucXorConfiguration.isWriteMatchedPeaks()) bwPks.write( psm.writeMatchedPks(modelingMapCID, modelingMapHCD) );
		}
		
		bw.close();
        if(LucXorConfiguration.isWriteMatchedPeaks()) bwPks.close();
		
	}


    // Function to compute the False Localization Rate of the PSMs
    private static void calcFLR() throws InterruptedException, ExecutionException {
        double maxDeltaScore = -1.0;
        FLR flr = new FLR();

        log.info("\nComputing False Localization Rate (FLR)");

        // Identify maxDeltaScore
        for(PSM psm : psmList) {
            if(psm.getDeltaScore() > maxDeltaScore) maxDeltaScore = psm.getDeltaScore();

            if(psm.getDeltaScore() > Constants.MIN_DELTA_SCORE) {
                if(psm.isDecoy()) flr.addDecoy(psm);
                else flr.addTarget(psm);
            }
        }

        flr.setMaxDeltaScore(maxDeltaScore);
        flr.prepArrays();

        flr.initializeTickMarks();
        flr.evalTickMarks(Constants.REAL);
        flr.evalTickMarks(Constants.DECOY);

        flr.calcBothFDRs();
        flr.setMinorMaps();
        flr.performMinorization();
        flr.assignFDRs(psmList);
//		flr.debugFLR();
    }


    private static void readInSpectra() throws IOException, IllegalStateException,
            FileParsingException, ExecutionException, InterruptedException {

        log.info("\nReading spectra from " + LucXorConfiguration.getSpectrumPath()
                .getCanonicalPath() + "  (" + LucXorConfiguration.getSpectrumPrefix().toUpperCase()
                + " format)");

        log.info("This can take a while so please be patient.");

        Map<String, List<Integer>> scanMap;

        scanMap = psmList.parallelStream()
                .filter( p -> {
                    String pathStr = LucXorConfiguration.getSpectrumPath() + "/" + p.getSrcFile();
                    File f = new File(pathStr);
                    return f.exists();
                }).collect(Collectors.groupingBy(PSM::getSrcFile,
                        Collectors.mapping(PSM::getScanNum, Collectors.toList())));

        if(LucXorConfiguration.getSpectrumPrefix().equalsIgnoreCase(MGF_TYPE)) {
            TIntObjectHashMap<Spectrum> curSpectra;

            for(String specFile : scanMap.keySet()) {
                curSpectra = read_mgf(specFile);
                String fn = new File(specFile).getName(); // get just the file name of specFile
                int assignedSpectraCtr = 0;

                // Assign the spectra to their respective PSMs
                for(PSM p : psmList) {
                    if(p.getSrcFile().equalsIgnoreCase(fn)) {
                        if(curSpectra.containsKey(p.getScanNum())) {
                            p.recordSpectra( curSpectra.get(p.getScanNum()) );
                            assignedSpectraCtr++;
                        }
                    }
                }
                log.info(fn + ": " + assignedSpectraCtr + " spectra read in.");
            }
        }

        // Read mzXML files
        if(LucXorConfiguration.getSpectrumPrefix().equalsIgnoreCase(MZXML_TYPE))
            readXMLFile(scanMap, LucXorConfiguration.getSpectrumPath().getAbsolutePath(),MZXML_TYPE);

        // Read mzML files
        if(LucXorConfiguration.getSpectrumPrefix().equalsIgnoreCase(MZML_TYPE))
            readXMLFile(scanMap, LucXorConfiguration.getSpectrumPath().getAbsolutePath(),
                    MZML_TYPE);
    }

    /**
     * Function reads in an MGF file and returns it as a HashMap
     * @param specFile Spectrum File
     * @return Spectrum Map
     * @throws IOException
     */
    private static TIntObjectHashMap<Spectrum> read_mgf(String specFile) throws
            IOException {
        TIntObjectHashMap<Spectrum> ret = new TIntObjectHashMap<>();

        File mgf = new File(specFile);
        BufferedReader br = new BufferedReader(new FileReader(mgf));
        String line;
        int scanNum = 0;
        Spectrum S;
        ArrayList<Double> mzAL = null, intensityAL = null;


        while( (line = br.readLine()) != null ) {

            if(line.length() < 2) continue;

            if(line.startsWith("END IONS")) {
                if( (null != mzAL) && (!mzAL.isEmpty()) ) {

                    int N = mzAL.size();
                    double[] mz = new double[ N ];
                    double[] I = new double[ N ];
                    for(short i = 0; i < N; i++) {
                        mz[i] = mzAL.get(i);
                        I[i] = intensityAL.get(i);
                    }
                    S = new Spectrum(mz, I);
                    ret.put(scanNum, S);

                    mzAL = null;
                    intensityAL = null;
                }
                scanNum = 0;
            }

            if(line.startsWith("BEGIN IONS")) {
                mzAL = new ArrayList<>(100);
                intensityAL = new ArrayList<>(100);
            }

            if(line.startsWith("TITLE=")) {
                int i = line.indexOf('.') + 1;
                int j = line.indexOf('.', i);
                String s = line.substring(i,j);
                scanNum = Integer.valueOf(s);
            }

            if(line.startsWith("CHARGE") || line.startsWith("PEPMASS")) continue;

            if( Character.isDigit( line.charAt(0) ) ) {
                String[] ary = line.split("\\s+");
                double mz = MathFunctions.roundDouble( Double.valueOf(ary[0]), 8 );
                double I  = MathFunctions.roundDouble( Double.valueOf(ary[1]), 8 );
                mzAL.add(mz);
                intensityAL.add(I);
            }
        }

        return ret;
    }

    private static void readXMLFile(Map<String, List<Integer>> scanMap, String
            pathSpectra, String fileType) throws FileParsingException, ExecutionException, InterruptedException {

        long timeLo = System.nanoTime();

        for(Map.Entry<String, List<Integer>> stringListEntry : scanMap.entrySet()) {
            String fileName = stringListEntry.getKey();
            String absolutePath = new File(pathSpectra + "/" + stringListEntry.getKey()).getAbsolutePath();

            LCMSDataSource<?> source;
            switch (fileType) {
                case MZML_TYPE:
                    source = new MZMLFile(absolutePath);
                    break;
                case MZXML_TYPE:
                    source = new MZXMLFile(absolutePath);
                    break;
                default:
                    throw new UnsupportedOperationException("file not supported");
            }
            // create data structure to hold scans and load all scans
            ScanCollectionDefault scans = new ScanCollectionDefault();
            scans.setDataSource(source);

            Index<?> index = source.parseIndex();
            index.getMapByNum();
            final int batchSize = 100;

            Iterator<Integer> it = stringListEntry.getValue().iterator();
            List<IScan> parsed = new ArrayList<>(stringListEntry.getValue().size());
            while (it.hasNext()) {
                List<Integer> scanNums = new ArrayList<>(batchSize);
                scanNums.add(it.next());
                while (it.hasNext() && scanNums.size() < batchSize) {
                    scanNums.add(it.next());
                }
                parsed.addAll(source.parse(scanNums));
                System.err.print("\r" + absolutePath + ":  " + parsed.size() + "... ");
            }

            // data index, can be used to locate scans by numbers or retention times at different ms levels

            // Create a new Async ForJoinPool
            ForkJoinPool forkJoinPool = new ForkJoinPool(LucXorConfiguration.getNumThreads(),
                    ForkJoinPool.defaultForkJoinWorkerThreadFactory, null, true);
            forkJoinPool.submit(() -> parsed.parallelStream().forEach(kv -> {

                final ISpectrum spectrum = kv.getSpectrum();
                int N = spectrum.getMZs().length;
                if(!(N == 0)) {
                    double[] mz = spectrum.getMZs();
                    double[] intensities = spectrum.getIntensities();
                    // If this happens, there is something wrong with the spectrum so skip it
                    if(mz.length == intensities.length) {
                        Spectrum X = new Spectrum(mz, intensities);
                        PSM psm = psmList.getByScanOrder(fileName, kv.getNum());
                        psm.recordSpectra(X);
                    }
                }
            })).get();
        }
        long timeHi = System.nanoTime();
        System.out.printf("\n\nLoading all the PSM spectra %d took %.1fs\n\n",
                psmList.size(), (timeHi - timeLo)/1e9f);


    }

    // Function assigns the global and local FLR for the current PSM from the FLRestimateMap
    private static void assignFLR() {

        ArrayList<Double> obsDeltaScores = new ArrayList<>(FLRestimateMap.keySet());
        Collections.sort(obsDeltaScores);  // sort them from low to high
        int N = obsDeltaScores.size();
        boolean assigned;
        for(PSM p : psmList) {
            double obs_ds = p.getDeltaScore();
            assigned = false;

            // iterate over the delta scores until you find the value closest to this one
            int i;
            for(i = 1; i < N; i++) {
                double curDS = obsDeltaScores.get(i);
                if(curDS > obs_ds) { // hit the limit, get the *previous* delta score
                    double[] d = FLRestimateMap.get( obsDeltaScores.get((i-1)) );
                    p.setGlobalFDR(d[0]);
                    p.setLocalFDR(d[1]);
                    assigned = true;
                    break;
                }
            }

            if(!assigned) { // very high scoring PSM
                double[] d = FLRestimateMap.get( obsDeltaScores.get((N-1)) );
                p.setGlobalFDR(d[0]);
                p.setLocalFDR(d[1]);
            }
        }
    }

    /**
     * Record the global and local FLR values estimated for all of the delta scores
     */
    private static void recordFLRestimates() {
        FLRestimateMap = new THashMap<>();

        for(PSM p : psmList) {
            if(p.isDecoy()) continue; // skip FLR data from decoys
            double[] d = new double[2];
            d[0] = p.getGlobalFDR();
            d[1] = p.getLocalFDR();
            FLRestimateMap.put( p.getDeltaScore(), d );
        }
    }

    private static Options generateToolOptions(){
        Options options = new Options();

        Option confFile = Option.builder("t")
                .longOpt("generation").hasArg()
                .required(false)
                .hasArg(false)
                .desc("generates a configuration file default")
                .build();
        options.addOption(confFile);

        Option config = Option.builder("c").longOpt("config")
                .argName("config")
                .hasArg(true)
                .required(true)
                .desc("config file to be used during the run").build();
        options.addOption(config);

        Option prideSpectra = Option.builder("p")
                .longOpt("prideFile")
                .hasArg(true)
                .required(false)
                .desc("use pride json file containing spectra")
                .build();
        options.addOption(prideSpectra);

        Option prideSpectraOutput = Option.builder("o")
                .longOpt("prideOutputFile")
                .hasArg(true)
                .required(false)
                .desc("output json file containing spectra and FLR")
                .build();
        options.addOption(prideSpectraOutput);

        Option algorithm = Option.builder("a")
                .longOpt("algorithm")
                .hasArg(true)
                .required(false)
                .desc("algorithm to be used by the tool HCD or CID")
                .build();
        options.addOption(algorithm);

        return options;
    }
}
