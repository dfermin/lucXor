/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import gnu.trove.map.TMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import lombok.extern.slf4j.Slf4j;
import org.xml.sax.SAXException;
import umich.ms.datatypes.LCMSData;
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scancollection.IScanCollection;
import umich.ms.datatypes.scancollection.ScanIndex;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.mzml.MZMLFile;
import umich.ms.fileio.filetypes.mzxml.MZXMLFile;

import javax.xml.parsers.ParserConfigurationException;
import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;
import java.util.zip.DataFormatException;

import static lucxor.Constants.*;

/**
 *
 * @author dfermin
 */
@Slf4j
public class Globals {

	static File spectrumPath = null;
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
	static double decoyMass;
	static double minMZ;
	static double max_num_permutations;
	static double precursorNLmass;
	static double ntermMass = 0d;
	static double ctermMass = 0d;
	static double minRelIntensity;
    static boolean writeMatchedPeaks;
	
	static PSMList psmList = new PSMList();

    static TMap<String, Double> TargetModMap = new THashMap<>(); // mods user wants to search for
	static TMap<String, Double> fixedModMap = new THashMap<>(); // fixed mods observed in data
	static TMap<String, Double> varModMap = new THashMap<>(); // variable mods observed in data
	static TMap<String, Double> nlMap = new THashMap<>(); // holds all neutral loss masses, k= list of amino acids v= NL mass
	static TMap<String, Double> decoyNLmap = new THashMap<>();
    static TMap<Double, double[]> FLRestimateMap = new THashMap<>(); // ary[0] = globalFLR, ary[1] = localFLR

	static TMap<Integer, ModelDataCID> modelingMap_CID = null;
	static TMap<Integer, ModelDataHCD> modelingMap_HCD = null;
	
	static void parseInputFile(String str) throws IOException {
		
		File inF = new File(str);
		if(!inF.exists()) { 
			System.err.print("\nERROR! Unable to open " + inF.getAbsolutePath() + "\n\n");
			System.exit(0);
		}
		
		debugMode = 0; // 0 means no debugging output
        runMode = 0; // calculate FLR and rescore PSMs without decoys (2 iterations)
		minRelIntensity = 0;
		minNumPSMsForModeling = 50;
        maxPepLen = 40;
        reduceNL = 0;
		numThreads = Runtime.getRuntime().availableProcessors();
		tsvHdr = 0; // default is no header
        writeMatchedPeaks = false; // true means you will generate the matched peaks

		BufferedReader br = new BufferedReader(new FileReader(inF));
		String line;
		while( (line = br.readLine()) != null) {
			if(line.startsWith("#")) continue;
			if(line.length() < 2) continue;
			
			if(line.startsWith("SPECTRUM_PATH")) {
				String s = Utils.parseInputLine(line);
				spectrumPath = new File(s).getCanonicalFile();
			}
			
			if(line.startsWith("SPECTRUM_SUFFIX")) {
				String s = Utils.parseInputLine(line);
				spectrumSuffix = s.toLowerCase();
			}
			
			if(line.startsWith("INPUT_DATA")) {
				String s = Utils.parseInputLine(line);
				inputFile = new File(s);
			}
			
			if(line.startsWith("OUTPUT_FILE")) {
				outputFile = Utils.parseInputLine(line);
			}
			
			if(line.startsWith("INPUT_TYPE")) {
				String s = Utils.parseInputLine(line);
				inputType = Integer.valueOf(s);
			}

            if(line.startsWith("MAX_PEP_LEN")) {
                String s = Utils.parseInputLine(line);
                maxPepLen = Integer.valueOf(s);
            }

			if(line.startsWith("MAX_NUM_PERM")) {
				String s = Utils.parseInputLine(line);
				max_num_permutations = Double.valueOf(s);
			}
			
			if(line.startsWith("MIN_NUM_PSMS_MODEL")) {
				String s = Utils.parseInputLine(line);
				minNumPSMsForModeling = Integer.valueOf(s);
			}
			
			if(line.startsWith("MS2_TOL") && !line.contains("_UNITS")) {
				String s = Utils.parseInputLine(line);
				ms2tol = Double.valueOf(s);
			}
			
			if(line.startsWith("MS2_TOL_UNITS")) {
				String s = Utils.parseInputLine(line);
				ms2tol_units = Integer.valueOf(s);
			}
			
			if(line.startsWith("ALGORITHM")) {
				String s = Utils.parseInputLine(line);
				scoringAlgorithm = Integer.valueOf(s);
			}

            if(line.startsWith("TSV_HEADER")) {
                String s = Utils.parseInputLine(line);
                tsvHdr = Integer.valueOf(s);
            }

			if(line.startsWith("REDUCE_PRECURSOR_NL")) {
				String s = Utils.parseInputLine(line);
				reduceNL = Integer.valueOf(s);
			}
			
			if(line.startsWith("PRECURSOR_NL_MASS_DIFF")) {
				String s = Utils.parseInputLine(line);
				precursorNLmass = Double.valueOf(s);
			}
			
			if(line.startsWith("SELECTION_METHOD")) {
				String s = Utils.parseInputLine(line);
				scoringMethod = Integer.valueOf(s);
			}
			
			if(line.startsWith("MODELING_SCORE_THRESHOLD")) {
				String s = Utils.parseInputLine(line);
				modelTH = Double.valueOf(s);
			}
			
			if(line.startsWith("MAX_CHARGE_STATE")) {
				String s = Utils.parseInputLine(line);
				maxChargeState = Integer.valueOf(s);
			}
			
			if(line.startsWith("NUM_THREADS")) {
				String s = Utils.parseInputLine(line);
				int x = Integer.valueOf(s);
				// We do minus 1 because 1 thread already goes to 
				// running the whole program
				if(x < 0) numThreads = 1;
				else if(x > 1) numThreads = (x - 1);
                else if(x == 0) numThreads = Runtime.getRuntime().availableProcessors();
                else numThreads = x;
			}
			
			if(line.startsWith("DEBUG_MODE")) {
				String s = Utils.parseInputLine(line);
				debugMode = Integer.valueOf(s);
			}

            if(line.startsWith("WRITE_MATCHED_PEAKS_FILE")) {
                String s = Utils.parseInputLine(line);
                if(s.equals("1")) writeMatchedPeaks = true;
            }

            if(line.startsWith("RUN_MODE")) {
                String s = Utils.parseInputLine(line);
                runMode = Integer.valueOf(s);
                if(runMode > 1) runMode = 0;
            }
			
			if(line.startsWith("SCORING_THRESHOLD")) {
				String s = Utils.parseInputLine(line);
				scoreTH = Double.valueOf(s);
			}
			
			if(line.startsWith("DECOY_MASS")) {
				String s = Utils.parseInputLine(line);
				decoyMass = Double.valueOf(s);
			}

            if(line.startsWith("DECOY_NL")) {
                String[] ary = parse_NL_line(line);
                String k = ary[0].substring(1);
                double m = Double.valueOf(ary[1]);
                decoyNLmap.put(k,m);
            }
			
			if(line.startsWith("NL")) {
				String[] ary = parse_NL_line(line);
				double m = Double.valueOf(ary[1]);
				nlMap.put(ary[0], m);
			}
			
			if(line.startsWith("MIN_MZ")) {
				String s = Utils.parseInputLine(line);
				minMZ = Double.valueOf(s);
			}

			if(line.startsWith("MOD_PEP_REP")) {
				String s = Utils.parseInputLine(line);
				peptideRepresentation = Integer.valueOf(s); 
			}
			
			if(line.startsWith("TARGET_MOD")) {
				String[] ary = parse_input_mod_line(line);
				double m = Double.valueOf(ary[1]);
				TargetModMap.put(ary[0].toUpperCase(), m);
			}


            if(line.startsWith("VAR_MOD")) {
                String[] ary = parse_input_mod_line(line);
                double m = Double.valueOf(ary[1]);
                varModMap.put(ary[0].toLowerCase(), m); // variable mods are lower case
            }

            if(line.startsWith("FIXED_MOD")) {
                String[] ary = parse_input_mod_line(line);
                double m = Double.valueOf(ary[1]);
                fixedModMap.put(ary[0].toUpperCase(), m);
            }


			// You only need to extract the VAR_MOD and FIXED_MOD values
			// from the input file if you are using TSV files
//			if(Globals.inputType == Constants.TSV) {
//				if(line.startsWith("VAR_MOD")) {
//					String[] ary = parse_input_mod_line(line);
//					double m = Double.valueOf(ary[1]);
//					varModMap.put(ary[0].toLowerCase(), m); // variable mods are lower case
//				}
//
//				if(line.startsWith("FIXED_MOD")) {
//					String[] ary = parse_input_mod_line(line);
//					double m = Double.valueOf(ary[1]);
//					fixedModMap.put(ary[0].toUpperCase(), m);
//				}
//			}

		}
		br.close();
		
		if( (null == outputFile) || (outputFile.isEmpty()) ) outputFile = "luciphor_results.tsv";
		
		String classStr = "";
		switch(scoringMethod) {
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
				System.err.print("\nERROR! Unknown scoring method: " + scoringMethod+ "\n\n");
				System.exit(0);
				break;
		}

        int NCPU = numThreads;
        if( (numThreads > 1) && (numThreads < Runtime.getRuntime().availableProcessors()) ) NCPU = (numThreads + 1);

		log.info("Spectrum Path:           " + spectrumPath.getAbsolutePath());
		log.info("Spectrum Suffix:         " + spectrumSuffix);
		log.info("Input file:              " + inputFile);
		log.info("Input type:              " + (inputType == Constants.PEPXML ? "pepXML" : "tsv"));
		log.info("MS2 tolerance:           " + ms2tol + (ms2tol_units == Constants.DALTONS ? " Da" : " ppm"));
		log.info("Luciphor Algorithm:      " + (scoringAlgorithm == Constants.CID ? "CID" : "HCD") );
		log.info("Classifying on:          " + classStr);
        log.info("Run Mode:                " + (runMode == 0 ? "Default" : "Report Decoys"));
		log.info("Num of Threads:          " + NCPU );
		log.info("Modeling Threshold:      " + modelTH);
		log.info("Scoring Threshold:       " + scoreTH);
		log.info("Permutation Limit:       " + max_num_permutations);
        log.info("Max peptide length:      " + maxPepLen);
		log.info("Min num PSMs for model:  " + minNumPSMsForModeling);
		log.info("Decoy Mass Adduct:       " + decoyMass);
		log.info("Max Charge State:        " + maxChargeState);
		log.info("Reduce NL:               " + (reduceNL == 0 ? "no" : "yes"));
		log.info("Output File:             " + outputFile);
        log.info("Write matched Peaks:     " + (writeMatchedPeaks ? "yes" : "no"));
        System.err.print("\n");

		if(debugMode != 0) {
			log.info("Debug mode:              " + debugMode + "  (Limiting to 1 CPU)\n");
		    Globals.numThreads = 1;
        }
		
		log.info("Mods to score:");
		for(String s : TargetModMap.keySet()) {
			log.info(s + "\t" + TargetModMap.get(s));
		}
		
		if(!nlMap.isEmpty()) {
			log.info("\nAllowed Neutral Losses:");
			for(String s : nlMap.keySet()) {
				log.info(s + "\t" + nlMap.get(s));
			}
            for(String s : decoyNLmap.keySet()) {
                log.info("<X>" + s + "\t" + decoyNLmap.get(s) + "  (Decoy NL)");
            }
		}

		// If the scores are going to be e-values, convert them to -log(evalues)
		if(scoringMethod == Constants.NEGLOGEXPECT) {
			double tmp = -1.0 * Math.log( modelTH );
			modelTH = tmp;

			tmp = -1.0 * Math.log(scoreTH);
			scoreTH = tmp;
		}
		
	}

	
	static String[] parse_input_mod_line(String line) {
		String[] ret = new String[2];
		
		char aa = 0;
		StringBuilder sb = new StringBuilder();
		double mass = 0d;
		int N = line.length();
		
		int b = line.indexOf("=") + 1;
		
		for(int i = b; i < N; i++) {
			char c = line.charAt(i);
			
			if(c == '#') break;
			if(c == ' ') continue;
			
			if(Character.isAlphabetic(c)) aa = c;
			else if( (c == '[') || (c == ']') ) aa = c;
			else sb.append(c);
		}
		
		mass = Double.valueOf( sb.toString() );
		
		ret[0] = Character.toString(aa);
		ret[1] = String.valueOf(mass);
		
		return ret;
	}
	
	
	static String[] parse_NL_line(String line) {
		String[] ret = new String[2];
		
		line = line.replaceAll("#", "");
		
		String[] tmp = line.split("\\s+");
				
		ret[0] = tmp[2] + tmp[3];
		ret[1] = tmp[4];
		
		return ret;
	}

	/*
	 * In Java 8, you can call to an initilization method
	 */
	static {
		// in case we need it, initialize the timestamp signature variable
		timeStamp = "";
		java.util.Date date = new java.util.Date();
		SimpleDateFormat sdf = new SimpleDateFormat("yyyyMMMdd-hh_mm_ss");
		timeStamp = sdf.format(date);
		sdf = new SimpleDateFormat("yyyMMMdd");
		dateStamp = sdf.format(date);
	}

	// Function to compute the False Localization Rate of the PSMs
	public static void calcFLR() throws InterruptedException, ExecutionException {
		double maxDeltaScore = -1.0;
		FLR flr = new FLR();
		ArrayList<FLR> flrAry = new ArrayList<>();

		log.info("\nComputing False Localization Rate (FLR)");

		// Identify maxDeltaScore
		for(PSM psm : Globals.psmList) {
			if(psm.getDeltaScore() > maxDeltaScore) maxDeltaScore = psm.getDeltaScore();

			if(psm.getDeltaScore() > Constants.MIN_DELTA_SCORE) {
				if(psm.isDecoy()) flr.decoyPSMs.add(psm);
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


	/**
	 * Function to incorporate information read in from input file into the
	 * AA_MASS_MAP object
	 */
	public static void loadUserMods() {

		// Assign the information you got from the input file to the 
		// relevant map. If none were given, then this should just proceed
		// without changing anything.
		for(String c : TargetModMap.keySet()) {
			double mass = AA_MASS_MAP.get(c) + TargetModMap.get(c);
			String symbol = c.toLowerCase();
			AA_MASS_MAP.put(symbol, mass);
		}
		
		for(String c : fixedModMap.keySet()) {

		    if( c.equalsIgnoreCase("[") ) {
                Globals.ntermMass = fixedModMap.get(c);
            }
            else if( c.equalsIgnoreCase("]") ) {
                Globals.ctermMass = fixedModMap.get(c);
            }
		    else {
                double mass = AA_MASS_MAP.get(c) + fixedModMap.get(c);
                String symbol = c.toUpperCase();
                AA_MASS_MAP.put(symbol, mass);
            }
		}
		
		for(String c : varModMap.keySet()) { // this will be a lower case key
			
			if( c.equalsIgnoreCase("[") ) {
				Globals.ntermMass = varModMap.get(c);
			}
			else if( c.equalsIgnoreCase("]") ) {
				Globals.ctermMass = varModMap.get(c);
			}
			else {
				double mass = AA_MASS_MAP.get(c.toUpperCase()) + varModMap.get(c);
				AA_MASS_MAP.put(c, mass);
			}
		}
		
		
		// Now add the decoy masses to the AA_MASS_MAP
		for(String c : DECOY_AA_MAP.keySet()) {
			String trueAA = DECOY_AA_MAP.get(c);
			
			if(varModMap.containsKey(trueAA)) continue;
			if(TargetModMap.containsKey(trueAA)) continue;
			
			double mass = AA_MASS_MAP.get(trueAA) + Globals.decoyMass;
			AA_MASS_MAP.put(c, mass);
		}

	}


	/**
	 * This function is only called if we are getting our modifications from
	 * a pepXML file.
	 */
	static void recordModsFromPepXML() {
		
		String alphabet = "ACDEFGHIKLMNPQRSTVWY";
		
		for(String c : fixedModMap.keySet()) {
			
			if(!alphabet.contains(c)) continue; // skip non-amino acid characters 
		
			double mass = AA_MASS_MAP.get(c) + fixedModMap.get(c);
			String symbol = c.toUpperCase();
			AA_MASS_MAP.put(symbol, mass);
		}
		
		// First filter the data in varModMap.
		// We want to remove non-standard amino acid characters and remove
		// the amino acids that are in our 'TargetModMap' variable.
		HashMap<String, Double> tmp = new HashMap<>(varModMap);
		varModMap.clear();
		
		for(String c : tmp.keySet()) {
			String C = c.toUpperCase();
			double mass = tmp.get(c);
			if(!alphabet.contains(C)) continue; // skip non-amino acid characters
			if(TargetModMap.containsKey(C)) continue; // skip target residues
			varModMap.put(c, mass);
		}
		tmp = null;
		
		for(String c : varModMap.keySet()) {
			String C = c.toUpperCase();
			double mass = AA_MASS_MAP.get(C) + varModMap.get(c);
			AA_MASS_MAP.put(c, mass);
		}
	}


	static void readInSpectra() throws IOException, IllegalStateException,
			FileParsingException {
		
		log.info("\nReading spectra from " + Globals.spectrumPath
				.getCanonicalPath() + "  (" + Globals.spectrumSuffix.toUpperCase()
				+ " format)");

		log.info("This can take a while so please be patient.");

		Map<String, List<Integer>> scanMap;
		
		scanMap = psmList.parallelStream()
				.filter( p -> {
					String pathStr = Globals.spectrumPath + "/" + p.getSrcFile();
					File f = new File(pathStr);
					return f.exists();
				}).collect(Collectors.groupingBy(PSM::getSrcFile,
						Collectors.mapping(PSM::getScanNum, Collectors.toList())));

		if(Globals.spectrumSuffix.equalsIgnoreCase(MGF_TYPE)) {
			TIntObjectHashMap<SpectrumClass> curSpectra = null;
			
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
		if(Globals.spectrumSuffix.equalsIgnoreCase(MZXML_TYPE))
			readMzXML(scanMap, spectrumPath.getAbsolutePath());

		// Read mzML files
		if(Globals.spectrumSuffix.equalsIgnoreCase(MZML_TYPE))
			readMzML(scanMap, spectrumPath.getAbsolutePath());
	}


    /****************
     * Function reads in spectral data from mzML files
     * @param scanMap
     */
    private static void readMzML(Map<String, List<Integer>> scanMap, String spectraPath) throws
			FileParsingException {

		long timeLo = System.nanoTime();

        // Iterate over the file names
        for(String fn : scanMap.keySet()) {
            String baseFN = new File(spectraPath + "/" + fn).getName();
            System.err.print("\n" + baseFN + ":  "); // beginning of info line

            int ctr = 0;
            int iter = 0;
            List<Integer> scanNums = scanMap.get(fn);
            Collections.sort(scanNums); // order scan numbers

            // read in the mzXML file
            String mzML_path = Globals.spectrumPath + "/" + baseFN;

			final MZMLFile curMZML = new MZMLFile(mzML_path);

            for(int scanNum : scanNums) {
                iter++;
                if( iter % 100 == 0 ) {
                    System.err.print("\r" + baseFN + ":  " + iter + "... ");
                    // beginning of info line
                }
				final IScan scan = curMZML.parseScan(scanNum, true);
				final ISpectrum spectrum = scan.getSpectrum();
				int N = spectrum.getMZs().length;
                if(N == 0) {
                    continue; // no valid spectrum for this scan number
                }

                double[] mz = spectrum.getMZs();
                double[] intensities = spectrum.getIntensities();

                // If this happens, there is something wrong with the spectrum so skip it
                if(mz.length != intensities.length) {
                    System.err.print(
                            "\nERROR:" + baseFN + " Scan: " + scanNum +
                            "\n# of mz values != # intensity values: " +
                            mz.length + " != " + intensities.length +
                            "\nSkipping this scan...\n"
                    );
                    continue;
                }

                SpectrumClass X = new SpectrumClass(mz, intensities);

                PSM psm = psmList.getByScanOrder(baseFN, scanNum);
                psm.recordSpectra(X);
                ctr++;
            }
			// end of file reading
            System.err.print("\r" + baseFN +  ":  " + ctr + " spectra read in.            ");
        }
		long timeHi = System.nanoTime();
        log.info("Loading took %.1fs", (timeHi - timeLo)/1e9f);

    }

//	/****************
//	 * Function reads in spectral data from mzML files
//	 * @param scanMap
//	 */
//	private static void readXMLFile(Map<String, List<Integer>> scanMap,
//									String spectrumPath, String fileType) throws FileParsingException {
//
//		scanMap.entrySet().stream().forEach( fileEntry -> {
//
//			AtomicInteger ctr = new AtomicInteger();
//			SimpleDateFormat fmt= new SimpleDateFormat("HH:mm:ss");
//
//			String fileName = fileEntry.getKey();
//			LCMSDataSource<?> source = null;
//
//
//			if(fileType.equalsIgnoreCase(MZML_TYPE))
//				source = new MZMLFile(spectrumPath + "/" + fileName);
//
//			if(fileType.equalsIgnoreCase(MZXML_TYPE))
//				source = new MZXMLFile(spectrumPath + "/" + fileName);
//
//			if(source == null)
//				new IOException("Error, file Type not supported");
//
//			// create data structure to hold scans and load all scans
//			ScanCollectionDefault scans = new ScanCollectionDefault();
//			scans.setDataSource(source);
//
//			try {
//				System.out.printf("Start loading whole file @ [%s]\n", fmt.format(new Date()));
//				long timeLo = System.nanoTime();
//
//				scans.loadData(LCMSDataSubset.WHOLE_RUN);
//				long timeHi = System.nanoTime();
//
//				System.out.printf("Done loading whole file @ [%s]\n", fmt.format(new Date()));
//				System.out.printf("Loading took %.1fs", (timeHi - timeLo)/1e9f);
//
//				// data index, can be used to locate scans by numbers or retention times at different ms levels
//				TreeMap<Integer, ScanIndex> index = scans.getMapMsLevel2index();
//
//				// iterate over MS2 scnas asynchronously, and calculate total intensity
//				ScanIndex ms2scansIndex = index.get(2);
//				if (ms2scansIndex == null || ms2scansIndex.getNum2scan().isEmpty())
//					throw new IllegalStateException("empty ms2 index");
//
//				ExecutorService exec = Executors
//						.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
//
//
//				Set<Map.Entry<Integer, IScan>> ms2scans = ms2scansIndex.getNum2scan().entrySet();
//
//				for (final Map.Entry<Integer, IScan> kv : ms2scansIndex.getNum2scan().entrySet()) {
//					exec.submit(() -> {
//						IScan scan = kv.getValue();
//
//						// assign this spectrum to it's PSM
//						// Todo: This is really slow Better to have a hashMap
//						for (PSM p : psmList) {
//							if (p.srcFile.equalsIgnoreCase(fileName) && (p.scanNum == kv.getKey())) {
//								final ISpectrum spectrum = scan.getSpectrum();
//								int N = spectrum.getMZs().length;
//								if(N == 0) {
//									continue; // no valid spectrum for this scan number
//								}
//
//								double[] mz = spectrum.getMZs();
//								double[] intensities = spectrum.getIntensities();
//
//								// If this happens, there is something wrong with the spectrum so skip it
//								if(mz.length != intensities.length) {
//									System.err.print(
//											"\nERROR:" + fileName + " Scan: " + kv.getValue() +
//													"\n# of mz values != # intensity values: " +
//													mz.length + " != " + intensities.length +
//													"\nSkipping this scan...\n"
//									);
//									continue;
//								}
//
//								SpectrumClass X = new SpectrumClass(mz, intensities);
//								p.recordSpectra(X);
//								ctr.getAndIncrement();
//								break;
//							}
//						}
//						System.err.print("\r" + fileName +  ":  " + ctr + " spectra read in.            ");
//					});
//				}
//
//			} catch (FileParsingException e) {
//				e.printStackTrace();
//			}
//		});
//
//	}


    static double getFragmentIonMass(String x, double z, double addl_mass) {
		double ret = 1.00728 * z;
		
		int start = x.indexOf(":") + 1;
		int stop  = x.length();
		
		if(x.contains("-")) stop = x.indexOf("-");
		
		for(int i = start; i < stop; i++) {
			String c = Character.toString( x.charAt(i) );
			
			if( AA_MASS_MAP.containsKey(c) ) {
				double mass = AA_MASS_MAP.get(c);
				ret += mass;
			}
		}
		
		ret += addl_mass; // y-ions have a water molecule extra
		
		return ret;
	}

	
	// Function reads in an MGF file and returns it as a HashMap
	private static TIntObjectHashMap<SpectrumClass > read_mgf(String specFile) throws IOException {
		TIntObjectHashMap<SpectrumClass > ret = new TIntObjectHashMap<>();
		
		File mgf = new File(specFile);
		BufferedReader br = new BufferedReader(new FileReader(mgf));
		String line;
		int scanNum = 0;
        SpectrumClass S = null;
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
                    S = new SpectrumClass(mz, I);
                    ret.put(scanNum, S);

                    mzAL = null;
                    intensityAL = null;
                    mz = null;
                    I = null;
				}
				S = null;
				scanNum = 0;
			}
					
			if(line.startsWith("BEGIN IONS")) {
                mzAL = new ArrayList<>();
                intensityAL = new ArrayList<>();
			}
			
			if(line.startsWith("TITLE=")) {
				int i = line.indexOf(".") + 1;
				int j = line.indexOf(".", i);
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

	/**
	 * Function to read spectra from mzXML file
	 * @param scanMap Scan Map
	 * @throws IllegalStateException
	 * @throws FileParsingException
	 */
	private static void readMzXML(Map<String, List<Integer>> scanMap, String pathSpectra) throws
			IllegalStateException,
			FileParsingException {
		
		// Iterate over the file names
		for(String fn : scanMap.keySet()) {
			String baseFN = new File(pathSpectra + "/" + fn).getName();
			System.err.print(baseFN + ":  "); // beginning of info line
			
			int ctr = 0;
			List<Integer> scanNums = scanMap.get(fn);
			Collections.sort(scanNums); // order the scan numbers

            int N = numThreads;
            if(numThreads > 1) N -= 1;

			final MZXMLFile mzxml = new MZXMLFile(fn, false);
			mzxml.setNumThreadsForParsing(N);
            mzxml.setParsingTimeout(60L); // 1 minute before it times out trying to read a file
			final LCMSData lcmsData = new LCMSData(mzxml);
			lcmsData.load(LCMSDataSubset.MS2_WITH_SPECTRA);
			final IScanCollection scans = lcmsData.getScans();
			final ScanIndex ms2ScanIndex = scans.getMapMsLevel2index().get(2);

			if( (ms2ScanIndex == null) || (ms2ScanIndex.getNum2scan().isEmpty()) ) {
				log.info("\nERROR: Globals.readMzXML(): Unable to read MS2 scans from '" + fn + "'\n");
				System.exit(0);
			}
			
			for(Map.Entry<Integer, IScan> num2scan : ms2ScanIndex.getNum2scan().entrySet()) {
				int scanNum = num2scan.getKey();
				IScan scan = num2scan.getValue();
				double[] mz = scan.getSpectrum().getMZs();
				double[] intensities = scan.getSpectrum().getIntensities();

                SpectrumClass curSpectrum = new SpectrumClass(mz, intensities);

				// assign this spectrum to it's PSM
				for(PSM p : psmList) {
					if( (p.getSrcFile().equalsIgnoreCase(baseFN)) && (p.getScanNum() == scanNum) ) {
						p.recordSpectra(curSpectrum);
						ctr++;
						break;
					}
				}
			}
			
			log.info(ctr + " spectra read in.");  // end of info line
		}
	}

	
	// Function returns the decoy 'key' from the decoy map for the given residue
	static String getDecoySymbol(char c) {
		String ret = "";
		String srcChar = Character.toString(c);
		
		for(String k : DECOY_AA_MAP.keySet()) {
			String v = DECOY_AA_MAP.get(k);
			
			if(v.equalsIgnoreCase(srcChar)) {
				ret = k;
				break;
			}
		}
		
		return ret;
	}

	
	// Function returns the TPP-formatted representation of the given single
	// character modification
	static String getTPPresidue(String c) {
		String ret = "";
		String orig = "";

        if(c.equalsIgnoreCase("[")) {
            int d = (int) Math.round(Globals.ntermMass) + 1; // adds a proton
            ret = "n[" + d + "]";
        }
        else if(c.equals("]")) {
            int d = (int) Math.round(Globals.ctermMass);
            ret += "c[" + d + "]";
        }
        else {
            int i = (int) Math.round(AA_MASS_MAP.get(c));

            if (isDecoyResidue(c)) {
                orig = DECOY_AA_MAP.get(c);
            } else orig = c.toUpperCase();

            ret = orig + "[" + i + "]";
        }
		
		return ret;
	}

	
	// Function returns true if the given sequence contains decoy characters
	static boolean isDecoySeq(String seq) {
		boolean ret = false;
		
		int score = 0;
		for(int i = 0; i < seq.length(); i++) {
			String c = Character.toString( seq.charAt(i) );
			if(c.equalsIgnoreCase("[")) continue; // n-term mod symbol
			if(c.equalsIgnoreCase("]")) continue; // c-term mod symbol
			
			if(isDecoyResidue(c)) score++;
		}
		
		if(score > 0) ret = true;
		
		return ret;
	}
	
	
	// Function returns true if the given residue is from the decoy list
	static boolean isDecoyResidue(String AA) {
		boolean ret = false;
		if(DECOY_AA_MAP.containsKey(AA)) ret = true;
		return ret;
	}
	
	
	// Function writes a sample input file for LucXor to disk
	static File writeTemplateInputFile() throws IOException {
		File outF = new File("luciphor2_input_template.txt");
		
		FileWriter fw = new FileWriter(outF.getAbsoluteFile());
		BufferedWriter bw = new BufferedWriter(fw);
		
		bw.write("## Input file for Luciphor2 (aka: LucXor).\n## Anything after a hash '#' is ignored\n");
		bw.write("## By default, these initial parameters are for performing a phosphorylation search\n\n");
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

        bw.write("WRITE_MATCHED_PEAKS_FILE = 0 ## Generate a tab-delimited file of all the matched peaks\n" +
                 "                             ## for the top 2 predictions of each spectra\n" +
                 "                             ## Useful for plotting spectra\n" +
                 "                             ## 0 = no, 1 = yes\n\n");

		bw.write("## Place here any FIXED modifications that were used in your search\n" +
				 "## This field is ONLY used for tab-delimited input\n" +
				 "## Syntax: FIXED_MOD = <RESIDUE> <MODIFICATION_MASS>\n" +
				 "## For N-terminal modifications use '[' for the residue character\n" +
				 "## For C-terminal modifications use ']' for the residue character\n" +
				 "FIXED_MOD = C 57.021464\n\n");
		
		bw.write("## Place here any VARIABLE modifications that might be encountered that you don't\n" +
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
				 "#NL = STE -H2O -18.01056    ## a correctly formatted example, (not actually recommended for phospho-searches)\n" +
				 "#NL = RKQN -NH3 -17.026548  ## another correctly formatted example, (again not recommended for phospho-searches)\n" +
				 "NL = sty -H3PO4 -97.97690\n\n");
		
		bw.write("DECOY_MASS = 79.966331  ## how much to add to an amino acid to make it a decoy\n\n");
		bw.write("## For handling the neutral loss from a decoy sequence.\n" +
                 "## The syntax for this is identical to that of the normal neutral losses given\n" +
                 "## above except that the residue is always 'X'\n" +
                 "## Syntax: DECOY_NL = X -<NEUTRAL_LOSS_MOLECULAR_FORMULA> <MASS_LOST>\n" +
                 "DECOY_NL = X -H3PO4 -97.97690\n\n");


		bw.write("MAX_CHARGE_STATE = 5 ## do not consider PSMs with a charge state above this value\n");
		bw.write("MAX_PEP_LEN = 40 ## restrict scoring to peptides with a length shorter than this value\n");
        bw.write("MAX_NUM_PERM = 16384 ## the maximum number of permutations a sequence can have\n\n");

        bw.write("SELECTION_METHOD = 0   ## 0 = Peptide Prophet probability (default)\n" +
				 "                       ## 1 = Mascot Ion Score\n" +
				 "                       ## 2 = -log(E-value)\n" +
				 "                       ## 3 = X!Tandem Hyperscore\n" +
                 "                       ## 4 = Sequest Xcorr\n\n");

		bw.write("MODELING_SCORE_THRESHOLD = 0.95 ## Minimum score a PSM needs to be considered for modeling\n");
		bw.write("                                ## The default assumes you are using SELECTION_METHOD=0\n");
		bw.write("                                ## If using SELECTION_METHOD=2 then set this value to your desired e-value it will be converted to -log(e-value) internally\n");
		bw.write("SCORING_THRESHOLD = 0    ## PSMs below this value will be discarded\n");
		bw.write("                         ## Again, if using SELECTION_METHOD=2 then set this value to your desired e-value it will be converted to -log(e-value) internally\n");
		bw.write("MIN_NUM_PSMS_MODEL = 50  ## The minimum number of PSMs you need for any charge state in order to build a model for it\n\n");
		

		bw.write("MOD_PEP_REP = 0 ## 0 = show single character modifications, 1 = show TPP-formatted modifications\n\n");
		
		bw.write("NUM_THREADS = 0 ## For multi-threading, zero = use all CPU found by JAVA\n\n");

        bw.write("RUN_MODE = 0 ## Determines how Luciphor will run.\n" +
                 "             ## 0 = Default: calculate FLR then rerun scoring without decoys (two iterations)\n" +
                 "             ## 1 = Report Decoys: calculate FLR but don't rescore PSMs, all decoy hits will be reported\n\n");

		bw.write("## This option can be used to help diagnose problems with Luciphor. Multi-threading is disabled in debug mode.\n" +
                "DEBUG_MODE = 0 ## 0 = default: turn off debugging\n" +
                 "               ## 1 = write peaks selected for modeling to disk\n" +
                 "               ## 2 = write the scores of all permutations for each PSM to disk\n" +
				 "               ## 3 = write the matched peaks for the top scoring permutation to disk\n" +
				 "               ## 4 = write HCD non-parametric models to disk (HCD-mode only option)\n\n");
		
		bw.close();
		return outF;
    }


    // Record the global and local FLR values estimated for all of the delta scores
    public static void recordFLRestimates() {
        FLRestimateMap = new THashMap<>();

        for(PSM p : Globals.psmList) {
            if(p.isDecoy()) continue; // skip FLR data from decoys
            double[] d = new double[2];
            d[0] = p.getGlobalFDR();
            d[1] = p.getLocalFDR();
            FLRestimateMap.put( p.getDeltaScore(), d );
        }
    }


    // Function assigns the global and local FLR for the current PSM from the FLRestimateMap
    public static void assignFLR() {

        ArrayList<Double> obsDeltaScores = new ArrayList<>(FLRestimateMap.keySet());
        Collections.sort(obsDeltaScores);  // sort them from low to high
        int N = obsDeltaScores.size();
        boolean assigned;
        for(PSM p : Globals.psmList) {
            double obs_ds = p.getDeltaScore();
            assigned = false;

            // iterate over the delta scores until you find the value closest to this one
            int i = 0;
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


    // This function prepares each PSM for the second iteration (the one after the FLR has been estimated)
    public static void clearPSMs() {
        for(PSM p : psmList) p.clearScores();
    }
}
