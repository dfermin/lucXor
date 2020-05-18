/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import gnu.trove.map.TMap;
import gnu.trove.map.hash.THashMap;
import lombok.extern.slf4j.Slf4j;
import lucxor.utils.Constants;
import lucxor.utils.Utils;

import java.io.*;
import java.util.*;

import static lucxor.utils.Constants.*;

/**
 *
 * @author dfermin
 */
@Slf4j
public class LucXorConfiguration {

	private static File SPECTRUM_PATH = null;
	private static String SPECTRUM_PREFIX = null;
    private static String matchedPkFile = null;
	private static File inputFile = null;
	private static String outputFile = null;
	private static int MS2TOL_UNITS;
	private static int inputType;
	private static int debugMode;
	private static int reduceNL;
	private static int peptideRepresentation;
	private static int scoringMethod;
	private static int scoringAlgorithm;
	private static int maxChargeState;
	private static int minNumPSMsForModeling;
	private static int numThreads;
    private static int runMode;
    private static int tsvHdr;
    private static int maxPepLen;
	private static double ms2tol;
	private static double modelTH;
	private static double scoreTH;
	private static double DECOY_MASS;
	private static double minMZ;
	private static double max_num_permutations;
	private static double precursorNLmass;
	private static double ntermMass = 0d;
	private static double ctermMass = 0d;
    private static boolean writeMatchedPeaks;

	// mods user wants to search for
	private static final TMap<String, Double> TARGET_MOD_MAP = new THashMap<>();
	// fixed mods observed in data
	private static final TMap<String, Double> FIXED_MOD_MAP = new THashMap<>();
	// variable mods observed in data
	private static final TMap<String, Double> VAR_MOD_MAP = new THashMap<>();
	// holds all neutral loss masses, k= list of amino acids v= NL mass
	private static final TMap<String, Double> NEUTRAL_LOSS_MAP = new THashMap<>();
	private static final TMap<String, Double> DECOY_NEUTRAL_LOSS_MAP = new THashMap<>();

	public static void parseConfigurationFile(String str) throws IOException {
		
		File inF = new File(str);
		if(!inF.exists()) { 
			System.err.print("\nERROR! Unable to open " + inF.getAbsolutePath() + "\n\n");
			System.exit(0);
		}
		
		debugMode = 0; // 0 means no debugging output
        runMode = 0; // calculate FLR and rescore PSMs without decoys (2 iterations)
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
				SPECTRUM_PATH = new File(s).getCanonicalFile();
			}
			
			if(line.startsWith("SPECTRUM_SUFFIX")) {
				String s = Utils.parseInputLine(line);
				SPECTRUM_PREFIX = s.toLowerCase();
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
				MS2TOL_UNITS = Integer.valueOf(s);
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
                else if(x == 0 || (numThreads > Runtime.getRuntime().availableProcessors()))
                	numThreads = Runtime.getRuntime().availableProcessors() - 1;
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
				DECOY_MASS = Double.valueOf(s);
			}

            if(line.startsWith("DECOY_NL")) {
                String[] ary = parseNLLine(line);
                String k = ary[0].substring(1);
                double m = Double.valueOf(ary[1]);
                DECOY_NEUTRAL_LOSS_MAP.put(k,m);
            }
			
			if(line.startsWith("NL")) {
				String[] ary = parseNLLine(line);
				double m = Double.valueOf(ary[1]);
				NEUTRAL_LOSS_MAP.put(ary[0], m);
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
				String[] ary = parseInputModLine(line);
				double m = Double.valueOf(ary[1]);
				TARGET_MOD_MAP.put(ary[0].toUpperCase(), m);
			}


            if(line.startsWith("VAR_MOD")) {
                String[] ary = parseInputModLine(line);
                double m = Double.valueOf(ary[1]);
                VAR_MOD_MAP.put(ary[0].toLowerCase(), m); // variable mods are lower case
            }

            if(line.startsWith("FIXED_MOD")) {
                String[] ary = parseInputModLine(line);
                double m = Double.valueOf(ary[1]);
                FIXED_MOD_MAP.put(ary[0].toUpperCase(), m);
            }


			// You only need to extract the VAR_MOD and FIXED_MOD values
			// from the input file if you are using TSV files
//			if(LucXorConfiguration.inputType == Constants.TSV) {
//				if(line.startsWith("VAR_MOD")) {
//					String[] ary = parseInputModLine(line);
//					double m = Double.valueOf(ary[1]);
//					VAR_MOD_MAP.put(ary[0].toLowerCase(), m); // variable mods are lower case
//				}
//
//				if(line.startsWith("FIXED_MOD")) {
//					String[] ary = parseInputModLine(line);
//					double m = Double.valueOf(ary[1]);
//					FIXED_MOD_MAP.put(ary[0].toUpperCase(), m);
//				}
//			}

		}
		br.close();
		
		if( (null == outputFile) || (outputFile.isEmpty()) )
			outputFile = "luciphor_results.tsv";
		
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

		log.info("Spectrum Path:           " + SPECTRUM_PATH.getAbsolutePath());
		log.info("Spectrum Suffix:         " + SPECTRUM_PREFIX);
		log.info("Input file:              " + inputFile);
		log.info("Input type:              " + (inputType == Constants.PEPXML ? "pepXML" : "tsv"));
		log.info("MS2 tolerance:           " + ms2tol + (MS2TOL_UNITS == Constants.DALTONS ? " Da" : " ppm"));
		log.info("Luciphor Algorithm:      " + (scoringAlgorithm == Constants.CID ? "CID" : "HCD") );
		log.info("Classifying on:          " + classStr);
        log.info("Run Mode:                " + (runMode == 0 ? "Default" : "Report Decoys"));
		log.info("Num of Threads:          " + NCPU );
		log.info("Modeling Threshold:      " + modelTH);
		log.info("Scoring Threshold:       " + scoreTH);
		log.info("Permutation Limit:       " + max_num_permutations);
        log.info("Max peptide length:      " + maxPepLen);
		log.info("Min num PSMs for model:  " + minNumPSMsForModeling);
		log.info("Decoy Mass Adduct:       " + DECOY_MASS);
		log.info("Max Charge State:        " + maxChargeState);
		log.info("Reduce NL:               " + (reduceNL == 0 ? "no" : "yes"));
		log.info("Output File:             " + outputFile);
        log.info("Write matched Peaks:     " + (writeMatchedPeaks ? "yes" : "no"));
        System.err.print("\n");

		if(debugMode != 0) {
			log.info("Debug mode:              " + debugMode + "  (Limiting to 1 CPU)\n");
		    LucXorConfiguration.numThreads = 1;
        }
		
		log.info("Mods to score:");
		for(Map.Entry<String, Double> stringDoubleEntry : TARGET_MOD_MAP.entrySet()) {
			log.info(stringDoubleEntry.getKey() + "\t" + stringDoubleEntry.getValue());
		}
		
		if(!NEUTRAL_LOSS_MAP.isEmpty()) {
			log.info("\nAllowed Neutral Losses:");
			for(Map.Entry<String, Double> stringDoubleEntry : NEUTRAL_LOSS_MAP.entrySet()) {
				log.info(stringDoubleEntry.getKey() + "\t" + stringDoubleEntry.getValue());
			}
            for(Map.Entry<String, Double> stringDoubleEntry : DECOY_NEUTRAL_LOSS_MAP.entrySet()) {
                log.info("<X>" + stringDoubleEntry.getKey() + "\t" + stringDoubleEntry.getValue() + "  (Decoy NL)");
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

	/**
	 * Parse a configuration line with the sutructre key = value
	 * @param line Configuration line
	 * @return values
	 */
	private static String[] parseInputModLine(String line) {
		String[] ret = new String[2];
		
		char aa = 0;
		StringBuilder sb = new StringBuilder(30);
		double mass;
		int N = line.length();
		
		int b = line.indexOf('=') + 1;
		
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

	/**
	 * Parse a line with the structure key=value
	 * @param line Configuration line
	 * @return List of values
	 */
	private static String[] parseNLLine(String line) {
		String[] ret = new String[2];
		
		line = line.replaceAll("#", "");
		
		String[] tmp = line.split("\\s+");
				
		ret[0] = tmp[2] + tmp[3];
		ret[1] = tmp[4];
		
		return ret;
	}

	/**
	 * Function to incorporate information read in from input file into the
	 * AA_MASS_MAP object
	 */
	public static void loadUserMods() {

		// Assign the information you got from the input file to the 
		// relevant map. If none were given, then this should just proceed
		// without changing anything.
		for(Map.Entry<String, Double> stringDoubleEntry : TARGET_MOD_MAP.entrySet()) {
			double mass = AA_MASS_MAP.get(stringDoubleEntry.getKey()) + stringDoubleEntry.getValue();
			String symbol = stringDoubleEntry.getKey().toLowerCase();
			AA_MASS_MAP.put(symbol, mass);
		}
		
		for(Map.Entry<String, Double> stringDoubleEntry : FIXED_MOD_MAP.entrySet()) {

		    if( stringDoubleEntry.getKey().equalsIgnoreCase("[") ) {
                LucXorConfiguration.ntermMass = stringDoubleEntry.getValue();
            }
            else if( stringDoubleEntry.getKey().equalsIgnoreCase("]") ) {
                LucXorConfiguration.ctermMass = stringDoubleEntry.getValue();
            }
		    else {
                double mass = AA_MASS_MAP.get(stringDoubleEntry.getKey()) + stringDoubleEntry.getValue();
                String symbol = stringDoubleEntry.getKey().toUpperCase();
                AA_MASS_MAP.put(symbol, mass);
            }
		}
		
		for(Map.Entry<String, Double> stringDoubleEntry : VAR_MOD_MAP.entrySet()) { // this will be a lower case key
			
			if( stringDoubleEntry.getKey().equalsIgnoreCase("[") ) {
				LucXorConfiguration.ntermMass = stringDoubleEntry.getValue();
			}
			else if( stringDoubleEntry.getKey().equalsIgnoreCase("]") ) {
				LucXorConfiguration.ctermMass = stringDoubleEntry.getValue();
			}
			else {
				double mass = AA_MASS_MAP.get(stringDoubleEntry.getKey().toUpperCase()) + stringDoubleEntry.getValue();
				AA_MASS_MAP.put(stringDoubleEntry.getKey(), mass);
			}
		}
		
		
		// Now add the decoy masses to the AA_MASS_MAP
		for(Map.Entry<String, String> stringStringEntry : DECOY_AA_MAP.entrySet()) {
			String trueAA = stringStringEntry.getValue();
			
			if(VAR_MOD_MAP.containsKey(trueAA)) continue;
			if(TARGET_MOD_MAP.containsKey(trueAA)) continue;
			
			double mass = AA_MASS_MAP.get(trueAA) + LucXorConfiguration.DECOY_MASS;
			AA_MASS_MAP.put(stringStringEntry.getKey(), mass);
		}

	}

//	/****************
//	 * Function reads in spectral data from mzML files
//	 * @param scanMap
//	 */
//	private static void readXMLFile(Map<String, List<Integer>> scanMap,
//									String SPECTRUM_PATH, String fileType) throws FileParsingException {
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
//				source = new MZMLFile(SPECTRUM_PATH + "/" + fileName);
//
//			if(fileType.equalsIgnoreCase(MZXML_TYPE))
//				source = new MZXMLFile(SPECTRUM_PATH + "/" + fileName);
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
//								Spectrum X = new Spectrum(mz, intensities);
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

	public static File getSpectrumPath() {
		return SPECTRUM_PATH;
	}

	public static String getSpectrumPrefix() {
		return SPECTRUM_PREFIX;
	}

	public static String getMatchedPkFile() {
		return matchedPkFile;
	}

	public static File getInputFile() {
		return inputFile;
	}

	public static String getOutputFile() {
		return outputFile;
	}

	public static int getMs2tolUnits() {
		return MS2TOL_UNITS;
	}

	public static int getInputType() {
		return inputType;
	}

	public static int getDebugMode() {
		return debugMode;
	}

	public static int getReduceNL() {
		return reduceNL;
	}

	public static int getPeptideRepresentation() {
		return peptideRepresentation;
	}

	public static int getScoringMethod() {
		return scoringMethod;
	}

	public static int getScoringAlgorithm() {
		return scoringAlgorithm;
	}

	public static int getMaxChargeState() {
		return maxChargeState;
	}

	public static int getMinNumPSMsForModeling() {
		return minNumPSMsForModeling;
	}

	public static int getNumThreads() {
		return numThreads;
	}

	public static int getRunMode() {
		return runMode;
	}

	public static int getTsvHdr() {
		return tsvHdr;
	}

	public static int getMaxPepLen() {
		return maxPepLen;
	}

	public static double getMs2tol() {
		return ms2tol;
	}

	public static double getModelTH() {
		return modelTH;
	}

	public static double getScoreTH() {
		return scoreTH;
	}

	public static double getMinMZ() {
		return minMZ;
	}

	public static double getMax_num_permutations() {
		return max_num_permutations;
	}

	public static double getPrecursorNLmass() {
		return precursorNLmass;
	}

	public static double getNtermMass() {
		return ntermMass;
	}

	public static double getCtermMass() {
		return ctermMass;
	}

	public static boolean isWriteMatchedPeaks() {
		return writeMatchedPeaks;
	}

	public static TMap<String, Double> getTargetModMap() {
		return TARGET_MOD_MAP;
	}

	public static TMap<String, Double> getFixedModMap() {
		return FIXED_MOD_MAP;
	}

	public static TMap<String, Double> getVarModMap() {
		return VAR_MOD_MAP;
	}

	public static TMap<String, Double> getNeutralLossMap() {
		return NEUTRAL_LOSS_MAP;
	}

	public static TMap<String, Double> getDecoyNeutralLossMap() {
		return DECOY_NEUTRAL_LOSS_MAP;
	}

	public static void setNTermMass(double modMass) {
		ntermMass = modMass;
	}

	public static void setCTermMass(double modMass) {
		ctermMass = modMass;
	}

	public static void clearVarMods() {
		VAR_MOD_MAP.clear();
	}

	public static void addVarMod(String key, double mass) {
		VAR_MOD_MAP.put(key, mass);
	}
}
