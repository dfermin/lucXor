/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import gnu.trove.map.TMap;
import gnu.trove.map.hash.THashMap;
import lucxor.utils.Constants;
import lucxor.utils.Utils;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;

import static lucxor.utils.Constants.*;

/**
 *
 * @author dfermin
 */

public class LucXorConfiguration {

	private static final org.slf4j.Logger log = LoggerFactory
			.getLogger(LucXorConfiguration.class);

	private static File SPECTRUM_PATH = null;
	private static String SPECTRUM_PREFIX = null;
    private static String MATCHED_PKFILE = null;
	private static File INPUT_FILE = null;
	private static String OUTPUT_FILE = null;
	private static int MS2TOL_UNITS;
	private static int INPUT_TYPE;
	private static int DEBUG_MODE;
	private static int REDUCENL;
	private static int PEPTIDE_REPRESENTATION;
	private static int SCORING_METHOD;
	private static int SCORING_ALGORITHM;
	private static int MAX_CHARGE_STATE;
	private static int MIN_NUM_PSM_MODELING;
	private static int NUM_THREADS;
    private static int RUN_MODE;
    private static int TSV_HDR;
    private static int MAX_PEP_LENGTH;
	private static double MS2TOL;
	private static double MODELTH;
	private static double SCORETH;
	private static double DECOY_MASS;
	private static double MIN_MZ;
	private static double MAX_NUM_PERMUTATIONS;
	private static double PRECURSOR_NL_MASS;
	private static double NTERM_MASS = 0d;
	private static double CTERM_MASS = 0d;
    private static boolean WRITE_MATCHED_PEAKS;

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
		
		DEBUG_MODE = 0; // 0 means no debugging output
        RUN_MODE = 0; // calculate FLR and rescore PSMs without decoys (2 iterations)
        MIN_NUM_PSM_MODELING = 50;
        MAX_PEP_LENGTH = 40;
        REDUCENL = 0;
		NUM_THREADS = Runtime.getRuntime().availableProcessors();
		TSV_HDR = 0; // default is no header
        WRITE_MATCHED_PEAKS = false; // true means you will generate the matched peaks

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
				INPUT_FILE = new File(s);
			}
			
			if(line.startsWith("OUTPUT_FILE")) {
				OUTPUT_FILE = Utils.parseInputLine(line);
			}
			
			if(line.startsWith("INPUT_TYPE")) {
				String s = Utils.parseInputLine(line);
				INPUT_TYPE = Integer.valueOf(s);
			}

            if(line.startsWith("MAX_PEP_LEN")) {
                String s = Utils.parseInputLine(line);
                MAX_PEP_LENGTH = Integer.valueOf(s);
            }

			if(line.startsWith("MAX_NUM_PERM")) {
				String s = Utils.parseInputLine(line);
				MAX_NUM_PERMUTATIONS = Double.valueOf(s);
			}
			
			if(line.startsWith("MIN_NUM_PSMS_MODEL")) {
				String s = Utils.parseInputLine(line);
				MIN_NUM_PSM_MODELING = Integer.valueOf(s);
			}
			
			if(line.startsWith("MS2_TOL") && !line.contains("_UNITS")) {
				String s = Utils.parseInputLine(line);
				MS2TOL = Double.valueOf(s);
			}
			
			if(line.startsWith("MS2_TOL_UNITS")) {
				String s = Utils.parseInputLine(line);
				MS2TOL_UNITS = Integer.valueOf(s);
			}
			
			if(line.startsWith("ALGORITHM")) {
				String s = Utils.parseInputLine(line);
				SCORING_ALGORITHM = Integer.valueOf(s);
			}

            if(line.startsWith("TSV_HEADER")) {
                String s = Utils.parseInputLine(line);
                TSV_HDR = Integer.valueOf(s);
            }

			if(line.startsWith("REDUCE_PRECURSOR_NL")) {
				String s = Utils.parseInputLine(line);
				REDUCENL = Integer.valueOf(s);
			}
			
			if(line.startsWith("PRECURSOR_NL_MASS_DIFF")) {
				String s = Utils.parseInputLine(line);
				PRECURSOR_NL_MASS = Double.valueOf(s);
			}
			
			if(line.startsWith("SELECTION_METHOD")) {
				String s = Utils.parseInputLine(line);
				SCORING_METHOD = Integer.valueOf(s);
			}
			
			if(line.startsWith("MODELING_SCORE_THRESHOLD")) {
				String s = Utils.parseInputLine(line);
				MODELTH = Double.valueOf(s);
			}
			
			if(line.startsWith("MAX_CHARGE_STATE")) {
				String s = Utils.parseInputLine(line);
				MAX_CHARGE_STATE = Integer.valueOf(s);
			}
			
			if(line.startsWith("NUM_THREADS")) {
				String s = Utils.parseInputLine(line);
				int x = Integer.valueOf(s);
				// We do minus 1 because 1 thread already goes to 
				// running the whole program
				if(x < 0) NUM_THREADS = 1;
                else if(x == 0 || (NUM_THREADS > Runtime.getRuntime().availableProcessors()))
                	NUM_THREADS = Runtime.getRuntime().availableProcessors() - 1;
                else NUM_THREADS = x;
			}
			
			if(line.startsWith("DEBUG_MODE")) {
				String s = Utils.parseInputLine(line);
				DEBUG_MODE = Integer.valueOf(s);
			}

            if(line.startsWith("WRITE_MATCHED_PEAKS_FILE")) {
                String s = Utils.parseInputLine(line);
                if(s.equals("1")) WRITE_MATCHED_PEAKS = true;
            }

            if(line.startsWith("RUN_MODE")) {
                String s = Utils.parseInputLine(line);
                RUN_MODE = Integer.valueOf(s);
                if(RUN_MODE > 1) RUN_MODE = 0;
            }
			
			if(line.startsWith("SCORING_THRESHOLD")) {
				String s = Utils.parseInputLine(line);
				SCORETH = Double.valueOf(s);
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
				MIN_MZ = Double.valueOf(s);
			}

			if(line.startsWith("MOD_PEP_REP")) {
				String s = Utils.parseInputLine(line);
				PEPTIDE_REPRESENTATION = Integer.valueOf(s);
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
//			if(LucXorConfiguration.INPUT_TYPE == Constants.TSV) {
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
		
		if( (null == OUTPUT_FILE) || (OUTPUT_FILE.isEmpty()) )
			OUTPUT_FILE = "luciphor_results.tsv";
		
		String classStr = "";
		switch(SCORING_METHOD) {
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
				System.err.print("\nERROR! Unknown scoring method: " + SCORING_METHOD + "\n\n");
				System.exit(0);
				break;
		}

        int NCPU = NUM_THREADS;
        if( (NUM_THREADS > 1) && (NUM_THREADS < Runtime.getRuntime().availableProcessors()) ) NCPU = (NUM_THREADS + 1);

		log.info("Spectrum Path:           " + SPECTRUM_PATH.getAbsolutePath());
		log.info("Spectrum Suffix:         " + SPECTRUM_PREFIX);
		log.info("Input file:              " + INPUT_FILE);
		log.info("Input type:              " + (INPUT_TYPE == Constants.PEPXML ? "pepXML" : "tsv"));
		log.info("MS2 tolerance:           " + MS2TOL + (MS2TOL_UNITS == Constants.DALTONS ? " Da" : " ppm"));
		log.info("Luciphor Algorithm:      " + (SCORING_ALGORITHM == Constants.CID ? "CID" : "HCD") );
		log.info("Classifying on:          " + classStr);
        log.info("Run Mode:                " + (RUN_MODE == 0 ? "Default" : "Report Decoys"));
		log.info("Num of Threads:          " + NCPU );
		log.info("Modeling Threshold:      " + MODELTH);
		log.info("Scoring Threshold:       " + SCORETH);
		log.info("Permutation Limit:       " + MAX_NUM_PERMUTATIONS);
        log.info("Max peptide length:      " + MAX_PEP_LENGTH);
		log.info("Min num PSMs for model:  " + MIN_NUM_PSM_MODELING);
		log.info("Decoy Mass Adduct:       " + DECOY_MASS);
		log.info("Max Charge State:        " + MAX_CHARGE_STATE);
		log.info("Reduce NL:               " + (REDUCENL == 0 ? "no" : "yes"));
		log.info("Output File:             " + OUTPUT_FILE);
        log.info("Write matched Peaks:     " + (WRITE_MATCHED_PEAKS ? "yes" : "no"));
        System.err.print("\n");

		if(DEBUG_MODE != 0) {
			log.info("Debug mode:              " + DEBUG_MODE + "  (Limiting to 1 CPU)\n");
		    LucXorConfiguration.NUM_THREADS = 1;
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
		if(SCORING_METHOD == Constants.NEGLOGEXPECT) {
			double tmp = -1.0 * Math.log(MODELTH);
			MODELTH = tmp;

			tmp = -1.0 * Math.log(SCORETH);
			SCORETH = tmp;
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
                LucXorConfiguration.NTERM_MASS = stringDoubleEntry.getValue();
            }
            else if( stringDoubleEntry.getKey().equalsIgnoreCase("]") ) {
                LucXorConfiguration.CTERM_MASS = stringDoubleEntry.getValue();
            }
		    else {
                double mass = AA_MASS_MAP.get(stringDoubleEntry.getKey()) + stringDoubleEntry.getValue();
                String symbol = stringDoubleEntry.getKey().toUpperCase();
                AA_MASS_MAP.put(symbol, mass);
            }
		}
		
		for(Map.Entry<String, Double> stringDoubleEntry : VAR_MOD_MAP.entrySet()) { // this will be a lower case key
			
			if( stringDoubleEntry.getKey().equalsIgnoreCase("[") ) {
				LucXorConfiguration.NTERM_MASS = stringDoubleEntry.getValue();
			}
			else if( stringDoubleEntry.getKey().equalsIgnoreCase("]") ) {
				LucXorConfiguration.CTERM_MASS = stringDoubleEntry.getValue();
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

	public static String getMatchedPkfile() {
		return MATCHED_PKFILE;
	}

	public static File getInputFile() {
		return INPUT_FILE;
	}

	public static String getOutputFile() {
		return OUTPUT_FILE;
	}

	public static int getMs2tolUnits() {
		return MS2TOL_UNITS;
	}

	public static int getInputType() {
		return INPUT_TYPE;
	}

	public static int getDebugMode() {
		return DEBUG_MODE;
	}

	public static int getREDUCENL() {
		return REDUCENL;
	}

	public static int getPeptideRepresentation() {
		return PEPTIDE_REPRESENTATION;
	}

	public static int getScoringMethod() {
		return SCORING_METHOD;
	}

	public static int getScoringAlgorithm() {
		return SCORING_ALGORITHM;
	}

	public static int getMaxChargeState() {
		return MAX_CHARGE_STATE;
	}

	public static int getMinNumPsmModeling() {
		return MIN_NUM_PSM_MODELING;
	}

	public static int getNumThreads() {
		return NUM_THREADS;
	}

	public static int getRunMode() {
		return RUN_MODE;
	}

	public static int getTsvHdr() {
		return TSV_HDR;
	}

	public static int getMaxPepLength() {
		return MAX_PEP_LENGTH;
	}

	public static double getMS2TOL() {
		return MS2TOL;
	}

	public static double getMODELTH() {
		return MODELTH;
	}

	public static double getSCORETH() {
		return SCORETH;
	}

	public static double getMinMz() {
		return MIN_MZ;
	}

	public static double getMaxNumPermutations() {
		return MAX_NUM_PERMUTATIONS;
	}

	public static double getPrecursorNlMass() {
		return PRECURSOR_NL_MASS;
	}

	public static double getNtermMass() {
		return NTERM_MASS;
	}

	public static double getCtermMass() {
		return CTERM_MASS;
	}

	public static boolean isWriteMatchedPeaks() {
		return WRITE_MATCHED_PEAKS;
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
		NTERM_MASS = modMass;
	}

	public static void setCTermMass(double modMass) {
		CTERM_MASS = modMass;
	}

	public static void clearVarMods() {
		VAR_MOD_MAP.clear();
	}

	public static void addVarMod(String key, double mass) {
		VAR_MOD_MAP.put(key, mass);
	}
}
