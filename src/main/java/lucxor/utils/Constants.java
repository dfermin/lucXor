/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor.utils;

import java.util.AbstractMap;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * This class contains a set of contants that are needed for the program.
 * @author dfermin
 * @author ypriverol
 */
public class Constants {

	public static final int PEPPROPHET = 0;
	public static final int MASCOTIONSCORE = 1;
	public static final int NEGLOGEXPECT = 2;
    public static final int XTDHYPERSCORE = 3;
    public static final int XCORR = 4;
	public static final int DALTONS = 0;
	public static final int PPM_UNITS = 1;
	public static final int PEPXML = 0;
	public static final int NTERM_MOD = -100;
	public static final int CTERM_MOD = 100;
	public static final int CID = 0;
	public static final int HCD = 1;
	public static final int SINGLE_CHAR = 0;
    public static final int DECOY = 0;
	public static final int REAL = 1;
	public static final int WRITE_MODEL_PKS = 1;
	public static final int WRITE_PERM_SCORES = 2;
	public static final int WRITE_SCORED_PKS = 3;
	public static final int WRITE_HCD_NONPARAM = 4;
    public static final int WRITE_ALL_MATCHED_PK_SCORES = 5;
    public static final int DEFAULT_RUN_MODE = 0;
    public static final int REPORT_DECOYS = 1;
    public static final int MIN_NUM_NEG_PKS = 50000;
	
	public static final double WATER = 18.01056;
	public static final double PROTON = 1.00728;
	public static final double PPM = 1.0 / 1000000.0;
	public static final double TINY_NUM = 1e-10; // represents a really small number
	public static final double MIN_DELTA_SCORE = 0.1;
    public static final double FUNCTION_TIME_LIMIT = 120; // how long a function call can go for (in seconds) before we kill it

	public static final String MGF_TYPE  = "MGF";
	public static final String MZML_TYPE = "MZML";
	public static final String MZXML_TYPE = "MZXML";

	public static final String PEPTIDE_PROPHET_HEADER  = "peptideprophet";

	public static final Map<String, Double> AA_MASS_MAP = Stream.of(
			new AbstractMap.SimpleEntry<>("A", 71.03711),
			new AbstractMap.SimpleEntry<>("R", 156.10111),
			new AbstractMap.SimpleEntry<>("N", 114.04293),
			new AbstractMap.SimpleEntry<>("D", 115.02694),
			new AbstractMap.SimpleEntry<>("C", 103.00919),
			new AbstractMap.SimpleEntry<>("E", 129.04259),
			new AbstractMap.SimpleEntry<>("Q", 128.05858),
			new AbstractMap.SimpleEntry<>("G", 57.02146),
			new AbstractMap.SimpleEntry<>("H", 137.05891),
			new AbstractMap.SimpleEntry<>("I", 113.08406),
			new AbstractMap.SimpleEntry<>("L", 113.08406),
			new AbstractMap.SimpleEntry<>("K", 128.09496),
			new AbstractMap.SimpleEntry<>("M", 131.04049),
			new AbstractMap.SimpleEntry<>("F", 147.06841),
			new AbstractMap.SimpleEntry<>("P", 97.05276),
			new AbstractMap.SimpleEntry<>("S", 87.03203),
			new AbstractMap.SimpleEntry<>("T", 101.04768),
			new AbstractMap.SimpleEntry<>("W", 186.07931),
			new AbstractMap.SimpleEntry<>("Y", 163.06333),
			new AbstractMap.SimpleEntry<>("V", 99.06841)
	).collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

	public static final Map<String, String> DECOY_AA_MAP = Stream.of(
					new AbstractMap.SimpleEntry<>("2", "A"),
					new AbstractMap.SimpleEntry<>("3", "R"),
					new AbstractMap.SimpleEntry<>("4", "N"),
					new AbstractMap.SimpleEntry<>("5", "D"),
					new AbstractMap.SimpleEntry<>("6", "C"),
					new AbstractMap.SimpleEntry<>("7", "E"),
					new AbstractMap.SimpleEntry<>("8", "Q"),
					new AbstractMap.SimpleEntry<>("9", "G"),
					new AbstractMap.SimpleEntry<>("0", "H"),
					new AbstractMap.SimpleEntry<>("@", "I"),
					new AbstractMap.SimpleEntry<>("#", "L"),
					new AbstractMap.SimpleEntry<>("$", "K"),
					new AbstractMap.SimpleEntry<>("%", "M"),
					new AbstractMap.SimpleEntry<>("&", "F"),
					new AbstractMap.SimpleEntry<>(";", "P"),
					new AbstractMap.SimpleEntry<>("?", "W"),
					new AbstractMap.SimpleEntry<>("~", "V"),
					new AbstractMap.SimpleEntry<>("^", "S"),
					new AbstractMap.SimpleEntry<>("*", "T"),
					new AbstractMap.SimpleEntry<>("=", "Y")
			).collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
}
