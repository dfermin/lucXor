/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucXor;

/**
 *
 * @author dfermin
 */
public class constants {
	public static final int PEPPROPHET = 0;
	public static final int MASCOTIONSCORE = 1;
	public static final int NEGLOGEXPECT = 2;
    public static final int XTDHYPERSCORE = 3;
    public static final int XCORR = 4;
	public static final int DALTONS = 0;
	public static final int PPM_UNITS = 1;
	public static final int PEPXML = 0;
	public static final int TSV = 1;
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


}
