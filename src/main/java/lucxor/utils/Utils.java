package lucxor.utils;

import static lucxor.utils.Constants.AA_MASS_MAP;
import static lucxor.utils.Constants.DECOY_AA_MAP;

public class Utils {

  /**
   * This function parse an string line and results
   * @param line line to be parse
   * @return Result string
   */
  public static String parseInputLine(String line) {
    String ret;
    StringBuilder sb = new StringBuilder();
    int N = line.length();

    int b = line.indexOf('=') + 1;

    for(int i = b; i < N; i++) {
      char c = line.charAt(i);
      if(c == '#') break;
      if(c == ' ') continue;

      sb.append(c);
    }
    ret = sb.toString();
    return ret;
  }

  /**
   * Generate an String Filename:Scan
   * @param srcFile Spectrum File
   * @param scanNum Scan number
   * @return  Filename:Scan
   */
  public static String generateIndex(String srcFile, long scanNum){
    return srcFile + ":" + scanNum;
  }

  /**
   * Get FragmentIon Mass
   * @param x String
   * @param z fragment
   * @param addl_mass
   * @return value
   */
  public static double getFragmentIonMass(String x, double z, double addl_mass) {
    double ret = 1.00728 * z;

    int start = x.indexOf(':') + 1;
    int stop  = x.length();

    if(x.contains("-")) stop = x.indexOf('-');

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

  /**
   * Function returns the decoy 'key' from the decoy map for the given residue
   * @param c
   * @return
   */
  public static String getDecoySymbol(char c) {
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

  /**
   *  Function returns the TPP-formatted representation of the given single
   *  character modification
   * @param c
   * @param ntermMass
   * @param ctermMass
   * @return
   */
  public static String getTPPresidue(String c, double ntermMass, double ctermMass) {
    String ret = "";
    String orig;

    if(c.equalsIgnoreCase("[")) {
      int d = (int) Math.round(ntermMass) + 1; // adds a proton
      ret = "n[" + d + "]";
    }
    else if(c.equals("]")) {
      int d = (int) Math.round(ctermMass);
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


  /**
   * Function returns true if the given sequence contains decoy characters
   * @param seq Sequence
   * @return if decoy true
   */
  public static boolean isDecoySeq(String seq) {
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


  /**
   * Function returns true if the given residue is from the decoy list
   * @param AA check is aminoa cid is decoy list
   * @return return true if decoy
   */
  public static boolean isDecoyResidue(String AA) {
    boolean ret = false;
    if(DECOY_AA_MAP.containsKey(AA))
      ret = true;
    return ret;
  }
}
