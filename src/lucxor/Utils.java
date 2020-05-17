package lucxor;

public class Utils {

  /**
   * This function parse an string line and results
   * @param line line to be parse
   * @return Result string
   */
  static String parseInputLine(String line) {
    String ret = "";
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
  public static String generateIndex(String srcFile, Integer scanNum){
    return srcFile + ":" + scanNum;
  }
}
