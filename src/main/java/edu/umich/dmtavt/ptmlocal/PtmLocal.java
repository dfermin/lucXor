package edu.umich.dmtavt.ptmlocal;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;

/**
 * @author Dmitry Avtonomov
 */
public class PtmLocal {

  public static void main(String[] args) throws IOException {

    if (args.length == 0 || Arrays.stream(args).anyMatch("-h"::equals)) {
      System.err.println(usage());
      System.exit(0);
    }

    // write template params file and exit
    if (args[0].equalsIgnoreCase("-t")) {
      final Path path = LucxorParams.writeTemplateInputFile();
      System.err.printf("\nPlease edit the input file: \n\t%s\nwith your favorite text editor\n\n", path);
      System.exit(0);
    }

    final LucxorParams params = LucxorParams.parseParameterFile(args[0]);

  }


  public static String usage() {
    StringBuilder sb = new StringBuilder();
    sb.append("USAGE: java -jar luciphor2.jar <input_file>\n\n");
    sb.append("\tGenerate a luciphor2 input file with: java -jar luciphor2.jar -t\n");
    sb.append("\tModify the input file to suit your needs and submit it to the program.\n");
    sb.append("\tExample: java -jar luciphor2.jar input_file_you_edited\n\n");
    return sb.toString();
  }
}
