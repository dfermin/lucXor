package edu.umich.dmtavt.ptmlocal;

import lucxor.LucxorParams;

/**
 * @author Dmitry Avtonomov
 */
public class PtmLocal {

  public static void main(String[] args) {

    if (args[0].equalsIgnoreCase("-t")) {
      LucxorParams.writeTemplateInputFile();
    }

  }

}
