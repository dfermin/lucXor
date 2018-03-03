/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;


/**
 * @author dfermin
 */
class PepXML extends DefaultHandler {


  private HashMap<String, Double> variableMods;
  private HashMap<String, Double> fixedMods;
  private String searchEngine;
  private final String AA_alphabet = "ACDEFGHIKLMNPQRSTVWY";
  private PSM curPSM = null;
  private boolean recordMods;
  private File inputPepXml = null;

  // Default constructor for this class
  public PepXML() {
    variableMods = new HashMap<>();
    fixedMods = new HashMap<>();
    recordMods = true; // changes to false after the end of first search_summary section
  }

  public void load(File inputPepXml) throws ParserConfigurationException, SAXException, IOException {
    if (this.inputPepXml != null)
      throw new IllegalStateException("Calling load() twice on PepXML instance is not allowed.");
    this.inputPepXml = inputPepXml;
    SAXParserFactory factory = SAXParserFactory.newInstance();
    factory.setValidating(false);
    SAXParser parser = factory.newSAXParser();
    parser.parse(inputPepXml, this); // this class itself is the default handler, hence the use of 'this'
  }

  public void startElement(String uri, String localName, String qName, Attributes attr) {

    if (qName.equalsIgnoreCase("search_summary")) {
      String se = attr.getValue("search_engine");
      if (se.startsWith("Hydra")) {
        searchEngine = "hydra";
      } else if (se.equalsIgnoreCase("MASCOT")) {
        searchEngine = "mascot";
      } else if (se.contains("Tandem")) {
        searchEngine = "tandem";
      } else if (se.equalsIgnoreCase("Comet")) {
        searchEngine = "comet";
      } else if (se.equalsIgnoreCase("Sequest")) {
        searchEngine = "sequest";
      }
    }


    if (qName.equalsIgnoreCase("aminoacid_modification")) {
      // This is for handling the modification description lines at the beginning of a pepXML file
      String aa = attr.getValue("aminoacid");
      double modMass = Double.valueOf(attr.getValue("massdiff"));
      String varStatus = attr.getValue("variable");

      if (AA_alphabet.contains(aa)) {
        // if this is a valid AA character that is a variable modification, record it
        // as a lower case character in the varModMap
        if (varStatus.equalsIgnoreCase("y")) {
          if (!Globals.varModMap.containsKey(aa.toLowerCase())) {
            Globals.varModMap.put(aa.toLowerCase(), modMass);
          }
        } else {
          if (!Globals.fixedModMap.containsKey(aa.toUpperCase())) {
            Globals.fixedModMap.put(aa.toUpperCase(), modMass);
          }
        }
      }


    } else if (qName.equalsIgnoreCase("terminal_modification")) {
      // for handling terminal modifications (this is found in the top portion of the pepXML file)
      String terminus = attr.getValue("terminus");
      double modMass = Double.valueOf(attr.getValue("massdiff"));

      String aa = "x";

      if (terminus.equalsIgnoreCase("n")) {
        Globals.ntermMass = modMass;
      }
      if (terminus.equalsIgnoreCase("c")) {
        Globals.ctermMass = modMass;
      }


    } else if (qName.equalsIgnoreCase("spectrum_query")) {
      curPSM = new PSM();

      curPSM.specId = attr.getValue("spectrum");
      curPSM.scanNum = Integer.valueOf(attr.getValue("start_scan"));
      curPSM.charge = Integer.valueOf(attr.getValue("assumed_charge"));


    } else if (qName.equalsIgnoreCase("search_hit")) {
      curPSM.origPep.peptide = attr.getValue("peptide");


    } else if (qName.equalsIgnoreCase("modification_info")) {
      if (attr.getLocalName(0).equalsIgnoreCase("mod_nterm_mass")) {
        double nterm = Double.valueOf(attr.getValue("mod_nterm_mass"));
        curPSM.modCoordMap.put(Constants.NTERM_MOD, nterm);
      } else if (attr.getLocalName(0).equalsIgnoreCase("mod_cterm_mass")) {
        double cterm = Double.valueOf(attr.getValue("mod_cterm_mass"));
        curPSM.modCoordMap.put(Constants.CTERM_MOD, cterm);
      }
    } else if (qName.equalsIgnoreCase("mod_aminoacid_mass")) {
      // record the mass of a modified amino acid residue
      int pos = Integer.valueOf(attr.getValue("position")) - 1; // ensures zero-based coordiantes
      double mass = Double.valueOf(attr.getValue("mass"));
      curPSM.modCoordMap.put(pos, mass);


    } else if (qName.equalsIgnoreCase("search_score")) {

      if (Globals.scoringMethod != Constants.PEPPROPHET) {
        for (int i = 0; i < attr.getLength() - 1; i++) {
          int j = i + 1;
          String k = attr.getValue(i);
          double score = Double.valueOf(attr.getValue(j));

          if (Globals.scoringMethod == Constants.NEGLOGEXPECT) {
            if (k.equalsIgnoreCase("expect")) {
              curPSM.PSMscore = -1.0 * Math.log(score);
            }
          }

          if (Globals.scoringMethod == Constants.MASCOTIONSCORE) {
            if (k.equalsIgnoreCase("ionscore")) {
              curPSM.PSMscore = score;
            }
          }

          if (Globals.scoringMethod == Constants.XTDHYPERSCORE) {
            if (k.equalsIgnoreCase("hyperscore")) {
              curPSM.PSMscore = score;
            }
          }

          if (Globals.scoringMethod == Constants.XCORR) {
            if (k.equalsIgnoreCase("xcorr")) {
              curPSM.PSMscore = score;
            }
          }
        }

      }


    } else if (qName.equalsIgnoreCase("peptideprophet_result")) {
      if (Globals.scoringMethod == Constants.PEPPROPHET) {
        curPSM.PSMscore = Double.valueOf(attr.getValue("probability"));
      }
    }
  }


  public void endElement(String uri, String localName, String qName) {

    if (qName.equalsIgnoreCase("search_summary")) { // record the AA modifications

      // The pepXML has the modifications for the search results multiple times
      // (once per spectral file searched). We only need to record these modifications once
      if (recordMods) {
        Globals.recordModsFromPepXML();
        recordMods = false;
      }
    } else if (qName.equalsIgnoreCase("search_hit")) { // end of record

      // skip PSMs with non-standard amino acid characters
      String x = curPSM.origPep.peptide;
      int numBadChars = 0;
      for (int i = 0; i < x.length(); i++) {
        String c = Character.toString(x.charAt(i));
        if (!"ACDEFGHIKLMNPQRSTVWY".contains(c)) {
          numBadChars++;
        }
      }

      // Skip PSMs that exceed the number of candidate permutations
      if (curPSM.origPep.getNumPerm() > Globals.max_num_permutations) {
        numBadChars = 100;
      }

      if (numBadChars == 0) {
        curPSM.process();
        if (curPSM.isKeeper) {
          Globals.PSM_list.add(curPSM);
        }
        curPSM = null;
      }
    }

  }


  /**************************************************************************
   // This is a critical function to the SAX parser.
   // It handles non-XML tag text when it's encountered
   */
  public void characters(char[] buffer, int start, int length) {
  }

}
