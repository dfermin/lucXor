/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.AnalysisResult;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.EngineType;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.ModAminoacidMass;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.ModificationInfo;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.MsmsPipelineAnalysis;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.MsmsRunSummary;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.NameValueType;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.PeptideprophetResult;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.SearchHit;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.SearchResult;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.SearchSummary;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.TerminalModification;


/**
 * @author dfermin
 */
class PepXML extends DefaultHandler {


  public static final String AA_CHARS_UPPER = "ACDEFGHIKLMNPQRSTVWY".toUpperCase();
  public static final RangeSet<Integer> AA_CHARS_UPPER_RS;
  public static final String AA_CHARS_LOWER = AA_CHARS_UPPER.toLowerCase();
  public static final RangeSet<Integer> AA_CHARS_LOWER_RS;
  public static final String PEPXML_ANALYSIS_PEPTIDEPROPHET = "peptideprophet";
  private HashMap<String, Double> variableMods;
  private HashMap<String, Double> fixedMods;
  private String searchEngine;
  private PSM curPSM = null;
  private boolean recordMods;
  private File inputPepXml = null;

  static {
    {
      RangeSet<Integer> rs = TreeRangeSet.create();
      for (char c : AA_CHARS_UPPER.toCharArray()) {
        int code = (int) c;
        rs.add(Range.closed(code, code));
      }
      AA_CHARS_UPPER_RS = rs;
    }
    {
      RangeSet<Integer> rs = TreeRangeSet.create();
      for (char c : AA_CHARS_LOWER.toCharArray()) {
        int code = (int) c;
        rs.add(Range.closed(code, code));
      }
      AA_CHARS_LOWER_RS = rs;
    }
  }

  // Default constructor for this class
  PepXML() {
    variableMods = new HashMap<>();
    fixedMods = new HashMap<>();
    recordMods = true; // changes to false after the end of first search_summary section
  }

  public void load(File inputPepXml)
      throws ParserConfigurationException, SAXException, IOException {
    if (this.inputPepXml != null) {
      throw new IllegalStateException("Calling load() twice on PepXML instance is not allowed.");
    }
    this.inputPepXml = inputPepXml;
    SAXParserFactory factory = SAXParserFactory.newInstance();
    factory.setValidating(false);
    SAXParser parser = factory.newSAXParser();
    parser.parse(inputPepXml,
        this); // this class itself is the default handler, hence the use of 'this'
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

      if (AA_CHARS_UPPER.contains(aa)) {
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

      try {
        // skip PSMs with non-standard amino acid characters
        String x = curPSM.origPep.peptide;
        for (int i = 0; i < x.length(); i++) {
          String c = Character.toString(x.charAt(i));
          if (!AA_CHARS_UPPER.contains(c)) {
            return;
          }
        }

        // Skip PSMs that exceed the number of candidate permutations
        if (curPSM.origPep.getNumPerm() > Globals.max_num_permutations) {
          return;
        }

        curPSM.init();
        if (curPSM.isKeeper) {
          Globals.PSM_list.add(curPSM);
        }

      } finally {
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

  public void load(MsmsPipelineAnalysis pepxml) {

    // gather info about search engines
    final List<EngineType> engineTypes = pepxml.getMsmsRunSummary().stream()
        .flatMap(msmsRunSummary -> msmsRunSummary.getSearchSummary().stream())
        .map(SearchSummary::getSearchEngine).distinct().collect(Collectors.toList());
    if (engineTypes.size() > 1)
      throw new IllegalStateException("More than 1 search engine type detected in pepxml, not allowed.");
    if (engineTypes.size() == 1)
      searchEngine = engineTypes.get(0).value().toLowerCase();

    // gather info about modifications
    // This is for handling the modification description lines at the beginning of a pepXML file
    pepxml.getMsmsRunSummary().stream()
        .map(MsmsRunSummary::getSearchSummary).flatMap(Collection::stream)
        .map(SearchSummary::getAminoacidModification).flatMap(Collection::stream)
        .forEach(aaMod -> {
          final String aa = aaMod.getAminoacid();
          final double massdiff = aaMod.getMassdiff();
          final String variable = aaMod.getVariable();
          if (AA_CHARS_UPPER.contains(aa)) {
            // if this is a valid AA character that is a variable modification, record it
            // as a lower case character in the varModMap
            if ("Y".equals(variable) || "y".equals(variable)) {
              Globals.varModMap.put(aa.toLowerCase(), massdiff);
            } else {
              Globals.fixedModMap.put(aa.toUpperCase(), massdiff);
            }
          }
        });


    // terminal modifications

    // N-term
    final List<Mod> nTermMods = pepxml.getMsmsRunSummary().stream()
        .map(MsmsRunSummary::getSearchSummary).flatMap(Collection::stream)
        .map(SearchSummary::getTerminalModification).flatMap(Collection::stream)
        .map(tm -> new Mod(tm.getTerminus(), equalsYy(tm.getProteinTerminus()),
            equalsYy(tm.getVariable()), tm.getMassdiff()))
        .distinct().filter(mod -> equalsXx("n", "N", mod.terminus))
        .collect(Collectors.toList());

    if (nTermMods.size() == 1) {
      Globals.ntermMass = nTermMods.get(0).massDiff;

    } else if (nTermMods.size() == 2) {
      final Mod tm1 = nTermMods.get(0);
      final Mod tm2 = nTermMods.get(1);
      if (tm1.isProtein ^ tm2.isProtein) {
        System.err.println("Two N-term mods found, but one as 'protein terminus' only, so it's ok");
        // If one is 'protein terminus == Y' and the other is 'protein terminus == N' - it's ok.
        // Otherwise it's not ok.
        Globals.ntermMassProt = tm1.isProtein ? tm1.massDiff : tm2.massDiff;
        Globals.ntermMass = tm1.isProtein ? tm2.massDiff : tm1.massDiff;

      } else {
        throw new IllegalStateException("More than one distinct N-term mod mass found with "
            + "the same site-specificity. Not allowed.");
      }

    } else if (!nTermMods.isEmpty()) {
      throw new IllegalStateException("More than two distinct N-term mod mass found, not allowed.");
    }


    // C-term
    final List<Mod> cTermMods = pepxml.getMsmsRunSummary().stream()
        .map(MsmsRunSummary::getSearchSummary).flatMap(Collection::stream)
        .map(SearchSummary::getTerminalModification).flatMap(Collection::stream)
        .map(tm -> new Mod(tm.getTerminus(), equalsYy(tm.getProteinTerminus()),
            equalsYy(tm.getVariable()), tm.getMassdiff()))
        .distinct().filter(mod -> equalsXx("c", "C", mod.terminus))
        .collect(Collectors.toList());

    if (cTermMods.size() == 1) {
      Globals.ctermMass = nTermMods.get(0).massDiff;

    } else if (cTermMods.size() == 2) {
      final Mod tm1 = cTermMods.get(0);
      final Mod tm2 = cTermMods.get(1);
      if (tm1.isProtein ^ tm2.isProtein) {
        System.err.println("Two C-term mods found, but one as 'protein terminus' only, so it's ok");
        // If one is 'protein terminus == Y' and the other is 'protein terminus == N' - it's ok.
        // Otherwise it's not ok.
        Globals.ctermMassProt = tm1.isProtein ? tm1.massDiff : tm2.massDiff;
        Globals.ctermMass = tm1.isProtein ? tm2.massDiff : tm1.massDiff;

      } else {
        throw new IllegalStateException("More than one distinct C-term mod mass found with "
            + "the same site-specificity. Not allowed.");
      }

    } else if (!cTermMods.isEmpty()) {
      throw new IllegalStateException("More than two distinct C-term mod mass found, not allowed.");
    }


    if (recordMods) {
      Globals.recordModsFromPepXML();
      recordMods = false;
    } else {
      throw new IllegalStateException("Called Globals.recordModsFromPepXML() twice, this should not happen");
    }

    // PSMs - going through "spectrum_query" entries
    pepxml.getMsmsRunSummary().stream()
        .map(MsmsRunSummary::getSpectrumQuery).flatMap(Collection::stream)
        .forEach(sq -> { // sq == spectrum_query
          if (sq.getSearchResult().isEmpty())
            return;
          if (sq.getSearchResult().size() > 1)
            throw new IllegalStateException("More than one 'search_result' found in 'spectrum_query'. Not supported.");

          final SearchResult sr = sq.getSearchResult().get(0);

          if (sr.getSearchHit().isEmpty())
            return;
          if (sr.getSearchHit().size() > 1)
            throw new IllegalStateException("More than one 'search_hit' found in 'search_result'. Not supported.");
          SearchHit sh = sr.getSearchHit().get(0);

          PSM psm = new PSM();
          psm.specId = sq.getSpectrum();
          if (sq.getStartScan() >= Integer.MAX_VALUE)
            throw new IllegalStateException("Found spectrum_query with scan number larger than Integer.MAX_VALUE. Not supported.");
          psm.scanNum = (int)sq.getStartScan();
          if (sq.getAssumedCharge() == null)
            throw new IllegalStateException("Found spectrum_query with assumed_charge == null. Not supported.");
          psm.charge = sq.getAssumedCharge();

          // skip peptides with non-standard AAs
          for (char c : sh.getPeptide().toCharArray()) {
            if (!AA_CHARS_UPPER_RS.contains((int)c))
              return;
          }

          psm.origPep.peptide = sh.getPeptide();

          final ModificationInfo modInfo = sh.getModificationInfo();
          if (modInfo != null) {
            if (modInfo.getModNtermMass() != null) {
              psm.modCoordMap.put(Constants.NTERM_MOD, modInfo.getModNtermMass());
            }
            if (modInfo.getModCtermMass() != null) {
              psm.modCoordMap.put(Constants.CTERM_MOD, modInfo.getModCtermMass());
            }

            // record masses of modified amino acid residues
            for (ModAminoacidMass modAa : modInfo.getModAminoacidMass()) {
              if (modAa.getPosition() == null)
                throw new IllegalStateException("Position of a 'mod_aminoacid_mass' was null. Not supported.");
              int pos = modAa.getPosition() - 1;
              double mass = modAa.getMass();
              psm.modCoordMap.put(pos, mass);
            }
          }


          // scores
          switch (Globals.scoringMethod) {

            case Constants.PEPPROPHET:
              boolean isFound = false;
              for (AnalysisResult ar : sh.getAnalysisResult()) {
                if (!PEPXML_ANALYSIS_PEPTIDEPROPHET.equals(ar.getAnalysis()))
                  continue;
                isFound = true;
                if (ar.getAny().isEmpty())
                  throw new IllegalStateException("Empty peptideprophet analysis not supported.");
                PeptideprophetResult ppr = (PeptideprophetResult)ar.getAny().get(0);
                psm.PSMscore = ppr.getProbability();
              }
              if (!isFound)
                throw new IllegalStateException("Scoring method was peptide prophet, but "
                    + "<analysis_result analysis=\"peptideprophet\"> not found for at least "
                    + "one 'search_hit' entry.");
              break;

            case Constants.NEGLOGEXPECT:
              psm.PSMscore = findSearchScore(sh, "expect");
              break;

            case Constants.MASCOTIONSCORE:
              psm.PSMscore = findSearchScore(sh, "ionscore");
              break;

            case Constants.XTDHYPERSCORE:
              psm.PSMscore = findSearchScore(sh, "hyperscore");
              break;

            case Constants.XCORR:
              psm.PSMscore = findSearchScore(sh, "xcorr");
              break;
          }
        });
  }

  private static boolean equalsYy(String s) {
    return "Y".equals(s) || "y".equals(s);
  }

  private static boolean equalsNn(String s) {
    return "N".equals(s) || "n".equals(s);
  }

  private static boolean equalsXx(String c1, String c2, String s) {
    return c1.equals(s) || c2.equals(s);
  }

  private static double findSearchScore(SearchHit sh, String scoreName) {
    for (NameValueType nv : sh.getSearchScore()) {
      if (nv.getName().equals(scoreName))
        return Double.parseDouble(nv.getValueStr());
    }
    throw new IllegalStateException(String.format("Requested search score '%s' not found for at least one search_hit", scoreName));
  }

  private class Mod {
    final String terminus;
    final boolean isProtein;
    final boolean isVariable;
    final double massDiff;

    public Mod(String terminus, boolean isProtein, boolean isVariable, double massDiff) {
      this.terminus = terminus;
      this.isProtein = isProtein;
      this.isVariable = isVariable;
      this.massDiff = massDiff;
    }

    @Override
    public boolean equals(Object o) {
      if (this == o) {
        return true;
      }
      if (o == null || getClass() != o.getClass()) {
        return false;
      }

      Mod mod = (Mod) o;

      if (isProtein != mod.isProtein) {
        return false;
      }
      if (isVariable != mod.isVariable) {
        return false;
      }
      if (Double.compare(mod.massDiff, massDiff) != 0) {
        return false;
      }
      return terminus != null ? terminus.equals(mod.terminus) : mod.terminus == null;
    }

    @Override
    public int hashCode() {
      int result;
      long temp;
      result = terminus != null ? terminus.hashCode() : 0;
      result = 31 * result + (isProtein ? 1 : 0);
      result = 31 * result + (isVariable ? 1 : 0);
      temp = Double.doubleToLongBits(massDiff);
      result = 31 * result + (int) (temp ^ (temp >>> 32));
      return result;
    }
  }
}
