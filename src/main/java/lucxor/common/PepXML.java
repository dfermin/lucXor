/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor.common;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import lucxor.utils.Constants;
import lucxor.LucXorConfiguration;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import static lucxor.utils.Constants.AA_MASS_MAP;


/**
 *
 * @author dfermin
 */
public class PepXML extends DefaultHandler {

	private final PSMList psmList;
	private String temp;
	private PSM curPSM = null;
	private boolean recordMods = true; // changes to false after the end of first search_summary section

	public static void readPepXMLFile(File inputXML, PSMList psmList) throws ParserConfigurationException, SAXException, IOException {
		new PepXML(inputXML, psmList);
	}
	
	private PepXML(File inputXML, PSMList psmList) throws ParserConfigurationException, SAXException, IOException {

		SAXParserFactory factory = SAXParserFactory.newInstance();
		factory.setValidating(false);
		SAXParser parser = factory.newSAXParser();
		this.psmList = psmList;
		parser.parse(inputXML, this); // this class itself is the default handler, hence the use of 'this'
	}
	
	public void startElement(String uri, String localName, String qName, Attributes attr) {
		temp = "";
		
		
		if(qName.equalsIgnoreCase("search_summary")) {
			String se = attr.getValue("search_engine");
			String searchEngine;
			if(se.startsWith("Hydra")) searchEngine = "hydra";
			else if(se.equalsIgnoreCase("MASCOT")) searchEngine = "mascot";
			else if(se.contains("Tandem")) searchEngine = "tandem";
			else if(se.equalsIgnoreCase("Comet")) searchEngine = "comet";
			else if(se.equalsIgnoreCase("Sequest")) searchEngine = "sequest";
		}

        // This is for handling the modification description lines at the beginning of a pepXML file
		if(qName.equalsIgnoreCase("aminoacid_modification")) {
			String aa = attr.getValue("aminoacid");
			double modMass = Double.valueOf( attr.getValue("massdiff") );
			String varStatus = attr.getValue("variable");

			String AA_alphabet = "ACDEFGHIKLMNPQRSTVWY";
			if( AA_alphabet.contains(aa) ) {
				// if this is a valid AA character that is a variable modification, record it
				// as a lower case character in the varModMap
				if(varStatus.equalsIgnoreCase("y")) {
					if(!LucXorConfiguration.getVarModMap().containsKey(aa.toLowerCase())) {
						LucXorConfiguration.getVarModMap().put(aa.toLowerCase(), modMass);
					}
				}
				else {
					if(!LucXorConfiguration.getFixedModMap().containsKey(aa.toUpperCase())) {
						LucXorConfiguration.getFixedModMap().put(aa.toUpperCase(), modMass);
					}
				}
			}
		}
		
		// for handling terminal modifications (this is found in the top portion of the pepXML file)
		if(qName.equalsIgnoreCase("terminal_modification")) {
			String terminus = attr.getValue("terminus");
			double modMass = Double.valueOf( attr.getValue("massdiff") );
			
			if(terminus.equalsIgnoreCase("n"))
				LucXorConfiguration.setNTermMass(modMass);
			if(terminus.equalsIgnoreCase("c"))
				LucXorConfiguration.setCTermMass(modMass);
		}
		
		if(qName.equalsIgnoreCase("spectrum_query")) {
			curPSM = new PSM();
			
			curPSM.setSpecId(attr.getValue("spectrum"));
			curPSM.setScanNum(Integer.valueOf( attr.getValue("start_scan") ));
			curPSM.setCharge(Integer.valueOf( attr.getValue("assumed_charge") ));
		}
		
		if(qName.equalsIgnoreCase("search_hit")) {
			curPSM.setPeptideSequence(attr.getValue("peptide"));
		}
		
		if(qName.equalsIgnoreCase("modification_info")) {
			if(attr.getLocalName(0).equalsIgnoreCase("mod_nterm_mass")) {
				double nterm = Double.valueOf(attr.getValue("mod_nterm_mass"));
				curPSM.getModCoordMap().put(Constants.NTERM_MOD, nterm);
			}
			else if(attr.getLocalName(0).equalsIgnoreCase("mod_cterm_mass")) {
				double cterm = Double.valueOf(attr.getValue("mod_cterm_mass"));
				curPSM.getModCoordMap().put(Constants.CTERM_MOD, cterm);
			}
		}
		
		// record the mass of a modified amino acid residue
		if(qName.equalsIgnoreCase("mod_aminoacid_mass")) {
			int pos = Integer.valueOf( attr.getValue("position") ) - 1; // ensures zero-based coordiantes
			double mass = Double.valueOf( attr.getValue("mass") );
			curPSM.getModCoordMap().put(pos, mass);
		}
		
		if(qName.equalsIgnoreCase("search_score")) {
			
			if(LucXorConfiguration.getScoringMethod() != Constants.PEPPROPHET) {
				for(int i = 0; i < attr.getLength() - 1; i++) {
					int j = i + 1;
					String k = attr.getValue(i);
					double score = Double.valueOf( attr.getValue(j) );
					
					if(LucXorConfiguration.getScoringMethod() == Constants.NEGLOGEXPECT) {
						if(k.equalsIgnoreCase("expect")) curPSM.setPSMscore(-1.0 * Math.log(score));
					}
					
					if(LucXorConfiguration.getScoringMethod() == Constants.MASCOTIONSCORE) {
						if(k.equalsIgnoreCase("ionscore")) curPSM.setPSMscore(score);
					}

                    if(LucXorConfiguration.getScoringMethod() == Constants.XTDHYPERSCORE) {
                        if(k.equalsIgnoreCase("hyperscore")) curPSM.setPSMscore(score);
                    }

                    if(LucXorConfiguration.getScoringMethod() == Constants.XCORR) {
                        if (k.equalsIgnoreCase("xcorr")) curPSM.setPSMscore(score);
                    }
				}
				
			}
		}
		
		if(qName.equalsIgnoreCase("peptideprophet_result")) {
			if(LucXorConfiguration.getScoringMethod() == Constants.PEPPROPHET)
				curPSM.setPSMscore(Double.valueOf(attr.getValue("probability")));
		}
	}


	public void endElement(String uri, String localName, String qName) {
		
		temp = null; // shouldn't need this anymore
		
		if(qName.equalsIgnoreCase("search_summary")) { // record the AA modifications
			
				// The pepXML has the modifications for the search results multiple times
				// (once per spectral file searched). We only need to record these modifications once
				if(recordMods) {
					recordModsFromPepXML();
					recordMods = false;
				}
		}
		
		
		if(qName.equalsIgnoreCase("search_hit")) { // end of record
			
			// skip PSMs with non-standard amino acid characters
			String x = curPSM.getPeptideSequence();
			int numBadChars = 0;
			for(int i = 0; i < x.length(); i++) {
				String c = Character.toString(x.charAt(i));
				if( !"ACDEFGHIKLMNPQRSTVWY".contains(c) ) numBadChars++;
			}
			
			// Skip PSMs that exceed the number of candidate permutations
			if(curPSM.getOrigPep().getNumPerm() > LucXorConfiguration.getMaxNumPermutations())
				numBadChars = 100;

			if(numBadChars == 0) {
                curPSM.process();
				if(curPSM.isKeeper()) psmList.add(curPSM);
                curPSM = null;
			}
		}

	}


	/**
	 * This function is only called if we are getting our modifications from
	 * a pepXML file.
	 */
	private static void recordModsFromPepXML() {

		String alphabet = "ACDEFGHIKLMNPQRSTVWY";

		for(Map.Entry<String, Double> stringDoubleEntry : LucXorConfiguration.getFixedModMap().entrySet()) {

			if(!alphabet.contains(stringDoubleEntry.getKey())) continue; // skip non-amino acid characters

			double mass = AA_MASS_MAP.get(stringDoubleEntry.getKey()) + stringDoubleEntry.getValue();
			String symbol = stringDoubleEntry.getKey().toUpperCase();
			AA_MASS_MAP.put(symbol, mass);
		}

		// First filter the data in VAR_MOD_MAP.
		// We want to remove non-standard amino acid characters and remove
		// the amino acids that are in our 'TARGET_MOD_MAP' variable.
		HashMap<String, Double> tmp = new HashMap<>(LucXorConfiguration.getVarModMap());
		LucXorConfiguration.clearVarMods();

		for(Map.Entry<String, Double> stringDoubleEntry : tmp.entrySet()) {
			String C = stringDoubleEntry.getKey().toUpperCase();
			double mass = stringDoubleEntry.getValue();
			if(!alphabet.contains(C)) continue; // skip non-amino acid characters
			if(LucXorConfiguration.getTargetModMap().containsKey(C)) continue; // skip target residues
			LucXorConfiguration.addVarMod(stringDoubleEntry.getKey(), mass);
		}

		for(Map.Entry<String, Double> stringDoubleEntry : LucXorConfiguration.getVarModMap().entrySet()) {
			String C = stringDoubleEntry.getKey().toUpperCase();
			double mass = AA_MASS_MAP.get(C) + stringDoubleEntry.getValue();
			AA_MASS_MAP.put(stringDoubleEntry.getKey(), mass);
		}
	}
	
	
	/**************************************************************************
	// This is a critical function to the SAX parser.
	// It handles non-XML tag text when it's encountered
	*/
	public void characters(char[] buffer, int start, int length) {
		temp = new String(buffer, start, length);
	}
	
}
