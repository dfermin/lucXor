/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor.common;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

import lucxor.utils.Constants;
import lucxor.LucXorConfiguration;
import org.slf4j.LoggerFactory;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.pepxml.PepXmlParser;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.*;

import static lucxor.utils.Constants.AA_MASS_MAP;
import static lucxor.utils.Constants.PEPPROPHET;


/**
 * This class allows to read a {@link PepXmlParser} file and retrieve the corresponding peptides
 * It contains some extra logic to validate if the input pepXML is valid.
 *
 * @author dfermin
 * @author ypriverol
 */
public class PepXML {

	private static final org.slf4j.Logger log = LoggerFactory.getLogger(PepXML.class);

	/**
	 * Read the {@link PepXML} using Stream
	 * @param inputXML input file
	 * @param psmList {@link PSMList} list
	 * @throws IOException IO exception
	 * @throws FileParsingException Parsing exception
	 */
	public static void readStreamFile(File inputXML, PSMList psmList) throws IOException, FileParsingException {
		Path path = Paths.get(inputXML.getAbsolutePath());

		try (final FileInputStream fis = new FileInputStream(path.toString())) {

			Iterator<MsmsRunSummary> it = PepXmlParser.parse(fis);
			boolean readMods = true;
			while (it.hasNext()) {

				MsmsRunSummary runSummary = it.next();

				// PArse PTMs information from Summary
				if(runSummary.getSearchSummary() != null && !runSummary.getSearchSummary().isEmpty()){
					List<SearchSummary> searchSummaryList = runSummary.getSearchSummary();
					PepXML.parseSummaryPTMs(searchSummaryList);
				}

				runSummary.getSpectrumQuery().forEach(query ->{
					long startScan = query.getStartScan();
					String spectrumId = query.getSpectrum();
					int charge = query.getAssumedCharge();
					query.getSearchResult().forEach(results->{

						List<SearchHit> hits = results.getSearchHit()
								.stream()
								.filter(hit -> hit.getHitRank() == 1)
								.collect(Collectors.toList());

						// If the number of rank=1 hits is higher than one, will be supported only for
						// peptidephrophet probability and
						if(hits.size() > 1 && LucXorConfiguration.getScoringMethod() != 0){
							log.info("Multiple peptides reported for the same spectrum, LucXor hasn't been tested for that");
							System.err.println("The current spectrum scan will be skip -- " + startScan);
						}else if(hits.size() > 1){
							// Returns only the hits from peptidephrophet
							hits = hits.stream().filter(hit -> {
								for(AnalysisResult analysisResult: hit.getAnalysisResult()){
									String analysisString = analysisResult.getAnalysis();
									if(analysisString.equalsIgnoreCase(Constants.PEPTIDE_PROPHET_HEADER))
										return true;
								}
								return false;
							}).collect(Collectors.toList());
						}
						if(hits.size() > 1){
							log.info("Multiple peptides reported for the same spectrum, LucXor hasn't been tested" +
									" for that");
							hits.forEach(hit -> System.err
									.println("Warning, more than one peptideprophet candidate -- "
											+ hit.getPeptide()));
						}

						hits.forEach(hit-> {
							PSM curPSM = new PSM();
							curPSM.setSpecId(spectrumId);
curPSM.setScanNum(Long.valueOf(startScan).intValue());
							curPSM.setCharge(charge);
							curPSM.setPeptideSequence(hit.getPeptide());
							List<NameValueType> scores = hit.getSearchScore();
							PepXML.parseScores(curPSM, LucXorConfiguration.getScoringMethod(), scores, hit.getAnalysisResult());
							if (hit.getModificationInfo() != null)
								PepXML.addPTMs(curPSM, hit.getModificationInfo());
							if(PepXML.isValidPSM(curPSM))
								psmList.add(curPSM);
						});

					});
				});

				// The pepXML has the modifications for the search results multiple times
				// (once per spectral file searched). We only need to record these modifications once
				if(readMods) {
					recordModsFromPepXML();
					readMods = false;
				}
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

	/**
	 * This function parse the score of the search_hit in pepXML
	 * @param psm {@link PSM}
	 * @param scores List of {@link NameValueType}
	 */
	public static PSM parseScores(PSM psm, int scoringMethod, List<NameValueType> scores,
								  List<AnalysisResult> results){

		if(scoringMethod != PEPPROPHET){
			scores.forEach( score-> {
				if(LucXorConfiguration.getScoringMethod() ==  Constants.NEGLOGEXPECT && (score.getName()
						.equalsIgnoreCase("expect")))
					psm.setPSMscore(-1.0 * Math.log(Double.parseDouble(score.getValueStr())));
				if(LucXorConfiguration.getScoringMethod() ==  Constants.MASCOTIONSCORE && (score.getName()
						.equalsIgnoreCase("ionscore")))
					psm.setPSMscore(Double.parseDouble(score.getValueStr()));
				if(LucXorConfiguration.getScoringMethod() ==  Constants.XTDHYPERSCORE && (score.getName()
						.equalsIgnoreCase("hyperscore")))
					psm.setPSMscore(Double.parseDouble(score.getValueStr()));
				if(LucXorConfiguration.getScoringMethod() ==  Constants.XCORR && (score.getName()
						.equalsIgnoreCase("xcorr")))
					psm.setPSMscore(Double.parseDouble(score.getValueStr()));
			});
		}else{
			for(AnalysisResult analysisResult: results){
				String analysisString = analysisResult.getAnalysis();
				if(analysisString.equalsIgnoreCase(Constants.PEPTIDE_PROPHET_HEADER)){
					List<Object> probabilities = analysisResult.getAny();
					probabilities.forEach( x -> {
						if (x instanceof PeptideprophetResult)
							psm.setPSMscore(((PeptideprophetResult)x).getProbability());
					});
				}
			}
		}
		return psm;
	}

	/**
	 * Check if one PSM is valid (do not contains aminoacid variants)
	 * @param psm {@link PSM}
	 * @return True if not variants are present
	 */
	public static boolean isValidPSM(PSM psm) {

		// skip PSMs with non-standard amino acid characters
		String x = psm.getPeptideSequence();
		int numBadChars = 0;
		for(int i = 0; i < x.length(); i++) {
			String c = Character.toString(x.charAt(i));
			if( !"ACDEFGHIKLMNPQRSTVWY".contains(c) ) numBadChars++;
		}

		// Skip PSMs that exceed the number of candidate permutations
		if(psm.getOrigPep().getNumPerm() > LucXorConfiguration.getMaxNumPermutations())
			numBadChars = 100;

		if(numBadChars == 0)
			psm.process();

		return psm.isKeeper();
	}

	/**
	 * Parse general summary modifications
	 * @param searchSumaryList
	 */
	public static void parseSummaryPTMs(List<SearchSummary> searchSumaryList){

		// for handling terminal modifications (this is found in the top portion of the pepXML file)
		searchSumaryList.forEach( x-> {
			List<TerminalModification> terminalMods = x.getTerminalModification();
			if(terminalMods != null && !terminalMods.isEmpty()){
				terminalMods.forEach( termMod -> {
					String terminus = termMod.getTerminus();
					double modMass = termMod.getMassdiff();
					if(terminus.equalsIgnoreCase("n"))
						LucXorConfiguration.setNTermMass(modMass);
					if(terminus.equalsIgnoreCase("c"))
						LucXorConfiguration.setCTermMass(modMass);
				});
			}
			List<AminoacidModification> aminoMods = x.getAminoacidModification();
			if(aminoMods != null && !aminoMods.isEmpty()){
				aminoMods.forEach( aminoMod -> {
					String aa = aminoMod.getAminoacid();
					double modMass = aminoMod.getMassdiff();
					String varStatus = aminoMod.getVariable();

					String AA_alphabet = "ACDEFGHIKLMNPQRSTVWY";
					if( AA_alphabet.contains(aa) ) {

						// if this is a valid AA character that is a variable modification, record it
						// as a lower case character in the varModMap
						if(varStatus.equalsIgnoreCase("y")) {
							if(!LucXorConfiguration.getVarModMap().containsKey(aa.toLowerCase())) {
								LucXorConfiguration.getVarModMap().put(aa.toLowerCase(), modMass);
							}
						}else {
							if(!LucXorConfiguration.getFixedModMap().containsKey(aa.toUpperCase())) {
								LucXorConfiguration.getFixedModMap().put(aa.toUpperCase(), modMass);
							}
						}
					}
				});
			}
		});


	}

	public static PSM addPTMs(PSM psm, ModificationInfo mods){

		// for handling terminal modifications (this is found in the top portion of the pepXML file)
		if(mods.getModCtermMass() != null){
			double cterm = mods.getModCtermMass();
            psm.getModCoordMap().put(Constants.CTERM_MOD, cterm);
		}
		if(mods.getModNtermMass() != null){
			double nterm = mods.getModNtermMass();
			psm.getModCoordMap().put(Constants.NTERM_MOD, nterm);
		}

		mods.getModAminoacidMass().forEach( massMod -> {
			int pos = massMod.getPosition() - 1; // ensures zero-based coordiantes
			double mass = massMod.getMass();
			psm.getModCoordMap().put(pos, mass);
		});

		return psm;
	}
}
