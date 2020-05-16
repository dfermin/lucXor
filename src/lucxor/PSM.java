/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import lombok.extern.slf4j.Slf4j;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author dfermin
 */
@Slf4j
class PSM {

	private String specId; // TPP-based name
	private String srcFile; // name of file from which spectrum is derived
	private int scanNum;
	private int charge;
	private double PSMscore; // score assigned to this PSM by original search engine
	private double deltaScore;
	private double localFDR;  // local false localization rate
	private double globalFDR; // global false localization rate
	private boolean isKeeper; // true means the PSM can be scored
	private boolean useForModel; // true means the PSM can be used for modeling
	private boolean isDecoy; // true means the top match for this PSM is a decoy
	private boolean isUnambiguous; // true means the number potential PTM sites is equal to the number
	                               // of reported PTM sites
	
	private TIntDoubleHashMap modCoordMap = null; // holds modified amino acid positions
	
	private THashMap<String, Double> posPermutationScoreMap = null;
	private THashMap<String, Double> negPermutationScoreMap = null;

	// matched and unmatched peaks
	private ArrayList<Peak> posPeaks = null;
    private ArrayList<Peak> negPeaks = null;

	private SpectrumClass PeakList = null;
	
	
	private Peptide origPep;
	private Peptide score1pep; // top scoring peptide permutation
	private Peptide score2pep; // 2nd best scoring peptide permutation


	/**
	 * Default counstructor for {@link PSM}
	 */
	PSM() {
		origPep = new Peptide();
		modCoordMap = new TIntDoubleHashMap();
		isKeeper = false;
		useForModel = false;
		isDecoy = false;
		isUnambiguous = false;
        PeakList = new SpectrumClass();
	}

	
	void process() {
		
		origPep.initialize(modCoordMap);
		origPep.setCharge(this.charge);
		
		// Determine if the PSM is to be kept and if it should be used for modeling
		int keepingScore = 0;
		if(origPep.getNumPPS() > 0) keepingScore++;
		if(origPep.getNumRPS() > 0) keepingScore++;
		if(PSMscore >= Globals.scoreTH) keepingScore++;
        if(charge <= Globals.maxChargeState) keepingScore++;
        if(origPep.getPepLen() <= Globals.maxPepLen) keepingScore++;
		
		if(origPep.getNumPPS() == origPep.getNumRPS()) isUnambiguous = true;
		
		if(keepingScore == 5) isKeeper = true;
		
		if(isKeeper && (PSMscore >= Globals.modelTH)) useForModel = true;
		
		
		if(isKeeper) {
			// We will reconstruct the specId here to remove left-padded zeros
			String suffix = "";
			String scanStr = Integer.toString(scanNum);

			if(Globals.inputType == Constants.PEPXML) {
				if(Globals.spectrumSuffix.equalsIgnoreCase("mzXML")) suffix = "mzXML";
				if(Globals.spectrumSuffix.equalsIgnoreCase("mzML")) suffix = "mzML";
				if(Globals.spectrumSuffix.equalsIgnoreCase("mgf")) suffix = "mgf";

                /*
                ** Some scans might have a name formatted like this:
                ** expt_file_name.c.14033.14033.2
                ** we need to be able to handle these cases
                 */
                Pattern r = Pattern.compile("(.+)\\.\\d+\\.\\d+\\.\\d$"); // pattern of a spectrum ID
                Matcher m = r.matcher(specId);

                if(m.find()) {
                    srcFile = m.group(1) + "." + suffix;
                    specId = m.group(1) + "." + scanNum + "." + scanNum + "." +
							charge;
                }
                else {
                    log.info("\nERROR: PSM.java::Unable to correctly format specId variable.\n");
                    System.exit(0);
                }
            }
			else { 
				// We need to construct a valid TPP-specId variable for 
				// this PSM instance because it was read in from a tab-delimited file
                Pattern r = Pattern.compile("(.+)\\.(mgf|mzXML|mzML)$"); // pattern of a spectrum ID
                Matcher m = r.matcher(srcFile);

                if(m.find()) {
                    specId = m.group(1) + "." + scanNum + "." +
							scanNum + "." +
							charge;
                }
                else {
                    log.info("\nERROR: PSM.java::Unable to correctly format specId variable.\n");
                    System.exit(0);
                }
			}
			
			origPep.buildIonLadders();
		}
		
	}

	void recordSpectra(SpectrumClass S) {
		PeakList = S;
		// reduce the impact of the neutral loss peak (if any)
		if(Globals.reduceNL == 1) reduceNLpeak();
		// normalize the peak intensities to the median
		PeakList.medianNormalizeSpectra();
	}


	/**
	 * If the user asked for it, reduce the intensity of the precursor neutral loss peak
 	 */
	private void reduceNLpeak() {
		double pepMHplus = Globals
				.getFragmentIonMass(origPep.getModPeptide(), 1.0, Constants.WATER);

		// the Globals.precursorNLmass is a negative number
		double NLmass = MathFunctions
				.roundDouble( (pepMHplus + Globals.precursorNLmass), 3 );

		double NLmz = MathFunctions
				.roundDouble( (NLmass / (double) charge), 3 );
		
		// Do not perform this NL reduction if the peptide doesn't contain
		// a serine or theronine
		if(!origPep.getPeptide().contains("S") && !origPep.getPeptide().contains("T")) return;

		// Find the most intense peak in the spectrum
        double maxPk_mz = PeakList.mz[ PeakList.maxI_index ];

		double mzDelta = MathFunctions
				.roundDouble( Math.abs( (maxPk_mz - NLmz) ), 5);

		double mzErr = Constants.PROTON;
		if(mzDelta <= mzErr) { // the max peak is a neutral loss peak

            double orig_maxI = PeakList.maxI; // record the original intensity of the peak
            int orig_maxI_index = PeakList.maxI_index;

            PeakList.raw_intensity[ PeakList.maxI_index ] = 0; // bring the peak down to zero

            // identify the second-most intense peak in the spectrum
            double new_maxI = 0d;
            int new_maxI_index = 0;
            for(int i = 0; i < PeakList.N; i++) {
                if( PeakList.raw_intensity[i] > new_maxI ) {
                    new_maxI = PeakList.raw_intensity[i];
                    new_maxI_index = i;
                }
            }
            PeakList.maxI = new_maxI;

			// recompute the relative peak intensities
            PeakList.calcRelativeIntensity();
            PeakList.raw_intensity[ orig_maxI_index ] = orig_maxI;

            // set the peak intensity to be 50%
            PeakList.rel_intensity[ orig_maxI_index ] = 50.0;

//          PeakList.maxI *= 0.5; // reduce the peak intensity by 50%
		}

	}
	
	
	// Function returns a map of the modified residues in the given modPeptide
	private TIntDoubleHashMap getModMap(String modPep) {
		TIntDoubleHashMap ret = new TIntDoubleHashMap();
		
		// record n-terminal modification (if any)
		if(modCoordMap.containsKey(Constants.NTERM_MOD)) {
            ret.put(Constants.NTERM_MOD, Globals.ntermMass);
		}

		// record c-terminal modification (if any)
		if(modCoordMap.containsKey(Constants.CTERM_MOD)) {
			ret.put(Constants.CTERM_MOD, Globals.ntermMass);
		}
		
		// record all variable modifications
		int startAt = 0;
		
		for(int i = startAt; i < modPep.length(); i++) {
			double mass = 0;
			char c = modPep.charAt(i);
			
			// record decoy-modified residues
			if( Globals.isDecoyResidue(Character.toString(c)) ) {
				mass = Constants.AA_MASS_MAP.get(Character.toString(c));
				ret.put(i, mass);
			}
			// record normal-modified residues
			else if( Character.isLowerCase(c) ) {				
				mass = Constants.AA_MASS_MAP.get(Character.toString(c));
				ret.put(i, mass);
			}
		}
		
		
		// Search the original sequence assigned to this PSM to identify 
		// any residues that are already modified for this peptide. These
		// are residues that cannot be used for building decoys
		for(int p : origPep.getNonTargetMods().keySet()) {
			if( (p == Constants.NTERM_MOD) || (p == Constants.CTERM_MOD)) continue;
			
			String aa = origPep.getNonTargetMods().get(p);
			double mass = Constants.AA_MASS_MAP.get(aa);
			ret.put(p, mass);
		}
		
		return ret;
	}


	
	// This function identifies all of the peaks in PeakList that can be matched to any fragment ion of any permutation
	// of the peptide sequence associated with this PSM
	void matchAllPeaks() {

        posPeaks = new ArrayList(PeakList.N);

		TIntDoubleHashMap mcp = new TIntDoubleHashMap(); // use this to record modifications for each permutation

        // This variable will hold all of the fragment ions that are possible for every
        // permutation of the peptide sequence assigned to this peptide
        THashMap<String, Double> all_ions = new THashMap<>();

		for(String pep : posPermutationScoreMap.keySet()) {
            mcp.clear();
            mcp = getModMap(pep); // get the modified residues for this sequence permutation
            Peptide curPep = new Peptide();
            curPep.initialize(origPep.getPeptide(), pep, charge, mcp);
            curPep.buildIonLadders();

            // get the ion ladder for this permutation and match it against the
            // spectrum assigned to this PSM
            all_ions.putAll(curPep.getIonLadder());
        }


        // The same observed peak could be assigned to multiple theoretical peaks.
        // We use this data object to collect all of these cases for later processing.
        // k = observed Peak MZ value
        // v = array list of all the theoretical peaks that match to it
        TDoubleObjectHashMap<ArrayList<Peak>> candMatchedPks = new TDoubleObjectHashMap<>(PeakList.N);


        // Iterate over each theoretical ion
        for(String theo_ion : all_ions.keySet()) {

            double theo_mz = all_ions.get(theo_ion);

            Peak matchedPk = getMatchedPeak(theo_ion, theo_mz);
            if (null == matchedPk) continue; // no match was found for this theoretical peak

            if (candMatchedPks.containsKey(matchedPk.getMz())) {
                ArrayList<Peak> ary = candMatchedPks.get(matchedPk.getMz());
                Peak newEntry = new Peak(matchedPk);
                ary.add(newEntry);
                candMatchedPks.put(matchedPk.getMz(), ary);
            } else {
                ArrayList<Peak> ary = new ArrayList<>();
                ary.add(matchedPk);
                candMatchedPks.put(matchedPk.getMz(), ary);
            }
        }

        // For each matched peak in 'candMatchedPks' now select the match that is the closest by m/z value
        for(double mz : candMatchedPks.keys()) {
            ArrayList<Peak> ary = candMatchedPks.get(mz);
            ary.sort(Peak.comparator_mz_abs_dist);
            posPeaks.add(ary.get(0));
        }
        candMatchedPks.clear();
        candMatchedPks = null;


		// This can happen if there are no matched peaks 
		// identified for this spectrum due to the exclusion of neutral loss ions
		// or due to a very small fragment ion tolerance window.
		// When it happens we just exit the function.
		if(posPeaks.isEmpty()) {
			useForModel = false;
			posPeaks.clear();
			return;
		}

		// Now that you have matched all the peaks you can for this spectrum
		// All the remaining peaks are noise, we put them in the negPeaks arraylist
		negPeaks = new ArrayList<>(PeakList.N);

        for(int i = 0; i < PeakList.N; i++) {

            Peak negPk = PeakList.getPeakClassInstance(i);

            if( !posPeaks.contains(negPk) ) {

                double finalDist = 0;

                // record the distance of this unmatched peak to the nearest peak in 'posPeaks'
                ArrayList<Double> dist = new ArrayList( posPeaks.size() );
                ArrayList<Double> absDist = new ArrayList( posPeaks.size() );

                int j = 0;
                for(Peak matchedPk : posPeaks) {
                    double d = matchedPk.getMz() - negPk.getMz();
                    dist.add(d);
                    absDist.add( Math.abs(d) );
                    j++;
                }

                Collections.shuffle(dist);
                Collections.sort(absDist); // sort the m/z distances from low to high
                finalDist = absDist.get(0);
                for(double d : dist) {
                    if(Math.abs(d) == finalDist) {
                        negPk.setMatched(false);
                        negPk.setDist(d);
                        negPeaks.add( negPk );
                        break; // leave loop
                    }
                }
            }
        }
	}


    /*******
     * This function returns a Peak object matched to the given theoretical mz value.
     * The function returns a null object if no match was found.
     *
     * @param theo_ion
     * @param theo_mz
     * @return
     */
    private Peak getMatchedPeak(String theo_ion, double theo_mz) {

        Peak ret = null;

        double matchErr = 0;
        double a = 0, b = 0;

        // Compute the fragment error tolerance that will be used
        if(Globals.ms2tol_units == Constants.PPM_UNITS) {
            double ppmErr = Globals.ms2tol / Constants.PPM;
            matchErr = theo_mz * ppmErr;
        }
        else matchErr = Globals.ms2tol;

        matchErr *= 0.5; // split in half

        a = theo_mz - matchErr;
        b = theo_mz + matchErr;


        // Iterate over the peaks in PeakList
        ArrayList<Peak> candMatches = new ArrayList<>();
        for(int i = 0; i < PeakList.N; i++) {

            if( (PeakList.mz[i] >= a) && (PeakList.mz[i] <= b) ) {
                Peak pk = PeakList.getPeakClassInstance(i);
                candMatches.add(pk);
            }
        }

        // At least one match was found for this theoretical mz value
        if( !candMatches.isEmpty() ) {

            // identify the most intense peak in candMatches
            candMatches.sort(Peak.comparator_intensity_hi2low);

            ret = candMatches.get(0);
            ret.setMatched(true);
            ret.setDist(ret.getMz() - theo_mz); // obs - expected
            ret.setMatchedIonStr(theo_ion);
            ret.setMatchedIonMZ(theo_mz);
        }

        return  ret;
    }
	
	
	// Function creates all of the permutations for the sequence assigned to this PSM
	void generatePermutations(int RN) {
		posPermutationScoreMap = new THashMap();
		posPermutationScoreMap = origPep.getPermutations(0);
		
		if(!isUnambiguous) {
            negPermutationScoreMap = new THashMap();
            if( RN == 0 ) { // RN = runMode is selected to create decoys
                negPermutationScoreMap = origPep.getPermutations(1);
            }
		}
	}

    /**
     * Function to fill in required variables when a thread instance of this class is killed due to
     * either running over time or when the executor thread pool is shutdown.
     */
    private void killThreadResults() {
        deltaScore = -1;
        score1pep = origPep;
        score2pep = origPep;
    }
	
	
	// Function scores every candidate sequence associated with this PSM
	void scorePermutations() throws IOException {
		TIntDoubleHashMap mcp = null;
		
		File debugF = null;
		FileWriter fw = null;
		BufferedWriter bw = null;
		String line;

        if(Thread.interrupted()) { // this code just prevents excessive run times
            killThreadResults();
            return;
        }

		if(Globals.debugMode == Constants.WRITE_PERM_SCORES) {
				debugF = new File("all_scores.debug");
				
				if(!debugF.exists()) {
					fw = new FileWriter(debugF);
					bw = new BufferedWriter(fw);
					bw.write("specId\torigModPep\tcurPermutation\tisDecoy\tscore\n");
				}
				else {
					fw = new FileWriter(debugF, true); // open for appending
					bw = new BufferedWriter(fw);
				}
		}

        if(Globals.debugMode == Constants.WRITE_ALL_MATCHED_PK_SCORES) {
            debugF = new File("all_matched_pks.debug");

            if(!debugF.exists()) {
                fw = new FileWriter(debugF);
                bw = new BufferedWriter(fw);
                bw.write("specId\tcurPermutation\tisDecoy\tionSeq\tmz\tintensity\tmzDist\tIscore\tDscore\tscore\n");
            }
            else {
                fw = new FileWriter(debugF, true); // open for appending
                bw = new BufferedWriter(fw);
            }
        }


        // Scores for forward sequences
		for(String curSeq : posPermutationScoreMap.keySet()) {

            mcp = this.getModMap(curSeq);
			
			Peptide curPep = new Peptide();
			curPep.initialize(origPep.getPeptide(), curSeq, charge, mcp);
			curPep.buildIonLadders();
			curPep.setNumPPS(origPep.getNumPPS());
			curPep.setNumRPS(origPep.getNumRPS());

            curPep.matchPeaks(PeakList); // match all the peaks you can for this peptide permutation

			if(Globals.scoringAlgorithm == Constants.CID) curPep.calcScoreCID();
			if(Globals.scoringAlgorithm == Constants.HCD) curPep.calcScoreHCD();
			
			if(Globals.debugMode == Constants.WRITE_PERM_SCORES) {
				line = specId + "\t" + origPep.getModPeptide() + "\t" +
					   curPep.getModPeptide() + "\t0\t" + curPep.getScore() + "\n";
				bw.write(line);
			}

            if(Globals.debugMode == Constants.WRITE_ALL_MATCHED_PK_SCORES) {
                for(Peak pk : curPep.getMatchedPeaks() ) {
                    line = specId + "\t" + curSeq + "\t" + Globals.isDecoySeq(curSeq) + "\t" +
                           pk.getMatchedIonStr() + "\t" + pk.getMz() + "\t" +
                           pk.getRelIntensity() + "\t" + pk.getDist() + "\t" +
                           pk.getIntensityScore() + "\t" + pk.getDistScore() + "\t" + pk.getScore() + "\n";
                    bw.write(line);
                }
                bw.write("\n");
            }

			posPermutationScoreMap.put(curSeq, curPep.getScore());
			curPep = null;
			mcp = null;
		}

		if(!isUnambiguous) { // at least 2 possible permutations exist for this PSM

			// Scores for decoy sequences
			for(String curSeq : negPermutationScoreMap.keySet()) {
				mcp = this.getModMap(curSeq);

				Peptide curPep = new Peptide();
				mcp = this.getModMap(curSeq);
				curPep.initialize(origPep.getModPeptide(), curSeq, charge, mcp);
				curPep.buildIonLadders();
				curPep.setNumPPS(origPep.getNumPPS());
				curPep.setNumPPS(origPep.getNumRPS());

                curPep.matchPeaks(PeakList);

				if(Globals.scoringAlgorithm == Constants.CID) curPep.calcScoreCID();
				if(Globals.scoringAlgorithm == Constants.HCD) curPep.calcScoreHCD();

				if(Globals.debugMode == Constants.WRITE_PERM_SCORES) {
					line = specId + "\t" + origPep.getModPeptide() + "\t" +
						curPep.getModPeptide() + "\t1\t" + curPep.getScore() + "\n";
					bw.write(line);
				}

                if(Globals.debugMode == Constants.WRITE_ALL_MATCHED_PK_SCORES) {
                    for(Peak pk : curPep.getMatchedPeaks() ) {
                        line = specId + "\t" + curSeq + "\t" + Globals.isDecoySeq(curSeq) + "\t" +
                                pk.getMatchedIonStr() + "\t" + pk.getMz() + "\t" +
                                pk.getRelIntensity() + "\t" + pk.getDist() + "\t" +
                                pk.getIntensityScore() + "\t" + pk.getDistScore() + "\t" + pk.getScore() + "\n";
                        bw.write(line);
                    }
                    bw.write("\n");
                }


				negPermutationScoreMap.put(curSeq, curPep.getScore());
				curPep = null;
				mcp = null;
			}
		}
		
		if(Globals.debugMode == Constants.WRITE_PERM_SCORES) bw.close();
        if(Globals.debugMode == Constants.WRITE_ALL_MATCHED_PK_SCORES) bw.close();


        // collect all of the scores into a single arraylist and select the top
		// two matches to compute the delta score
		ArrayList<Double> allScores = new ArrayList();
		for(String curSeq : posPermutationScoreMap.keySet()) {
			double d = posPermutationScoreMap.get(curSeq);
			allScores.add(d);
		}
		
		if(!isUnambiguous) {
			for(String curSeq : negPermutationScoreMap.keySet()) {
				double d = negPermutationScoreMap.get(curSeq);
				allScores.add(d);
			}
		}
		
		Collections.sort(allScores); // low to high
		Collections.reverse(allScores);// high to low
		
		double score1 = MathFunctions.roundDouble(allScores.get(0), 6);
		double score2 = 0;
		if(!isUnambiguous) score2 = MathFunctions.roundDouble(allScores.get(1), 6);
		
		String pep1 = "";
		String pep2 = "";
		int numAssigned = 0;
		
		if(!isUnambiguous) {
			// Find the permutations that have the values stored in score1 and score2
			// We'll first try to find the top scores among the non-decoy 
			// permutations
			for(Entry<String, Double> e : posPermutationScoreMap.entrySet()) {
				String curSeq = e.getKey();
				double x = e.getValue();
				double d = MathFunctions.roundDouble(x, 6);
				
				if( (d == score1) && (pep1.isEmpty()) ) { 
					pep1 = curSeq;
					numAssigned++;
				}
				else if( (d == score2) && (pep2.isEmpty()) ) {
					pep2 = curSeq;
					numAssigned++;
				}
				
				if(numAssigned == 2) break;
			}
			
			
			// if this is true, then you need to search among the decoys
			if(numAssigned != 2) {
				for(Entry<String, Double> e : negPermutationScoreMap.entrySet()) {
					String curSeq = e.getKey();
					double x = e.getValue();
					double d = MathFunctions.roundDouble(x, 6);
					
					if( (d == score1) && (pep1.isEmpty()) ) { 
						pep1 = curSeq;
						numAssigned++;
					}
					else if( (d == score2) && (pep2.isEmpty()) ) {
						pep2 = curSeq;
						numAssigned++;
					}

					if(numAssigned == 2) break;
				}
			}
		}
		else { // special case for unambiguous PSMs
			for(Entry<String, Double> e : posPermutationScoreMap.entrySet()) {
				String curSeq = e.getKey();
				pep1 = curSeq;
			}
		}
		
		
		// Record the best match into score1
		mcp = this.getModMap(pep1);
		score1pep = new Peptide();
		score1pep.initialize(origPep.getPeptide(), pep1, charge, mcp);
		score1pep.setScore(score1);
		isDecoy = score1pep.isDecoyPep(); // determine if this PSM is a decoy hit
		mcp = null;
		
		if(isUnambiguous) {
			// This is an unambiguous case.
			// Make a duplicate of score1pep and assign it to score2pep.
			//mcp = new HashMap();
			mcp = this.getModMap(pep1);
			score2pep = new Peptide();
			score2pep.initialize(origPep.getPeptide(), pep1, charge, mcp);
			score2pep.setScore(0);
			mcp = null;
			deltaScore = score1;
		}
		else {
			// Record second best match
			//mcp = new HashMap();
			mcp = this.getModMap(pep2);
			score2pep = new Peptide();
			score2pep.initialize(origPep.getPeptide(), pep2, charge, mcp);
			score2pep.setScore(score2);
			mcp = null;
			deltaScore = score1 - score2;
		}
	}

	
	// Function returns results as a string
	String getResults() {

		DecimalFormat df = new DecimalFormat("#.####");
		
		String ret = specId + "\t";

		if(Globals.peptideRepresentation == Constants.SINGLE_CHAR) {
			ret += origPep.getPeptide() + "\t";
			ret += score1pep.getModPeptide() + "\t";
			ret += score2pep.getModPeptide() + "\t";
			
		}
		else {
			ret += origPep.getModPepTPP() + "\t";
			ret += score1pep.getModPepTPP() + "\t";
			ret += score2pep.getModPepTPP() + "\t";
		}
		
		ret += origPep.getNumPPS() + "\t" + origPep.getNumRPS() + "\t";
		
		ret += df.format(PSMscore) + "\t";

        if(Globals.runMode == Constants.REPORT_DECOYS)
		    ret += (score1pep.isDecoyPep() ? 1 : 0) + "\t" + (score2pep.isDecoyPep() ? 1 : 0) + "\t";
		
		ret += df.format(deltaScore) + "\t";
		
		ret += df.format(score1pep.getScore()) + "\t"
			 + df.format(score2pep.getScore()) + "\t";

		ret += df.format(globalFDR) + "\t" + df.format(localFDR);
		ret += "\n";
		
		return ret;
	}


	/**
	 * This function will write the scored peaks for the top permutation assigned
	 * to this PSM. The data is written in a TAB-delimited format.
	 * @throws IOException
	 */
	void debugWriteScoredPeaks() throws IOException {

		DecimalFormat df = new DecimalFormat("#.#####");
		int numDecimals = 4;
		
		String tsvFileName = "";
		if(Globals.inputType == Constants.PEPXML) tsvFileName = specId + ".tsv";
		else {
			int i = srcFile.lastIndexOf(Globals.spectrumSuffix);
			String k = String.format("%05d", scanNum) + Integer.toString(charge);
			tsvFileName = srcFile.substring(0, i) + k + ".tsv";
		}
		
		// With the next 4 lines you create output directory if necessary
		String outDirName = "debug_scored_peaks." + Globals.dateStamp;
		String outFileStr = outDirName + "/" + tsvFileName;
		
		File outF = new File(outFileStr);
		outF.getParentFile().mkdirs();
		
		FileWriter fw = new FileWriter(outF.getAbsoluteFile());
		BufferedWriter bw = new BufferedWriter(fw);
		
		bw.write("num\tmodPeptide\tfragmentIon\tobs m/z\tmatchDist\trelIntensity\t" +
				 "normIntensity\tDscore\tIscore\tfinalScore\n");

		score1pep.buildIonLadders();
        score1pep.matchPeaks(PeakList);

		if(Globals.scoringAlgorithm == Constants.CID) score1pep.calcScoreCID();
		if(Globals.scoringAlgorithm == Constants.HCD) score1pep.calcScoreHCD();
		
		ArrayList<Peak> mPK = score1pep.getMatchedPeaks();
		mPK.sort(Peak.comparator_mz);
		for(Peak pk : mPK) {
			if(!pk.isMatched()) continue;
			
			String mz = df.format(MathFunctions.roundDouble(pk.getMz(), 4));
			String relI = df.format(MathFunctions.roundDouble(pk.getRelIntensity(), 4));
			String normI = df.format(MathFunctions.roundDouble(pk.getNormIntensity(), 4));
			String Iscore = df.format(MathFunctions.roundDouble(pk.getIntensityScore(), 4));
			String Dscore = df.format(MathFunctions.roundDouble(pk.getDistScore(), 4));
			String score = df.format(MathFunctions.roundDouble(pk.getScore(), 4));
			String dist = df.format(MathFunctions.roundDouble(pk.getDist(), 4));
			
			bw.write("0\t" + score1pep.getModPeptide() + "\t" + pk.getMatchedIonStr() + "\t" + mz + "\t" + dist + "\t" +
					 relI + "\t" + normI + "\t" + Dscore + "\t" + Iscore + "\t" +
					 score + "\n");
		}
		

		score2pep.buildIonLadders();
		score2pep.matchPeaks(PeakList);
        if(Globals.scoringAlgorithm == Constants.CID) score2pep.calcScoreCID();
		if(Globals.scoringAlgorithm == Constants.HCD) score2pep.calcScoreHCD();

		
		mPK.clear();
		mPK = score2pep.getMatchedPeaks();
		mPK.sort(Peak.comparator_mz);
		for(Peak pk : mPK) {
			if(!pk.isMatched()) continue;
			
			String mz = df.format(MathFunctions.roundDouble(pk.getMz(), 4));
			String relI = df.format(MathFunctions.roundDouble(pk.getRelIntensity(), 4));
			String normI = df.format(MathFunctions.roundDouble(pk.getNormIntensity(), 4));
			String Iscore = df.format(MathFunctions.roundDouble(pk.getIntensityScore(), 4));
			String Dscore = df.format(MathFunctions.roundDouble(pk.getDistScore(), 4));
			String score = df.format(MathFunctions.roundDouble(pk.getScore(), 4));
			String dist = df.format(MathFunctions.roundDouble(pk.getDist(), 4));
			
			bw.write("1\t" + score2pep.getModPeptide() + "\t" + pk.getMatchedIonStr() + "\t" + mz + "\t" + dist + "\t" +
					 relI + "\t" + normI + "\t" + Dscore + "\t" + Iscore + "\t" +
					 score + "\n");
		}
		
		
		bw.close();
	}


	/**
	 * This function returns the matched peaks for the top 2 permutations assigned
	 * to this PSM. The data is returned in a TAB-delimited format.
	 * @return Line String
	 */
	String writeMatchedPks() {

		DecimalFormat df = new DecimalFormat("#.####");
		int numDecimals = 4;

        String ret = "";

        score1pep.buildIonLadders();
        score1pep.matchPeaks(PeakList);

        if(Globals.scoringAlgorithm == Constants.CID) score1pep.calcScoreCID();
        if(Globals.scoringAlgorithm == Constants.HCD) score1pep.calcScoreHCD();

        ArrayList<Peak> mPK = score1pep.getMatchedPeaks();
        mPK.sort(Peak.comparator_mz);
        for(Peak pk : mPK) {
            if(!pk.isMatched()) continue;

            String mz = df.format(MathFunctions.roundDouble(pk.getMz(), numDecimals));
            String relI = df.format(MathFunctions.roundDouble(pk.getRelIntensity(), numDecimals));
            String Iscore = df.format(MathFunctions.roundDouble(pk.getIntensityScore(), numDecimals));
            String Dscore = df.format(MathFunctions.roundDouble(pk.getDistScore(), numDecimals));
            String score = df.format(MathFunctions.roundDouble(pk.getDistScore(), numDecimals));

            ret += specId + "\t1\t" + score1pep.getModPeptide() + "\t" + pk.getMatchedIonStr() + "\t" +
                   mz + "\t" + relI + "\t" + Dscore + "\t" + Iscore + "\t" +
                   score + "\n";
        }


        score2pep.buildIonLadders();
        score2pep.matchPeaks(PeakList);
        if(Globals.scoringAlgorithm == Constants.CID)
        	score2pep.calcScoreCID();
        if(Globals.scoringAlgorithm == Constants.HCD)
        	score2pep.calcScoreHCD();


        mPK.clear();
        mPK = score2pep.getMatchedPeaks();
        mPK.sort(Peak.comparator_mz);
        for(Peak pk : mPK) {
            if(!pk.isMatched()) continue;

            String mz = df.format(MathFunctions.roundDouble(pk.getMz(), 4));
            String relI = df.format(MathFunctions.roundDouble(pk.getRelIntensity(), 4));
            String Iscore = df.format(MathFunctions.roundDouble(pk.getIntensityScore(), 4));
            String Dscore = df.format(MathFunctions.roundDouble(pk.getDistScore(), 4));
            String score = df.format(MathFunctions.roundDouble(pk.getScore(), 4));

            ret += specId + "\t2\t" + score2pep.getModPeptide() + "\t" + pk.getMatchedIonStr() + "\t" +
                    mz + "\t" + relI + "\t" + Dscore + "\t" + Iscore + "\t" +
                    score + "\n";
        }

        return(ret);
    }

	/**
	 * This function prepares the PSM to be scored again *AFTER* the FLR has been estimated.
	 */
	public void clearScores() {

        deltaScore = 0;
        globalFDR = Double.NaN;
        localFDR = Double.NaN;
        isDecoy = false;
        score1pep = null;
        score2pep = null;
        posPermutationScoreMap = null;
        negPermutationScoreMap = null;

    }

	public String getSpecId() {
		return specId;
	}

	public String getSrcFile() {
		return srcFile;
	}

	public int getScanNum() {
		return scanNum;
	}

	public int getCharge() {
		return charge;
	}

	public double getPSMscore() {
		return PSMscore;
	}

	public double getDeltaScore() {
		return deltaScore;
	}

	public double getLocalFDR() {
		return localFDR;
	}

	public double getGlobalFDR() {
		return globalFDR;
	}

	public void setLocalFDR(double localFDR) {
		this.localFDR = localFDR;
	}

	public void setGlobalFDR(double globalFDR) {
		this.globalFDR = globalFDR;
	}

	public boolean isKeeper() {
		return isKeeper;
	}

	public boolean isUseForModel() {
		return useForModel;
	}

	public boolean isDecoy() {
		return isDecoy;
	}

	public boolean isUnambiguous() {
		return isUnambiguous;
	}

	public TIntDoubleHashMap getModCoordMap() {
		return modCoordMap;
	}

	public THashMap<String, Double> getPosPermutationScoreMap() {
		return posPermutationScoreMap;
	}

	public THashMap<String, Double> getNegPermutationScoreMap() {
		return negPermutationScoreMap;
	}

	public ArrayList<Peak> getPosPeaks() {
		return posPeaks;
	}

	public ArrayList<Peak> getNegPeaks() {
		return negPeaks;
	}

	public SpectrumClass getPeakList() {
		return PeakList;
	}

	public Peptide getOrigPep() {
		return origPep;
	}

	public String getPeptideSequence(){
		return origPep.getPeptide();
	}

	public Peptide getScore1pep() {
		return score1pep;
	}

	public Peptide getScore2pep() {
		return score2pep;
	}

	public void setSpecId(String specId) {
		this.specId = specId;
	}

	public void setSrcFile(String srcFile) {
		this.srcFile = srcFile;
	}

	public void setScanNum(int scanNum) {
		this.scanNum = scanNum;
	}

	public void setCharge(int charge) {
		this.charge = charge;
	}

	public void setPSMscore(double PSMscore) {
		this.PSMscore = PSMscore;
	}

	public void setDeltaScore(double deltaScore) {
		this.deltaScore = deltaScore;
	}

	public void setKeeper(boolean keeper) {
		isKeeper = keeper;
	}

	public void setUseForModel(boolean useForModel) {
		this.useForModel = useForModel;
	}

	public void setDecoy(boolean decoy) {
		isDecoy = decoy;
	}

	public void setUnambiguous(boolean unambiguous) {
		isUnambiguous = unambiguous;
	}

	public void setPeakList(SpectrumClass peakList) {
		PeakList = peakList;
	}

	public void setOrigPep(Peptide origPep) {
		this.origPep = origPep;
	}

	public void setPeptideSequence(String sequence){
		this.origPep.setPeptide(sequence);
	}

	public void setScore1pep(Peptide score1pep) {
		this.score1pep = score1pep;
	}

	public void setScore2pep(Peptide score2pep) {
		this.score2pep = score2pep;
	}
}
