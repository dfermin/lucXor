/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor.common;

import gnu.trove.map.TMap;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import lombok.Builder;
import lombok.extern.slf4j.Slf4j;
import lucxor.algorithm.ModelDataCID;
import lucxor.algorithm.ModelDataHCD;
import lucxor.utils.Constants;
import lucxor.LucXorConfiguration;
import lucxor.utils.MathFunctions;
import lucxor.utils.Utils;
import org.slf4j.LoggerFactory;

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
@Builder
@Slf4j
public class PSM {


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
	private TIntDoubleHashMap modCoordMap; // holds modified amino acid positions
	
	private TMap<String, Double> posPermutationScoreMap = null;
	private TMap<String, Double> negPermutationScoreMap = null;

	// matched and unmatched peaks
	private List<Peak> posPeaks = null;
    private List<Peak> negPeaks = null;

	private Spectrum PeakList;

	private Peptide origPep;
	private Peptide score1pep; // top scoring peptide permutation
	private Peptide score2pep; // 2nd best scoring peptide permutation


	/**
	 * Default constructor for {@link PSM}
	 */
	public PSM() {
		origPep = new Peptide();
		modCoordMap = new TIntDoubleHashMap();
		isKeeper = false;
		useForModel = false;
		isDecoy = false;
		isUnambiguous = false;
        PeakList = new Spectrum();
	}

	public PSM(String specId, String srcFile, int scanNum, int charge, double PSMscore, double deltaScore,
			   double localFDR, double globalFDR, boolean isKeeper, boolean useForModel, boolean isDecoy,
			   boolean isUnambiguous, TIntDoubleHashMap modCoordMap, TMap<String, Double> posPermutationScoreMap,
			   TMap<String, Double> negPermutationScoreMap, List<Peak> posPeaks, List<Peak> negPeaks,
			   Spectrum peakList, Peptide origPep, Peptide score1pep, Peptide score2pep) {
		this.specId = specId;
		this.srcFile = srcFile;
		this.scanNum = scanNum;
		this.charge = charge;
		this.PSMscore = PSMscore;
		this.deltaScore = deltaScore;
		this.localFDR = localFDR;
		this.globalFDR = globalFDR;
		this.isKeeper = isKeeper;
		this.useForModel = useForModel;
		this.isDecoy = isDecoy;
		this.isUnambiguous = isUnambiguous;
		this.modCoordMap = modCoordMap;
		this.posPermutationScoreMap = posPermutationScoreMap;
		this.negPermutationScoreMap = negPermutationScoreMap;
		this.posPeaks = posPeaks;
		this.negPeaks = negPeaks;
		PeakList = peakList;
		this.origPep = origPep;
		this.score1pep = score1pep;
		this.score2pep = score2pep;
	}

	public void process() {

		origPep.initialize(modCoordMap);
		origPep.setCharge(this.charge);
		
		// Determine if the PSM is to be kept and if it should be used for modeling
		int keepingScore = 0;
		if(origPep.getNumPPS() > 0) keepingScore++;
		if(origPep.getNumRPS() > 0) keepingScore++;
		if(PSMscore >= LucXorConfiguration.getSCORETH()) keepingScore++;
        if(charge <= LucXorConfiguration.getMaxChargeState()) keepingScore++;
        if(origPep.getPepLen() <= LucXorConfiguration.getMaxPepLength()) keepingScore++;
		
		if(origPep.getNumPPS() == origPep.getNumRPS()) isUnambiguous = true;
		
		if(keepingScore == 5) isKeeper = true;
		
		if(isKeeper && (PSMscore >= LucXorConfiguration.getMODELTH()))
			useForModel = true;
		
		
		if(isKeeper) {
			// We will reconstruct the specId here to remove left-padded zeros
			String suffix = "";

			if(LucXorConfiguration.getInputType() == Constants.PEPXML) {
				if(LucXorConfiguration.getSpectrumPrefix().equalsIgnoreCase(Constants.MZXML_TYPE))
					suffix = Constants.MZXML_TYPE;
				if(LucXorConfiguration.getSpectrumPrefix().equalsIgnoreCase(Constants.MZML_TYPE))
					suffix = Constants.MZML_TYPE;
				if(LucXorConfiguration.getSpectrumPrefix().equalsIgnoreCase(Constants.MGF_TYPE))
					suffix = Constants.MGF_TYPE;
				if(LucXorConfiguration.getSpectrumPrefix().equalsIgnoreCase(Constants.PRIDE_TYPE))
					suffix = Constants.PRIDE_TYPE;

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

	public void recordSpectra(Spectrum S) {
		PeakList = S;
		// reduce the impact of the neutral loss peak (if any)
		if(LucXorConfiguration.getREDUCENL() == 1) reduceNLpeak();
		// normalize the peak intensities to the median
		PeakList.medianNormalizeSpectra();
	}


	/**
	 * If the user asked for it, reduce the intensity of the precursor neutral loss peak
 	 */
	private void reduceNLpeak() {
		double pepMHplus = Utils.getFragmentIonMass(origPep.getModPeptide(), 1.0, Constants.WATER);

		// the LucXorConfiguration.precursorNLmass is a negative number
		double NLmass = MathFunctions
				.roundDouble( (pepMHplus + LucXorConfiguration.getPrecursorNlMass()), 3 );

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
            for(int i = 0; i < PeakList.N; i++) {
                if( PeakList.raw_intensity[i] > new_maxI ) {
                    new_maxI = PeakList.raw_intensity[i];
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
            ret.put(Constants.NTERM_MOD, LucXorConfiguration.getNtermMass());
		}

		// record c-terminal modification (if any)
		if(modCoordMap.containsKey(Constants.CTERM_MOD)) {
			ret.put(Constants.CTERM_MOD, LucXorConfiguration.getNtermMass());
		}
		
		// record all variable modifications
		int startAt = 0;
		
		for(int i = startAt; i < modPep.length(); i++) {
			double mass;
			char c = modPep.charAt(i);
			
			// record decoy-modified residues
			if( Utils.isDecoyResidue(Character.toString(c)) ) {
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
	public void matchAllPeaks() {

        posPeaks = new ArrayList<>(PeakList.N);

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
        for(Entry<String, Double> stringDoubleEntry : all_ions.entrySet()) {

            double theo_mz = stringDoubleEntry.getValue();

            Peak matchedPk = getMatchedPeak(stringDoubleEntry.getKey(), theo_mz);
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

                double finalDist;

                // record the distance of this unmatched peak to the nearest peak in 'posPeaks'
                ArrayList<Double> dist = new ArrayList<>( posPeaks.size() );
                ArrayList<Double> absDist = new ArrayList<>( posPeaks.size() );

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

        double matchErr;
        double a, b;

        // Compute the fragment error tolerance that will be used
        if(LucXorConfiguration.getMs2tolUnits() == Constants.PPM_UNITS) {
            double ppmErr = LucXorConfiguration.getMS2TOL() / Constants.PPM;
            matchErr = theo_mz * ppmErr;
        }
        else matchErr = LucXorConfiguration.getMS2TOL();

        matchErr *= 0.5; // split in half

        a = theo_mz - matchErr;
        b = theo_mz + matchErr;


        // Iterate over the peaks in PeakList
        ArrayList<Peak> candMatches = new ArrayList<>(PeakList.N/2);
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


	/**
	 * Function creates all of the permutations for the sequence assigned to this PSM
	 * @param RN iterations
	 */
	public void generatePermutations(int RN) {
		posPermutationScoreMap = new THashMap<>();
		posPermutationScoreMap = origPep.getPermutations(0);
		
		if(!isUnambiguous) {
            negPermutationScoreMap = new THashMap<>();
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


	/**
	 * Function scores every candidate sequence associated with this PSM
	 * @throws IOException
	 * @param modelingMapHCD
	 * @param modelingMapCID
	 */
	public void scorePermutations(TMap<Integer, ModelDataHCD> modelingMapHCD,
								  TMap<Integer, ModelDataCID> modelingMapCID) throws IOException {
		TIntDoubleHashMap mcp;
		
		File debugF;
		FileWriter fw;
		BufferedWriter bw = null;
		String line;

        if(Thread.interrupted()) { // this code just prevents excessive run times
            killThreadResults();
            return;
        }

		if(LucXorConfiguration.getDebugMode() == Constants.WRITE_PERM_SCORES) {
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

        if(LucXorConfiguration.getDebugMode() == Constants.WRITE_ALL_MATCHED_PK_SCORES) {
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
		for(Entry<String, Double> curSeq : posPermutationScoreMap.entrySet()) {

            mcp = this.getModMap(curSeq.getKey());
			
			Peptide curPep = new Peptide();
			curPep.initialize(origPep.getPeptide(), curSeq.getKey(), charge, mcp);
			curPep.buildIonLadders();
			curPep.setNumPPS(origPep.getNumPPS());
			curPep.setNumRPS(origPep.getNumRPS());

            curPep.matchPeaks(PeakList); // match all the peaks you can for this peptide permutation

			if(LucXorConfiguration.getScoringAlgorithm() == Constants.CID)
				curPep.calcScoreCID(modelingMapCID);
			if(LucXorConfiguration.getScoringAlgorithm() == Constants.HCD)
				curPep.calcScoreHCD(modelingMapHCD);
			
			if(LucXorConfiguration.getDebugMode() == Constants.WRITE_PERM_SCORES) {
				line = specId + "\t" + origPep.getModPeptide() + "\t" +
					   curPep.getModPeptide() + "\t0\t" + curPep.getScore() + "\n";
				bw.write(line);
			}

            if(LucXorConfiguration.getDebugMode() == Constants.WRITE_ALL_MATCHED_PK_SCORES) {
                for(Peak pk : curPep.getMatchedPeaks() ) {
                    line = specId + "\t" + curSeq + "\t" + Utils.isDecoySeq(curSeq.getKey()) + "\t" +
                           pk.getMatchedIonStr() + "\t" + pk.getMz() + "\t" +
                           pk.getRelIntensity() + "\t" + pk.getDist() + "\t" +
                           pk.getIntensityScore() + "\t" + pk.getDistScore() + "\t" + pk.getScore() + "\n";
                    bw.write(line);
                }
                bw.write("\n");
            }
			posPermutationScoreMap.put(curSeq.getKey(), curPep.getScore());
		}

		if(!isUnambiguous) { // at least 2 possible permutations exist for this PSM

			// Scores for decoy sequences
			for(Entry<String, Double> curSeq : negPermutationScoreMap.entrySet()) {

				Peptide curPep = new Peptide();
				mcp = this.getModMap(curSeq.getKey());
				curPep.initialize(origPep.getModPeptide(), curSeq.getKey(), charge, mcp);
				curPep.buildIonLadders();
				curPep.setNumPPS(origPep.getNumPPS());
				curPep.setNumPPS(origPep.getNumRPS());

                curPep.matchPeaks(PeakList);

				if(LucXorConfiguration.getScoringAlgorithm() == Constants.CID)
					curPep.calcScoreCID(modelingMapCID);
				if(LucXorConfiguration.getScoringAlgorithm() == Constants.HCD)
					curPep.calcScoreHCD(modelingMapHCD);

				if(LucXorConfiguration.getDebugMode() == Constants.WRITE_PERM_SCORES) {
					line = specId + "\t" + origPep.getModPeptide() + "\t" +
						curPep.getModPeptide() + "\t1\t" + curPep.getScore() + "\n";
					bw.write(line);
				}

                if(LucXorConfiguration.getDebugMode() == Constants.WRITE_ALL_MATCHED_PK_SCORES) {
                    for(Peak pk : curPep.getMatchedPeaks() ) {
                        line = specId + "\t" + curSeq + "\t" + Utils.isDecoySeq(curSeq.getKey()) + "\t" +
                                pk.getMatchedIonStr() + "\t" + pk.getMz() + "\t" +
                                pk.getRelIntensity() + "\t" + pk.getDist() + "\t" +
                                pk.getIntensityScore() + "\t" + pk.getDistScore() + "\t" + pk.getScore() + "\n";
                        bw.write(line);
                    }
                    bw.write("\n");
                }
				negPermutationScoreMap.put(curSeq.getKey(), curPep.getScore());
			}
		}
		
		if(LucXorConfiguration.getDebugMode() == Constants.WRITE_PERM_SCORES)
			bw.close();
        if(LucXorConfiguration.getDebugMode() == Constants.WRITE_ALL_MATCHED_PK_SCORES)
        	bw.close();


        // collect all of the scores into a single arraylist and select the top
		// two matches to compute the delta score
		List<Double> allScores = new ArrayList<>(posPermutationScoreMap.size());
		for(Entry<String, Double> stringDoubleEntry : posPermutationScoreMap.entrySet()) {
			double d = stringDoubleEntry.getValue();
			allScores.add(d);
		}
		
		if(!isUnambiguous) {
			for(Entry<String, Double> stringDoubleEntry : negPermutationScoreMap.entrySet()) {
				double d = stringDoubleEntry.getValue();
				allScores.add(d);
			}
		}
		
		Collections.sort(allScores); // low to high
		Collections.reverse(allScores);// high to low

		double score1 = MathFunctions.roundDouble(allScores.get(0), 6);
		double score2 = 0;
		if(!isUnambiguous)
			score2 = MathFunctions.roundDouble(allScores.get(1), 6);
		
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
				pep1 = e.getKey();
			}
		}
		
		
		// Record the best match into score1
		mcp = this.getModMap(pep1);
		score1pep = new Peptide();
		score1pep.initialize(origPep.getPeptide(), pep1, charge, mcp);
		score1pep.setScore(score1);
		isDecoy = score1pep.isDecoyPep(); // determine if this PSM is a decoy hit

		if(isUnambiguous) {
			// This is an unambiguous case.
			// Make a duplicate of score1pep and assign it to score2pep.
			//mcp = new HashMap();
			mcp = this.getModMap(pep1);
			score2pep = new Peptide();
			score2pep.initialize(origPep.getPeptide(), pep1, charge, mcp);
			score2pep.setScore(0);
			deltaScore = score1;
		}
		else {
			// Record second best match
			//mcp = new HashMap();
			mcp = this.getModMap(pep2);
			score2pep = new Peptide();
			score2pep.initialize(origPep.getPeptide(), pep2, charge, mcp);
			score2pep.setScore(score2);
			deltaScore = score1 - score2;
		}
	}


	/**
	 * Function returns results as a string
	 * @return
	 */
	public String getResults() {

		DecimalFormat df = new DecimalFormat("#.####");
		
		String ret = specId + "\t";

		if(LucXorConfiguration.getPeptideRepresentation() == Constants.SINGLE_CHAR) {
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

        if(LucXorConfiguration.getRunMode() == Constants.REPORT_DECOYS)
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
	public void debugWriteScoredPeaks(String dateStamp, TMap<Integer, ModelDataHCD> modelingMapHCD, TMap<Integer, ModelDataCID> modelingMapCID) throws IOException {

		DecimalFormat df = new DecimalFormat("#.#####");
		int numDecimals = 4;
		
		String tsvFileName;
		if(LucXorConfiguration.getInputType() == Constants.PEPXML) tsvFileName = specId + ".tsv";
		else {
			int i = srcFile.lastIndexOf(LucXorConfiguration.getSpectrumPrefix());
			String k = String.format("%05d", scanNum) + charge;
			tsvFileName = srcFile.substring(0, i) + k + ".tsv";
		}
		
		// With the next 4 lines you create output directory if necessary
		String outDirName = "debug_scored_peaks." + dateStamp;
		String outFileStr = outDirName + "/" + tsvFileName;
		
		File outF = new File(outFileStr);
		outF.getParentFile().mkdirs();
		
		FileWriter fw = new FileWriter(outF.getAbsoluteFile());
		BufferedWriter bw = new BufferedWriter(fw);
		
		bw.write("num\tmodPeptide\tfragmentIon\tobs m/z\tmatchDist\trelIntensity\t" +
				 "normIntensity\tDscore\tIscore\tfinalScore\n");

		score1pep.buildIonLadders();
        score1pep.matchPeaks(PeakList);

		if(LucXorConfiguration.getScoringAlgorithm() == Constants.CID) score1pep.calcScoreCID(modelingMapCID);
		if(LucXorConfiguration.getScoringAlgorithm() == Constants.HCD) score1pep.calcScoreHCD(modelingMapHCD);
		
		List<Peak> mPK = score1pep.getMatchedPeaks();
		mPK.sort(Peak.comparator_mz);
		for(Peak pk : mPK) {
			if(!pk.isMatched()) continue;
			
			String mz = df.format(MathFunctions.roundDouble(pk.getMz(), numDecimals));
			String relI = df.format(MathFunctions.roundDouble(pk.getRelIntensity(), numDecimals));
			String normI = df.format(MathFunctions.roundDouble(pk.getNormIntensity(), numDecimals));
			String Iscore = df.format(MathFunctions.roundDouble(pk.getIntensityScore(), numDecimals));
			String Dscore = df.format(MathFunctions.roundDouble(pk.getDistScore(), numDecimals));
			String score = df.format(MathFunctions.roundDouble(pk.getScore(), numDecimals));
			String dist = df.format(MathFunctions.roundDouble(pk.getDist(), numDecimals));
			
			bw.write("0\t" + score1pep.getModPeptide() + "\t" + pk.getMatchedIonStr() + "\t" + mz + "\t" + dist + "\t" +
					 relI + "\t" + normI + "\t" + Dscore + "\t" + Iscore + "\t" +
					 score + "\n");
		}
		

		score2pep.buildIonLadders();
		score2pep.matchPeaks(PeakList);
        if(LucXorConfiguration.getScoringAlgorithm() == Constants.CID)
        	score2pep.calcScoreCID(modelingMapCID);
		if(LucXorConfiguration.getScoringAlgorithm() == Constants.HCD)
			score2pep.calcScoreHCD(modelingMapHCD);

		
		mPK.clear();
		mPK = score2pep.getMatchedPeaks();
		mPK.sort(Peak.comparator_mz);
		for(Peak pk : mPK) {
			if(!pk.isMatched()) continue;
			
			String mz = df.format(MathFunctions.roundDouble(pk.getMz(), numDecimals));
			String relI = df.format(MathFunctions.roundDouble(pk.getRelIntensity(), numDecimals));
			String normI = df.format(MathFunctions.roundDouble(pk.getNormIntensity(), numDecimals));
			String Iscore = df.format(MathFunctions.roundDouble(pk.getIntensityScore(), numDecimals));
			String Dscore = df.format(MathFunctions.roundDouble(pk.getDistScore(), numDecimals));
			String score = df.format(MathFunctions.roundDouble(pk.getScore(), numDecimals));
			String dist = df.format(MathFunctions.roundDouble(pk.getDist(), numDecimals));
			
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
	 * @param modelingMapCID CID model
	 * @param modelingMapHCD HCD Moddel
	 */
	public String writeMatchedPks(TMap<Integer, ModelDataCID> modelingMapCID, TMap<Integer, ModelDataHCD> modelingMapHCD) {

		DecimalFormat df = new DecimalFormat("#.####");
		int numDecimals = 4;

        StringBuilder ret = new StringBuilder();

        score1pep.buildIonLadders();
        score1pep.matchPeaks(PeakList);

        if(LucXorConfiguration.getScoringAlgorithm() == Constants.CID) score1pep.calcScoreCID(modelingMapCID);
        if(LucXorConfiguration.getScoringAlgorithm() == Constants.HCD) score1pep.calcScoreHCD(modelingMapHCD);

        List<Peak> mPK = score1pep.getMatchedPeaks();
        mPK.sort(Peak.comparator_mz);
        for(Peak pk : mPK) {
            if(!pk.isMatched()) continue;

            String mz = df.format(MathFunctions.roundDouble(pk.getMz(), numDecimals));
            String relI = df.format(MathFunctions.roundDouble(pk.getRelIntensity(), numDecimals));
            String Iscore = df.format(MathFunctions.roundDouble(pk.getIntensityScore(), numDecimals));
            String Dscore = df.format(MathFunctions.roundDouble(pk.getDistScore(), numDecimals));
            String score = df.format(MathFunctions.roundDouble(pk.getDistScore(), numDecimals));

            ret.append(specId).append("\t1\t").append(score1pep.getModPeptide())
					.append("\t").append(pk.getMatchedIonStr()).append("\t")
					.append(mz).append("\t").append(relI).append("\t").append(Dscore)
					.append("\t").append(Iscore).append("\t").append(score).append("\n");
        }


        score2pep.buildIonLadders();
        score2pep.matchPeaks(PeakList);
        if(LucXorConfiguration.getScoringAlgorithm() == Constants.CID)
        	score2pep.calcScoreCID(modelingMapCID);
        if(LucXorConfiguration.getScoringAlgorithm() == Constants.HCD)
        	score2pep.calcScoreHCD(modelingMapHCD);


        mPK.clear();
        mPK = score2pep.getMatchedPeaks();
        mPK.sort(Peak.comparator_mz);
        for(Peak pk : mPK) {
            if(!pk.isMatched()) continue;

            String mz = df.format(MathFunctions.roundDouble(pk.getMz(), numDecimals));
            String relI = df.format(MathFunctions.roundDouble(pk.getRelIntensity(), numDecimals));
            String Iscore = df.format(MathFunctions.roundDouble(pk.getIntensityScore(), numDecimals));
            String Dscore = df.format(MathFunctions.roundDouble(pk.getDistScore(), numDecimals));
            String score = df.format(MathFunctions.roundDouble(pk.getScore(), numDecimals));

            ret.append(specId).append("\t2\t").append(score2pep.getModPeptide()).append("\t").append(pk.getMatchedIonStr()).append("\t").append(mz).append("\t").append(relI).append("\t").append(Dscore).append("\t").append(Iscore).append("\t").append(score).append("\n");
        }

        return(ret.toString());
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

	public TIntDoubleHashMap getModCoordMap() {
		if(modCoordMap == null){
			modCoordMap = new TIntDoubleHashMap();
		}
		return modCoordMap;
	}

	public List<Peak> getPosPeaks() {
		return posPeaks;
	}

	public List<Peak> getNegPeaks() {
		return negPeaks;
	}

	public Peptide getOrigPep() {
		return origPep;
	}

	public String getPeptideSequence(){
		return origPep.getPeptide();
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

	public void setPeptideSequence(String sequence){
		if(this.origPep == null){
			origPep = new Peptide();
		}

		this.origPep.setPeptide(sequence);
	}

	@Override
	public String toString() {
		return "PSM{" +
				"specId='" + specId + '\'' +
				", srcFile='" + srcFile + '\'' +
				", scanNum=" + scanNum +
				", charge=" + charge +
				", origPep=" + origPep +
				'}';
	}
}
