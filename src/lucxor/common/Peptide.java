/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor.common;

import java.util.*;
import java.util.Map.Entry;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntDoubleMap;
import gnu.trove.map.TMap;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import lucxor.*;
import lucxor.algorithm.ModelDataCID;
import lucxor.algorithm.ModelDataHCD;
import lucxor.utils.Constants;
import lucxor.utils.MathFunctions;
import lucxor.utils.Utils;
import org.apache.commons.math3.util.FastMath;

/**
 *
 * @author dfermin
 */
public class Peptide {

	private String peptide;
	private String modPeptide;
	private int charge = 0;
	private TIntDoubleMap modPosMap = new TIntDoubleHashMap(); // total mass of each modified amino acid
	private final TMap<Integer, String> nonTargetMods = new THashMap<>();// holds the mod coordinates for non-target mods
	private TMap<String, Double> bIons = null, yIons = null;

	private int pepLen = 0;  // peptide's length
	private int numRPS; // number of _R_eported _P_hospho_S_ites
	private int numPPS; // number of _P_otential _P_hospho_S_ites

	private double numPermutations = 0;
	private double numDecoyPermutations = 0;
	private double score = 0;

	private List<Peak> matchedPeaks;

	// this initialization function is called by the PSM class when reading
	// a new entry from a pepXML file.
	void initialize(TIntDoubleHashMap modCoordMap) {

		pepLen = peptide.length();
		peptide = peptide.toUpperCase();

		// remove any modifications from modCoordMap that are fixed modifications

		for(int pos : modCoordMap.keys()) {
			double mass = modCoordMap.get(pos);

			if( (pos == Constants.NTERM_MOD) || (pos == Constants.CTERM_MOD) ) modPosMap.put(pos, mass);
			else {
				String C = Character.toString(peptide.charAt(pos));
				String c = C.toLowerCase();

				if( LucXorConfiguration.getFixedModMap().containsKey(C)) {
				}
				else if( LucXorConfiguration.getTargetModMap().containsKey(C) ) modPosMap.put(pos, mass);
				else if( LucXorConfiguration.getVarModMap().containsKey(c)) nonTargetMods.put(pos, c); // record position of mod.
			}
		}

		// construct the single-character representation of all the modifications
		// in this peptide sequence

		modPeptide = "";

		for(int i = 0; i < pepLen; i++) {
			String AA = Character.toString( peptide.charAt(i) );

			if(nonTargetMods.containsKey(i)) {
				AA = nonTargetMods.get(i);
			}
			else if(modPosMap.containsKey(i)) { // this residue is a target site
				String c = AA;
				if( !Utils.isDecoyResidue(c) ) AA = c.toLowerCase();
			}

			modPeptide += AA;
		}


		// Determine how many potential target-mod sites this peptide has
		numPPS = 0;
		for(int i = 0; i < pepLen; i++) {
			String AA = Character.toString(peptide.charAt(i));
			if(LucXorConfiguration.getTargetModMap().containsKey(AA)) numPPS++;
		}

		// Determine how many reported target-mod sites this peptide has
		numRPS = 0;
		for(int i = 0; i < pepLen; i++) {
			String aa = Character.toString(modPeptide.charAt(i));
			int score = 0;
			if(LucXorConfiguration.getTargetModMap().containsKey(aa.toUpperCase())) score++;
			if(Character.isLowerCase(aa.charAt(0))) score++;

			if(score == 2) numRPS++;
		}


		numPermutations = MathFunctions.combinatorial(numPPS, numRPS);
		calcNumDecoyPermutations();

		boolean is_unambiguous = false;
		if(numPPS == numRPS) is_unambiguous = true;
	}

	// This initialize function is really for declaring a new permutation of a peptide sequence
	// In this function call, the origPepSeq string is the original *MODIFIED* string
	void initialize(String origPep, String curSeq, int z, TIntDoubleHashMap mcp) {

		peptide = origPep;
		pepLen = peptide.length();
		charge = z;
		modPosMap = mcp;

		// deal with the N-terminus and C-terminus values
		if(curSeq.startsWith("[")) {
			modPosMap.put(Constants.NTERM_MOD, LucXorConfiguration.getNtermMass());
			modPeptide = curSeq;
		}
		else modPeptide = curSeq;


		if(curSeq.endsWith("]")) {
			modPosMap.put(Constants.CTERM_MOD, LucXorConfiguration.getCtermMass());
		}
	}


	// Function returns the modified peptide assigned to this object in "TPP format"
	public String getModPepTPP() {
		String ret = "";
		int pepLen_local = pepLen;

		// Append N-term modification (if any is present)
		if(modPeptide.startsWith("[")) {
			int d = (int) Math.round(LucXorConfiguration.getNtermMass());
			ret = "n[" + d + "]";
			pepLen_local += 1; // need this to capture last residue in this function
		}

		// Now label peptide residues
		for(int i = 0; i < pepLen_local; i++) {
			String aa = Character.toString( modPeptide.charAt(i) );

			if(aa.equals("[") || aa.equals("]")) continue;

			// this is a modified amino acid
			if(Character.isLowerCase(modPeptide.charAt(i))) {
				ret += Utils.getTPPresidue(aa, LucXorConfiguration.getNtermMass(), LucXorConfiguration.getCtermMass());
			}
			else ret += aa;
		}

		// Append C-term modification (if any is present)
		if(modPeptide.endsWith("]")) {
			int d = (int) Math.round(LucXorConfiguration.getCtermMass());
			ret += "c[" + d + "]";
		}

		return ret;
	}



	void buildIonLadders() {
		bIons = new THashMap<>();
		yIons = new THashMap<>();


		String b, y;
		double bm, ym; // ion mass
		double bmz, ymz; // ion mass divided by charge z
		double ntermM = 0d, ctermM = 0d;

		final int minLen = 2; // minimum number of residues a fragment must contain

		if(modPosMap.containsKey(Constants.NTERM_MOD)) ntermM = LucXorConfiguration.getNtermMass();
		if(modPosMap.containsKey(Constants.CTERM_MOD)) ctermM = LucXorConfiguration.getCtermMass();

		for(double Z = 1.0; Z < (double) charge; Z++) {
			for(int i = 1; i < pepLen; i++) { // starting at 1 ensures our shortest fragment ion is 2 char long

				String prefix = modPeptide.substring(0,i);
				String suffix = modPeptide.substring(i);

				String prefix_len = Integer.toString(prefix.length());
				String suffix_len = Integer.toString(suffix.length());


				if(prefix.length() >= minLen) {
					b = "b" + prefix_len + ":" + prefix;

					bm = Utils.getFragmentIonMass(b, Z, 0.0) + ntermM;
					bmz = bm / Z;

					if(Z > 1.0) b += "/+" + (int) Z;

					if(bmz > LucXorConfiguration.getMinMZ()) {
						bIons.put(b, bmz);

						// Determine if this ion sequence can under go a neutral loss
						// if it can, record the neutral loss variant
						recordNLIons(b, Z, bm);
					}

				}

				if( (suffix.length() >= minLen) && !suffix.equalsIgnoreCase(modPeptide)) {
					y = "y" + suffix_len + ":" + suffix;
					ym = Utils.getFragmentIonMass(y, Z, Constants.WATER) + ctermM;
					ymz = ym / Z;

					if(Z > 1.0) y += "/+" + (int) Z;

					if(ymz > LucXorConfiguration.getMinMZ()) {
						yIons.put(y, ymz);

						// Determine if this ion sequence can under go a neutral loss
						// if it can, record the neutral loss variant
						recordNLIons(y, Z, ym);
					}
				}
			}
		} // end iteration over charge state
	}



	// Function determines if the passed ion sequence can under go a neutral loss
	// If it can, the NL is recorded
	private void recordNLIons(String ion, double z, double orig_ion_mass) {

		int start = ion.indexOf(':') + 1;
		int end = ion.length();
		if(ion.contains("/")) end = ion.indexOf('/') - 1;

		String x = ion.substring(start, end);

		if( !Utils.isDecoySeq(x) ) {
			for(Entry<String, Double> stringDoubleEntry : LucXorConfiguration.getNeutralLossMap().entrySet()) {

				String candResidues = stringDoubleEntry.getKey().substring( 0, stringDoubleEntry.getKey().indexOf('-') );
				int numCandRes = 0;

				for(int i = 0; i < x.length(); i++){
					String residue = Character.toString( x.charAt(i) );
					if(candResidues.contains(residue)) numCandRes++;
				}

				if(numCandRes > 0) { // this ion contains residues that can result in a NL
					double nl_mass = stringDoubleEntry.getValue();
					String NL_tag = stringDoubleEntry.getKey().substring( stringDoubleEntry.getKey().indexOf('-') );

					String nl_str = ion + NL_tag;
					double mass = orig_ion_mass + nl_mass;
					double mz = MathFunctions.roundDouble((mass / z), 4);

					if(mz > LucXorConfiguration.getMinMZ()) {
						if(ion.startsWith("b")) bIons.put(nl_str, mz);
						if(ion.startsWith("y")) yIons.put(nl_str, mz);
					}
				}
			} // end loop LucXorConfiguration.nlMap
		}
		else {
			// The 'ion' string variable must contain at least 1 decoy residue modification
			// if you got to this point of the code.
			for(Entry<String, Double> e : LucXorConfiguration.getDecoyNeutralLossMap().entrySet()) {
				String NL_tag = e.getKey();
				double nl_mass = e.getValue();

				String nl_str = ion + NL_tag;
				double mass = orig_ion_mass + nl_mass;
				double mz = MathFunctions.roundDouble((mass / z), 4);

				if(mz > LucXorConfiguration.getMinMZ()) {
					if(ion.startsWith("b")) bIons.put(nl_str, mz);
					if(ion.startsWith("y")) yIons.put(nl_str, mz);
				}
			}
		}
	}


	private void calcNumDecoyPermutations() {
		double ctr = 0;
		double p = (double) numPPS;
		double k = (double) numRPS;

		for(int i = 0; i < pepLen; i++) {
			String aa = Character.toString( peptide.charAt(i) );
			if( !LucXorConfiguration.getTargetModMap().containsKey(aa) ) ctr++;
		}

		if(ctr >= p) { // you have enough non-target reidues to make a decoy peptide
			numDecoyPermutations = MathFunctions.combinatorial(ctr, k);
		}
		else numDecoyPermutations = 0;
	}


	public double getNumPerm() {
		return (numPermutations + numDecoyPermutations);
	}


	// function will return a map where the key is each possible permutation
	// of the peptide and the value is the score this permutation recieves later on
	TMap<String, Double> getPermutations(int permType) {
		THashMap<String, Double> ret = new THashMap<>();
		TIntList candModSites = new TIntArrayList();
		List<TIntList> x;
		int i;

		if(permType == 0) { // generate forward (positive) permutations

			// Get candidate "true" sites to modify
			for(i = 0; i < pepLen; i++) {
				String aa = Character.toString( peptide.charAt(i) );
				if( LucXorConfiguration.getTargetModMap().containsKey(aa) ) candModSites.add(i);
			}

			// For the given candidate sites that can undergo modifications,
			// generate all possible permutations involving them.
			x = MathFunctions.getAllCombinations(candModSites, numRPS);

			for(TIntList a : x) {
				String modPep = "";
				if(modPosMap.containsKey(Constants.NTERM_MOD)) modPep = "[";

				for(i = 0; i < pepLen; i++) {
					String aa = Character.toString( peptide.charAt(i) ).toLowerCase();

					if(a.contains(i)) { // site to be modified
						modPep += aa; //Character.toString( Character.toLowerCase(peptide.charAt(i)) );
					}
					else if( nonTargetMods.containsKey(i) ) {
						modPep += aa;
					}
					else modPep += aa.toUpperCase();
				}

				if(modPosMap.containsKey(Constants.CTERM_MOD)) modPep = "]";

				ret.put(modPep, 0d);
			}
		}
		else { // generate decoy (negative) permutations

			for(i = 0; i < pepLen; i++) {
				String AA = Character.toString( peptide.charAt(i) );
				String aa = AA.toLowerCase();
				int score = 0;
				if( !LucXorConfiguration.getTargetModMap().containsKey(AA))
					score++;

				if( !this.nonTargetMods.containsKey(aa) )
					score++;

				if(score == 2) candModSites.add(i);
			}

			// For the given candidate sites that can undergo modifications,
			// generate all possible permutations involving them.
			x = MathFunctions.getAllCombinations(candModSites, numRPS);

			for(TIntList a : x) {
				String modPep = "";
				if(modPosMap.containsKey(Constants.NTERM_MOD)) modPep = "[";

				for(i = 0; i < pepLen; i++) {
					String aa = Character.toString( peptide.charAt(i) ).toLowerCase();

					if(a.contains(i)) { // site to be modified
						String decoyChar = Utils.getDecoySymbol( peptide.charAt(i) );
						modPep += decoyChar;
					}
					else if( nonTargetMods.containsKey(i) ) {
						modPep += aa;
					}
					else modPep += aa.toUpperCase();
				}

				if(modPosMap.containsKey(Constants.CTERM_MOD)) modPep = "]";
				ret.put(modPep, 0d);
			}
		}

		return ret;
	}


	// Function returns the ion ladder for this peptide sequence
	public THashMap<String, Double> getIonLadder() {
		THashMap<String, Double> ret = new THashMap<>();

		ret.putAll(bIons);
		ret.putAll(yIons);

		return ret;
	}


	// Function scores the current peptide
	void calcScoreCID(TMap<Integer, ModelDataCID> modelingMapCID) {

		if(matchedPeaks.isEmpty()) {
			score = 0;
		}
		else {
			ModelDataCID MD = modelingMapCID.get(this.charge);

			// Now compute the scores for these peaks
			double intensityM = 0;
			double intensityU;
			double distM = 0;
			double distU;
			double Iscore;
			double Dscore;
			score = 0;

			for(Peak pk : matchedPeaks) {

				intensityU = MathFunctions
						.logGaussianProb(MD.getMu_int_U(), MD.getVar_int_U(), pk.getNormIntensity());
				distU = MathFunctions
						.logGaussianProb(MD.getMu_dist_U(), MD.getVar_dist_U(), pk.getDist());

				if(pk.getMatchedIonStr().startsWith("b")) {
					intensityM = MathFunctions
							.logGaussianProb(MD.getMu_int_B(), MD.getVar_int_B(), pk.getNormIntensity());
					distM = MathFunctions
							.logGaussianProb(MD.getMu_dist_B(), MD.getVar_dist_B(), pk.getDist());
				}

				if(pk.getMatchedIonStr().startsWith("y")) {
					intensityM = MathFunctions
							.logGaussianProb(MD.getMu_int_Y(), MD.getVar_int_Y(), pk.getNormIntensity());
					distM = MathFunctions
							.logGaussianProb(MD.getMu_dist_Y(), MD.getVar_dist_Y(), pk.getDist());
				}

				Iscore = intensityM - intensityU;
				Dscore = distM - distU;

				double intense_wt = 1.0 / ( 1.0 + FastMath.exp(-Iscore) );
				double x;

				if(Double.isNaN(Dscore) || Double.isInfinite(Dscore)) x = 0;
				else x = intense_wt * Dscore;

				if(x < 0) x = 0; // this prevents the score from going negative

				pk.setScore(x);
				pk.setDistScore(Dscore);
				pk.setIntensityScore(Iscore);

				score += x;
			}
		}
	}



	// Function returns true if the modPeptide assigned to this object is a decoy
	boolean isDecoyPep() {
		boolean ret = false;

		int score = 0;
		for(int i = 0; i < pepLen; i++) {
			String aa = Character.toString( modPeptide.charAt(i) );
			if(Utils.isDecoyResidue(aa))
				score++;
		}

		if(score > 0) ret = true;

		return ret;
	}


	/****************
	 * Function computes the score for the HCD algorithm
	 */
	void calcScoreHCD(TMap<Integer, ModelDataHCD> modelingMapHCD) {

		if(matchedPeaks.isEmpty()) {
			score = 0;
		}
		else {
			ModelDataHCD MD = modelingMapHCD.get(this.charge);

			// we double the error window size for decoys. Otherwise they may not get matched peaks
			isDecoyPep();

			// Now compute the scores for these peaks
			double intensityM;
			double intensityU;
			double distM;
			double distU;
			double Iscore;
			double Dscore;
			score = 0;

			for(Peak pk : matchedPeaks) {
				intensityU = MD.getLogNPdensityInt('n', pk.getNormIntensity());
				distU = 0; // log of uniform distribution between -1 and 1 is zero

				char ionType = pk.getMatchedIonStr().charAt(0);
				intensityM = MD.getLogNPdensityInt(ionType, pk.getNormIntensity());
				distM = MD.getLogNPdensityDistPos(pk.getDist());

				Iscore = intensityM - intensityU;
				Dscore = distM - distU;

				if( Double.isNaN(Iscore) || Double.isInfinite(Iscore)) Iscore = 0;
				if( Double.isNaN(Dscore) || Double.isInfinite(Dscore)) Dscore = 0;

				double x = Iscore + Dscore;
				if(x < 0) x = 0; // this prevents the score from going negative

				pk.setScore(x);
				pk.setDistScore(Dscore);
				pk.setIntensityScore(Iscore);

				score += x;
			}
		}
	}


	/**********
	 * Function identifies all of the peaks in the passed spectrum that can be matched to the
	 * theoretical ions for this peptide
	 * @param obsPeakList
	 */
	public void matchPeaks(Spectrum obsPeakList) {
//        if( !LucXorConfiguration.modelingMapCID.containsKey(this.charge) ) {
//            log.info("\nError! " + peptide + "/+" + this.charge +
//                    ": a CID Model does not exist for this charge state!\nExiting now.\n");
//            System.exit(0);
//        }

		double matchErr;
		double a, b;

		matchedPeaks = new ArrayList<>();


		int N = bIons.size() + yIons.size() ;

		// The same observed peak might be matchable to multiple theoretical ions for this peptide.
		// This HashMap is used to keep track of the best possible match for an observed peak
		// k = observed peak's m/z value
		// v = the best match for this observed peak based upon the abs(mz_dist)
		TDoubleObjectHashMap<Peak> bestMatchMap = new TDoubleObjectHashMap<>(N);


		// Try to match y-ions first, they tend to be of higher intensity which is what we want
		for(Entry<String, Double> stringDoubleEntry : yIons.entrySet()) {
			double theo_mz  = stringDoubleEntry.getValue();

			if(LucXorConfiguration.getMs2tolUnits() == Constants.PPM_UNITS) {
				double ppmErr = LucXorConfiguration.getMs2tol() / Constants.PPM;
				matchErr = theo_mz * ppmErr;
			}
			else matchErr = LucXorConfiguration.getMs2tol();

			matchErr *= 0.5; // split in half
			a = MathFunctions.roundDouble((theo_mz - matchErr), 4);
			b = MathFunctions.roundDouble((theo_mz + matchErr), 4);

			// Hold candidate matching peaks here
			ArrayList<Peak> cand = new ArrayList<>();
			for(int i = 0; i < obsPeakList.N; i++) {
				double obsMZ = MathFunctions.roundDouble( obsPeakList.mz[i], 4 );

				if( (obsMZ >= a) && (obsMZ <= b) ) {
					cand.add( obsPeakList.getPeakClassInstance(i) );
				}
			}

			if(!cand.isEmpty()) {
				// Take the most intense peak as the best match
				cand.sort(Peak.comparator_intensity_hi2low);
				Peak pk = cand.get(0);

				pk.setMatched(true);
				pk.setMatchedIonStr(stringDoubleEntry.getKey());
				pk.setDist(pk.getMz() - theo_mz);

				// Check to see if you have already assigned this peak to a theoretical ion
				if(bestMatchMap.containsKey(pk.getMz())) {
					Peak oldPK = bestMatchMap.get(pk.getMz());
					if(Math.abs(oldPK.getDist()) > (Math.abs(pk.getDist()))) bestMatchMap.put(pk.getMz(), pk);
				}
				else bestMatchMap.put(pk.getMz(), pk);
			}
			cand.clear();
		}


		// Now try to match b-ions
		for(Entry<String, Double> stringDoubleEntry : bIons.entrySet()) {
			double theo_mz  = stringDoubleEntry.getValue();

			if(LucXorConfiguration.getMs2tolUnits() == Constants.PPM_UNITS) {
				double ppmErr = LucXorConfiguration.getMs2tol() / Constants.PPM;
				matchErr = theo_mz * ppmErr;
			}
			else matchErr = LucXorConfiguration.getMs2tol();

			matchErr *= 0.5; // split in half
			a = MathFunctions.roundDouble((theo_mz - matchErr), 4);
			b = MathFunctions.roundDouble((theo_mz + matchErr), 4);

			// Hold candidate matching peaks here
			ArrayList<Peak> cand = new ArrayList<>();
			for(int i = 0; i < obsPeakList.N; i++) {
				double obsMZ = MathFunctions.roundDouble( obsPeakList.mz[i], 4 );

				if( (obsMZ >= a) && (obsMZ <= b) ) {
					cand.add( obsPeakList.getPeakClassInstance(i) );
				}
			}

			if(!cand.isEmpty()) {
				// Take the most intense peak as the best match
				cand.sort(Peak.comparator_intensity_hi2low);
				Peak pk = cand.get(0);

				pk.setMatched(true);
				pk.setMatchedIonStr(stringDoubleEntry.getKey());
				pk.setDist(pk.getMz() - theo_mz);

				// Check to see if you have already assigned this peak to a theoretical ion
				if(bestMatchMap.containsKey(pk.getMz())) {
					Peak oldPK = bestMatchMap.get(pk.getMz());
					if(Math.abs(oldPK.getDist()) > (Math.abs(pk.getDist()))) bestMatchMap.put(pk.getMz(), pk);
				}
				else bestMatchMap.put(pk.getMz(), pk);
			}
			cand.clear();
		}

		// MM holds the best match for each observed peak.
		for(double mz : bestMatchMap.keys()) matchedPeaks.add( bestMatchMap.get(mz) );
		bestMatchMap.clear();
	}

	public String getPeptide() {
		return peptide;
	}

	public void setPeptide(String peptide) {
		this.peptide = peptide;
	}

	public String getModPeptide() {
		return modPeptide;
	}

	public void setCharge(int charge) {
		this.charge = charge;
	}

	public TMap<Integer, String> getNonTargetMods() {
		return nonTargetMods;
	}

	public int getPepLen() {
		return pepLen;
	}

	public int getNumRPS() {
		return numRPS;
	}

	public void setNumRPS(int numRPS) {
		this.numRPS = numRPS;
	}

	public int getNumPPS() {
		return numPPS;
	}

	public void setNumPPS(int numPPS) {
		this.numPPS = numPPS;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}

	public List<Peak> getMatchedPeaks() {
		return matchedPeaks;
	}

}
