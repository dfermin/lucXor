/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import java.util.*;
import java.util.Map.Entry;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import org.apache.commons.math3.util.FastMath;

/**
 *
 * @author dfermin
 */
class Peptide {

	private String peptide;
	private String modPeptide;
	private int charge = 0;
	private TIntDoubleHashMap modPosMap = new TIntDoubleHashMap(); // total mass of each modified amino acid
	private THashMap<Integer, String> nonTargetMods = new THashMap();// holds the mod coordinates for non-target mods
	private THashMap<String, Double> bIons = null, yIons = null;

	private int pepLen = 0;  // peptide's length
	private int numRPS; // number of _R_eported _P_hospho_S_ites
	private int numPPS; // number of _P_otential _P_hospho_S_ites

	private boolean is_unambiguous = false;

	private double numPermutations = 0;
	private double numDecoyPermutations = 0;
	private double score = 0;

	private ArrayList<Peak> matchedPeaks;

	// this initialization function is called by the PSM class when reading
	// a new entry from a pepXML file.
	void initialize(TIntDoubleHashMap modCoordMap) {

		pepLen = peptide.length();

		String origPepStr = peptide.toUpperCase();
		peptide = origPepStr;

		// remove any modifications from modCoordMap that are fixed modifications

		for(int pos : modCoordMap.keys()) {
			double mass = modCoordMap.get(pos);

			if( (pos == Constants.NTERM_MOD) || (pos == Constants.CTERM_MOD) ) modPosMap.put(pos, mass);
			else {
				String C = Character.toString(peptide.charAt(pos));
				String c = C.toLowerCase();

				if( Globals.fixedModMap.containsKey(C)) continue;
				else if( Globals.TargetModMap.containsKey(C) ) modPosMap.put(pos, mass);
				else if( Globals.varModMap.containsKey(c)) nonTargetMods.put(pos, c); // record position of mod.
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
				if( !Globals.isDecoyResidue(c) ) AA = c.toLowerCase();
			}

			modPeptide += AA;
		}


		// Determine how many potential target-mod sites this peptide has
		numPPS = 0;
		for(int i = 0; i < pepLen; i++) {
			String AA = Character.toString(peptide.charAt(i));
			if(Globals.TargetModMap.containsKey(AA)) numPPS++;
		}

		// Determine how many reported target-mod sites this peptide has
		numRPS = 0;
		for(int i = 0; i < pepLen; i++) {
			String aa = Character.toString(modPeptide.charAt(i));
			int score = 0;
			if(Globals.TargetModMap.containsKey(aa.toUpperCase())) score++;
			if(Character.isLowerCase(aa.charAt(0))) score++;

			if(score == 2) numRPS++;
		}


		numPermutations = MathFunctions.combinatorial(numPPS, numRPS);
		calcNumDecoyPermutations();

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
			modPosMap.put(Constants.NTERM_MOD, Globals.ntermMass);
			modPeptide = curSeq;
		}
		else modPeptide = curSeq;


		if(curSeq.endsWith("]")) {
			modPosMap.put(Constants.CTERM_MOD, Globals.ctermMass);
		}
	}


	// Function returns the modified peptide assigned to this object in "TPP format"
	public String getModPepTPP() {
		String ret = "";
		int pepLen_local = pepLen;

		// Append N-term modification (if any is present)
		if(modPeptide.startsWith("[")) {
			int d = (int) Math.round(Globals.ntermMass);
			ret = "n[" + String.valueOf( d ) + "]";
			pepLen_local += 1; // need this to capture last residue in this function
		}

		// Now label peptide residues
		for(int i = 0; i < pepLen_local; i++) {
			String aa = Character.toString( modPeptide.charAt(i) );

			if(aa.equals("[") || aa.equals("]")) continue;

			// this is a modified amino acid
			if(Character.isLowerCase(modPeptide.charAt(i))) {
				ret += Globals.getTPPresidue(aa);
			}
			else ret += aa;
		}

		// Append C-term modification (if any is present)
		if(modPeptide.endsWith("]")) {
			int d = (int) Math.round(Globals.ctermMass);
			ret += "c[" + String.valueOf(d) + "]";
		}

		return ret;
	}



	void buildIonLadders() {
		bIons = new THashMap();
		yIons = new THashMap();


		String b = null, y = null;
		double bm = 0d, ym = 0d; // ion mass
		double bmz = 0d, ymz = 0d; // ion mass divided by charge z
		double ntermM = 0d, ctermM = 0d;

		final int minLen = 2; // minimum number of residues a fragment must contain

		if(modPosMap.containsKey(Constants.NTERM_MOD)) ntermM = Globals.ntermMass;
		if(modPosMap.containsKey(Constants.CTERM_MOD)) ctermM = Globals.ctermMass;

		for(double Z = 1.0; Z < (double) charge; Z++) {
			for(int i = 1; i < pepLen; i++) { // starting at 1 ensures our shortest fragment ion is 2 char long

				String prefix = modPeptide.substring(0,i);
				String suffix = modPeptide.substring(i);

				String prefix_len = Integer.toString(prefix.length());
				String suffix_len = Integer.toString(suffix.length());


				if(prefix.length() >= minLen) {
					b = "b" + prefix_len + ":" + prefix;

					bm = Globals.getFragmentIonMass(b, Z, 0.0) + ntermM;
					bmz = bm / Z;

					if(Z > 1.0) b += "/+" + Integer.toString((int)Z);

					if(bmz > Globals.minMZ) {
						bIons.put(b, bmz);

						// Determine if this ion sequence can under go a neutral loss
						// if it can, record the neutral loss variant
						recordNLIons(b, Z, bm);
					}

				}

				if( (suffix.length() >= minLen) && !suffix.equalsIgnoreCase(modPeptide)) {
					y = "y" + suffix_len + ":" + suffix;
					ym = Globals.getFragmentIonMass(y, Z, Constants.WATER) + ctermM;
					ymz = ym / Z;

					if(Z > 1.0) y += "/+" + Integer.toString((int)Z);

					if(ymz > Globals.minMZ) {
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

		int start = ion.indexOf(":") + 1;
		int end = ion.length();
		if(ion.contains("/")) end = ion.indexOf("/") - 1;

		String x = ion.substring(start, end);

		if( !Globals.isDecoySeq(x) ) {
			for(String nlr : Globals.nlMap.keySet()) {

				String candResidues = nlr.substring( 0, nlr.indexOf("-") );
				int numCandRes = 0;

				for(int i = 0; i < x.length(); i++){
					String residue = Character.toString( x.charAt(i) );
					if(candResidues.contains(residue)) numCandRes++;
				}

				if(numCandRes > 0) { // this ion contains residues that can result in a NL
					double nl_mass = Globals.nlMap.get(nlr);
					String NL_tag = nlr.substring( nlr.indexOf("-") );

					String nl_str = ion + NL_tag;
					double mass = orig_ion_mass + nl_mass;
					double mz = MathFunctions.roundDouble((mass / z), 4);

					if(mz > Globals.minMZ) {
						if(ion.startsWith("b")) bIons.put(nl_str, mz);
						if(ion.startsWith("y")) yIons.put(nl_str, mz);
					}
				}
			} // end loop Globals.nlMap
		}
		else {
			// The 'ion' string variable must contain at least 1 decoy residue modification
			// if you got to this point of the code.
			for(Entry<String, Double> e : Globals.decoyNLmap.entrySet()) {
				String NL_tag = e.getKey();
				double nl_mass = e.getValue();

				String nl_str = ion + NL_tag;
				double mass = orig_ion_mass + nl_mass;
				double mz = MathFunctions.roundDouble((mass / z), 4);

				if(mz > Globals.minMZ) {
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
			if( !Globals.TargetModMap.containsKey(aa) ) ctr++;
		}

		if(ctr >= p) { // you have enough non-target reidues to make a decoy peptide
			numDecoyPermutations = MathFunctions.combinatorial(ctr, k);
		}
		else numDecoyPermutations = 0;
	}


	double getNumPerm() {
		return (numPermutations + numDecoyPermutations);
	}


	void printIons() {
		for(String p : bIons.keySet()) {
			double m = MathFunctions.roundDouble(bIons.get(p), 4);
			System.out.println(modPeptide + "\t" + p + "\t" + m);
		}

		for(String p : yIons.keySet()) {
			double m = MathFunctions.roundDouble(yIons.get(p), 4);
			System.out.println(modPeptide + "\t" + p + "\t" + m);
		}
	}


	// function will return a map where the key is each possible permutation
	// of the peptide and the value is the score this permutation recieves later on
	THashMap<String, Double> getPermutations(int permType) {
		THashMap<String, Double> ret = new THashMap();
		TIntArrayList candModSites = new TIntArrayList();
		ArrayList<TIntArrayList> x = new ArrayList();
		int i = 0;

		if(permType == 0) { // generate forward (positive) permutations

			// Get candidate "true" sites to modify
			for(i = 0; i < pepLen; i++) {
				String aa = Character.toString( peptide.charAt(i) );
				if( Globals.TargetModMap.containsKey(aa) ) candModSites.add(i);
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
				if( !Globals.TargetModMap.containsKey(AA) ) score++;

				if( !this.nonTargetMods.containsKey(aa) ) score++;

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
						String decoyChar = Globals.getDecoySymbol( peptide.charAt(i) );
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
		THashMap<String, Double> ret = new THashMap();

		ret.putAll(bIons);
		ret.putAll(yIons);

		return ret;
	}


	// Function scores the current peptide
	void calcScoreCID() {

		if(matchedPeaks.isEmpty()) {
			score = 0;
		}
		else {
			ModelDataCID MD = Globals.modelingMap_CID.get(this.charge);

			// Now compute the scores for these peaks
			double intensityM = 0;
			double intensityU = 0;
			double distM = 0;
			double distU = 0;
			double Iscore = 0d;
			double Dscore = 0d;
			score = 0;

			for(Peak pk : matchedPeaks) {

				intensityU = MathFunctions
						.log_gaussianProb(MD.mu_int_U, MD.var_int_U, pk.getNormIntensity());
				distU = MathFunctions
						.log_gaussianProb(MD.mu_dist_U, MD.var_dist_U, pk.getDist());

				if(pk.getMatchedIonStr().startsWith("b")) {
					intensityM = MathFunctions
							.log_gaussianProb(MD.mu_int_B, MD.var_int_B, pk.getNormIntensity());
					distM = MathFunctions
							.log_gaussianProb(MD.mu_dist_B, MD.var_dist_B, pk.getDist());
				}

				if(pk.getMatchedIonStr().startsWith("y")) {
					intensityM = MathFunctions
							.log_gaussianProb(MD.mu_int_Y, MD.var_int_Y, pk.getNormIntensity());
					distM = MathFunctions
							.log_gaussianProb(MD.mu_dist_Y, MD.var_dist_Y, pk.getDist());
				}

				Iscore = intensityM - intensityU;
				Dscore = distM - distU;

				double intense_wt = 1.0 / ( 1.0 + FastMath.exp(-Iscore) );
				double x = 0d;

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
			if(Globals.isDecoyResidue(aa)) score++;
		}

		if(score > 0) ret = true;

		return ret;
	}


	/****************
	 * Function computes the score for the HCD algorithm
	 */
	void calcScoreHCD() {

		if(matchedPeaks.isEmpty()) {
			score = 0;
		}
		else {
			ModelDataHCD MD = Globals.modelingMap_HCD.get(this.charge);

			double matchErr = 0;
			double decoyPadding = 1;
			double a = 0, b = 0;

			// we double the error window size for decoys. Otherwise they may not get matched peaks
			if(isDecoyPep()) decoyPadding = 2.0;


			// Now compute the scores for these peaks
			double intensityM = 0;
			double intensityU = 0;
			double distM = 0;
			double distU = 0;
			double Iscore = 0d;
			double Dscore = 0d;
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
	public void matchPeaks(SpectrumClass obsPeakList) {
//        if( !Globals.modelingMap_CID.containsKey(this.charge) ) {
//            log.info("\nError! " + peptide + "/+" + this.charge +
//                    ": a CID Model does not exist for this charge state!\nExiting now.\n");
//            System.exit(0);
//        }

		double matchErr = 0;
		double a = 0, b = 0;

		matchedPeaks = new ArrayList<>();


		int N = bIons.size() + yIons.size() ;

		// The same observed peak might be matchable to multiple theoretical ions for this peptide.
		// This HashMap is used to keep track of the best possible match for an observed peak
		// k = observed peak's m/z value
		// v = the best match for this observed peak based upon the abs(mz_dist)
		TDoubleObjectHashMap<Peak> bestMatchMap = new TDoubleObjectHashMap<>(N);


		// Try to match y-ions first, they tend to be of higher intensity which is what we want
		for(String theo_ion : yIons.keySet()) {
			double theo_mz  = yIons.get(theo_ion);

			if(Globals.ms2tol_units == Constants.PPM_UNITS) {
				double ppmErr = Globals.ms2tol / Constants.PPM;
				matchErr = theo_mz * ppmErr;
			}
			else matchErr = Globals.ms2tol;

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
				pk.setMatchedIonStr(theo_ion);
				pk.setDist(pk.getMz() - theo_mz);

				// Check to see if you have already assigned this peak to a theoretical ion
				if(bestMatchMap.containsKey(pk.getMz())) {
					Peak oldPK = bestMatchMap.get(pk.getMz());
					if(Math.abs(oldPK.getDist()) > (Math.abs(pk.getDist()))) bestMatchMap.put(pk.getMz(), pk);
					oldPK = null;
				}
				else bestMatchMap.put(pk.getMz(), pk);
			}
			cand.clear();
		}


		// Now try to match b-ions
		for(String theo_ion : bIons.keySet()) {
			double theo_mz  = bIons.get(theo_ion);

			if(Globals.ms2tol_units == Constants.PPM_UNITS) {
				double ppmErr = Globals.ms2tol / Constants.PPM;
				matchErr = theo_mz * ppmErr;
			}
			else matchErr = Globals.ms2tol;

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
				pk.setMatchedIonStr(theo_ion);
				pk.setDist(pk.getMz() - theo_mz);

				// Check to see if you have already assigned this peak to a theoretical ion
				if(bestMatchMap.containsKey(pk.getMz())) {
					Peak oldPK = bestMatchMap.get(pk.getMz());
					if(Math.abs(oldPK.getDist()) > (Math.abs(pk.getDist()))) bestMatchMap.put(pk.getMz(), pk);
					oldPK = null;
				}
				else bestMatchMap.put(pk.getMz(), pk);
			}
			cand.clear();
		}

		// MM holds the best match for each observed peak.
		for(double mz : bestMatchMap.keys()) matchedPeaks.add( bestMatchMap.get(mz) );
		bestMatchMap.clear();
		bestMatchMap = null;
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

	public void setModPeptide(String modPeptide) {
		this.modPeptide = modPeptide;
	}

	public int getCharge() {
		return charge;
	}

	public void setCharge(int charge) {
		this.charge = charge;
	}

	public TIntDoubleHashMap getModPosMap() {
		return modPosMap;
	}

	public void setModPosMap(TIntDoubleHashMap modPosMap) {
		this.modPosMap = modPosMap;
	}

	public THashMap<Integer, String> getNonTargetMods() {
		return nonTargetMods;
	}

	public void setNonTargetMods(THashMap<Integer, String> nonTargetMods) {
		this.nonTargetMods = nonTargetMods;
	}

	public THashMap<String, Double> getbIons() {
		return bIons;
	}

	public void setbIons(THashMap<String, Double> bIons) {
		this.bIons = bIons;
	}

	public THashMap<String, Double> getyIons() {
		return yIons;
	}

	public void setyIons(THashMap<String, Double> yIons) {
		this.yIons = yIons;
	}

	public int getPepLen() {
		return pepLen;
	}

	public void setPepLen(int pepLen) {
		this.pepLen = pepLen;
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

	public boolean isIs_unambiguous() {
		return is_unambiguous;
	}

	public void setIs_unambiguous(boolean is_unambiguous) {
		this.is_unambiguous = is_unambiguous;
	}

	public double getNumPermutations() {
		return numPermutations;
	}

	public void setNumPermutations(double numPermutations) {
		this.numPermutations = numPermutations;
	}

	public double getNumDecoyPermutations() {
		return numDecoyPermutations;
	}

	public void setNumDecoyPermutations(double numDecoyPermutations) {
		this.numDecoyPermutations = numDecoyPermutations;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}

	public ArrayList<Peak> getMatchedPeaks() {
		return matchedPeaks;
	}

	public void setMatchedPeaks(ArrayList<Peak> matchedPeaks) {
		this.matchedPeaks = matchedPeaks;
	}
}
