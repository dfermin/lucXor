/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucXor;

import com.google.common.collect.ArrayListMultimap;

import java.util.*;
import java.util.Map.Entry;

import gnu.trove.impl.hash.TIntDoubleHash;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import gnu.trove.procedure.TIntDoubleProcedure;
import org.apache.commons.math3.util.FastMath;

/**
 *
 * @author dfermin
 */
class Peptide {
	String peptide;
	String modPeptide;
	int charge;
	TIntDoubleHashMap modPosMap = null; // total mass of each modified amino acid
	THashMap<Integer, String> nonTargetMods = null; // holds the mod coordinates for non-target mods
	THashMap<String, Double> b_ions = null, y_ions = null;
	
	int pepLen = 0;  // peptide's length
	int numRPS; // number of _R_eported _P_hospho_S_ites
	int numPPS; // number of _P_otential _P_hospho_S_ites
	
	boolean is_unambiguous = false;
	
	double numPermutations = 0;
	double numDecoyPermutations = 0;
	double score = 0;
	
	ArrayList<PeakClass> matchedPeaks;
	
	Peptide() {
		modPosMap = new TIntDoubleHashMap();
		nonTargetMods = new THashMap();
		charge = 0;
		numPermutations = 0;
		score = 0;
	}

	
	// this initialization function is called by the PSM class when reading
	// a new entry from a pepXML file.
	void initialize(TIntDoubleHashMap modCoordMap) {
		
		pepLen = peptide.length();
		
		String origPepStr = peptide.toUpperCase();
		peptide = origPepStr; 
		
		// remove any modifications from modCoordMap that are fixed modifications

        for(int pos : modCoordMap.keys()) {
			double mass = modCoordMap.get(pos);

			if( (pos == constants.NTERM_MOD) || (pos == constants.CTERM_MOD) ) modPosMap.put(pos, mass);
			else {
				String C = Character.toString(peptide.charAt(pos));
				String c = C.toLowerCase();
				
				if( globals.fixedModMap.containsKey(C)) continue;
				else if( globals.targetModMap.containsKey(C) ) modPosMap.put(pos, mass);
				else if( globals.varModMap.containsKey(c)) nonTargetMods.put(pos, c); // record position of mod.
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
				if( !globals.isDecoyResidue(c) ) AA = c.toLowerCase();
			}
			
			modPeptide += AA;
		}

		
		// Determine how many potential target-mod sites this peptide has
		numPPS = 0;
		for(int i = 0; i < pepLen; i++) {
			String AA = Character.toString(peptide.charAt(i));
			if(globals.targetModMap.containsKey(AA)) numPPS++;
		}
		
		// Determine how many reported target-mod sites this peptide has
		numRPS = 0;
		for(int i = 0; i < pepLen; i++) {
			String aa = Character.toString(modPeptide.charAt(i));
			int score = 0;
			if(globals.targetModMap.containsKey(aa.toUpperCase())) score++;
			if(Character.isLowerCase(aa.charAt(0))) score++;
			
			if(score == 2) numRPS++;
		}
		
		
		numPermutations = globals.SF.combinatorial(numPPS, numRPS);
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
			modPosMap.put(constants.NTERM_MOD, globals.ntermMass);
			modPeptide = curSeq;
		}
		else modPeptide = curSeq;
		
		
		if(curSeq.endsWith("]")) {
			modPosMap.put(constants.CTERM_MOD, globals.ctermMass);
		}
	}

	
	// Function returns the modified peptide assigned to this object in "TPP format"
	public String getModPepTPP() {
		String ret = "";
		
		// Append N-term modification (if any is present)
		if(modPeptide.startsWith("[")) {
			int d = (int) Math.round(globals.ntermMass);
			ret = "n[" + String.valueOf( d ) + "]=";
		}
		
		
		// Now label peptide residues
		for(int i = 0; i < pepLen; i++) {
			String aa = Character.toString( modPeptide.charAt(i) );
			
			if(modPosMap.containsKey(i)) {
				ret += globals.getTPPresidue(aa);
			}
			else ret += aa;
		}
		
		// Append C-term modification (if any is present)
		if(modPeptide.endsWith("]")) {
			int d = (int) Math.round(globals.ctermMass);
			ret += "=c[" + String.valueOf(d) + "]";
		}
		
		return ret;
	}

	
	
	void build_ion_ladders() {
		b_ions = new THashMap();
		y_ions = new THashMap();

		
		String b = null, y = null;
		double bm = 0d, ym = 0d; // ion mass
		double bmz = 0d, ymz = 0d; // ion mass divided by charge z
		double ntermM = 0d, ctermM = 0d;
		
		final int minLen = 2; // minimum number of residues a fragment must contain
		
		if(modPosMap.containsKey(constants.NTERM_MOD)) ntermM = globals.ntermMass;
		if(modPosMap.containsKey(constants.CTERM_MOD)) ctermM = globals.ctermMass;
		
		for(double Z = 1.0; Z < (double) charge; Z++) {
			for(int i = 1; i < pepLen; i++) { // starting at 1 ensures our shortest fragment ion is 2 char long
				
				String prefix = modPeptide.substring(0,i);
				String suffix = modPeptide.substring(i);

				String prefix_len = Integer.toString(prefix.length());
				String suffix_len = Integer.toString(suffix.length());
				
				
				if(prefix.length() >= minLen) {
					b = "b" + prefix_len + ":" + prefix;
					
					bm = globals.getFragmentIonMass(b, Z, 0.0) + ntermM;
					bmz = bm / Z;
					
					if(Z > 1.0) b += "/+" + Integer.toString((int)Z);
					
					if(bmz > globals.minMZ) {
						b_ions.put(b, bmz);

						// Determine if this ion sequence can under go a neutral loss
						// if it can, record the neutral loss variant
						record_NL_ions(b, Z, bm);
					}
					
				}

				if( (suffix.length() >= minLen) && !suffix.equalsIgnoreCase(modPeptide)) {
					y = "y" + suffix_len + ":" + suffix;
					ym = globals.getFragmentIonMass(y, Z, constants.WATER) + ctermM;
					ymz = ym / Z;
					
					if(Z > 1.0) y += "/+" + Integer.toString((int)Z);
					
					if(ymz > globals.minMZ) {
						y_ions.put(y, ymz);

						// Determine if this ion sequence can under go a neutral loss
						// if it can, record the neutral loss variant
                        record_NL_ions(y, Z, ym);
					}
				}
			}
		} // end iteration over charge state
	}
	
	
	
	// Function determines if the passed ion sequence can under go a neutral loss
	// If it can, the NL is recorded
	private void record_NL_ions(String ion, double z, double orig_ion_mass) {
		
		int start = ion.indexOf(":") + 1;
		int end = ion.length();
		if(ion.contains("/")) end = ion.indexOf("/") - 1;
		
		String x = ion.substring(start, end);

        if( !globals.isDecoySeq(x) ) {
            for(String nlr : globals.nlMap.keySet()) {

                String candResidues = nlr.substring( 0, nlr.indexOf("-") );
                int numCandRes = 0;

                for(int i = 0; i < x.length(); i++){
                    String residue = Character.toString( x.charAt(i) );
                    if(candResidues.contains(residue)) numCandRes++;
                }

                if(numCandRes > 0) { // this ion contains residues that can result in a NL
                    double nl_mass = globals.nlMap.get(nlr);
                    String NL_tag = nlr.substring( nlr.indexOf("-") );

                    String nl_str = ion + NL_tag;
                    double mass = orig_ion_mass + nl_mass;
                    double mz = globals.round_dbl((mass / z), 4);

                    if(mz > globals.minMZ) {
                        if(ion.startsWith("b")) b_ions.put(nl_str, mz);
                        if(ion.startsWith("y")) y_ions.put(nl_str, mz);
                    }
                }
            } // end loop globals.nlMap
        }
        else {
            // The 'ion' string variable must contain at least 1 decoy residue modification
            // if you got to this point of the code.
            for(Entry<String, Double> e : globals.decoyNLmap.entrySet()) {
                String NL_tag = e.getKey();
                double nl_mass = e.getValue();

                String nl_str = ion + NL_tag;
                double mass = orig_ion_mass + nl_mass;
                double mz = globals.round_dbl((mass / z), 4);

                if(mz > globals.minMZ) {
                    if(ion.startsWith("b")) b_ions.put(nl_str, mz);
                    if(ion.startsWith("y")) y_ions.put(nl_str, mz);
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
			if( !globals.targetModMap.containsKey(aa) ) ctr++;
		}
		
		if(ctr >= p) { // you have enough non-target reidues to make a decoy peptide
			numDecoyPermutations = globals.SF.combinatorial(ctr, k);
		}
		else numDecoyPermutations = 0;
	}
	
	
	double getNumPerm() {
		return (numPermutations + numDecoyPermutations);
	}
	
	
	void printIons() {
		for(String p : b_ions.keySet()) {
			double m = globals.round_dbl(b_ions.get(p), 4);
			System.out.println(modPeptide + "\t" + p + "\t" + m);
		}
		
		for(String p : y_ions.keySet()) {
			double m = globals.round_dbl(y_ions.get(p), 4);
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
				if( globals.targetModMap.containsKey(aa) ) candModSites.add(i);
			}
			
			// For the given candidate sites that can undergo modifications,
			// generate all possible permutations involving them.
			x = globals.SF.getAllCombinations(candModSites, numRPS);
			
			for(TIntList a : x) {
				String modPep = "";
				if(modPosMap.containsKey(constants.NTERM_MOD)) modPep = "[";
				
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
				
				if(modPosMap.containsKey(constants.CTERM_MOD)) modPep = "]";
				
				ret.put(modPep, 0d);
			}
		}
		else { // generate decoy (negative) permutations
			
			for(i = 0; i < pepLen; i++) {
				String AA = Character.toString( peptide.charAt(i) );
				String aa = AA.toLowerCase();
                int score = 0;
                if( !globals.targetModMap.containsKey(AA) ) score++;

                if( !this.nonTargetMods.containsKey(aa) ) score++;

                if(score == 2) candModSites.add(i);
			}
			
			// For the given candidate sites that can undergo modifications,
			// generate all possible permutations involving them.
			x = globals.SF.getAllCombinations(candModSites, numRPS);
			
			for(TIntList a : x) {
				String modPep = "";
				if(modPosMap.containsKey(constants.NTERM_MOD)) modPep = "[";
				
				for(i = 0; i < pepLen; i++) {
					String aa = Character.toString( peptide.charAt(i) ).toLowerCase();
					
					if(a.contains(i)) { // site to be modified
						String decoyChar = globals.getDecoySymbol( peptide.charAt(i) );
						modPep += decoyChar;
					}
					else if( nonTargetMods.containsKey(i) ) {
						modPep += aa;
					}
					else modPep += aa.toUpperCase();
				}
				
				if(modPosMap.containsKey(constants.CTERM_MOD)) modPep = "]";
				ret.put(modPep, 0d);
			}
		}
		
		return ret;
	}

	
	// Function returns the ion ladder for this peptide sequence
	public THashMap<String, Double> getIonLadder() {
		THashMap<String, Double> ret = new THashMap();
		
		ret.putAll(b_ions);
		ret.putAll(y_ions);
		
		return ret;
	}

	
	// Function scores the current peptide
	void calcScore_CID() {

		if(matchedPeaks.isEmpty()) {
			score = 0;
		}
		else {
            ModelData_CID MD = globals.modelingMap_CID.get(this.charge);

            // Now compute the scores for these peaks
            double intensityM = 0;
            double intensityU = 0;
            double distM = 0;
            double distU = 0;
			double Iscore = 0d;
            double Dscore = 0d;
			score = 0;

            for(PeakClass pk : matchedPeaks) {

                intensityU = globals.SF.log_gaussianProb(MD.mu_int_U, MD.var_int_U, pk.norm_intensity);
                distU = globals.SF.log_gaussianProb(MD.mu_dist_U, MD.var_dist_U, pk.dist);

                if(pk.matchedIonStr.startsWith("b")) {
                    intensityM = globals.SF.log_gaussianProb(MD.mu_int_B, MD.var_int_B, pk.norm_intensity);
                    distM = globals.SF.log_gaussianProb(MD.mu_dist_B, MD.var_dist_B, pk.dist);
                }

                if(pk.matchedIonStr.startsWith("y")) {
                    intensityM = globals.SF.log_gaussianProb(MD.mu_int_Y, MD.var_int_Y, pk.norm_intensity);
                    distM = globals.SF.log_gaussianProb(MD.mu_dist_Y, MD.var_dist_Y, pk.dist);
                }

                Iscore = intensityM - intensityU;
                Dscore = distM - distU;

				double intense_wt = 1.0 / ( 1.0 + FastMath.exp(-Iscore) );
				double x = 0d;

				if(Double.isNaN(Dscore) || Double.isInfinite(Dscore)) x = 0;
				else x = intense_wt * Dscore;

                if(x < 0) x = 0; // this prevents the score from going negative

				pk.score = x;
                pk.distScore = Dscore;
                pk.intensityScore = Iscore;

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
			if(globals.isDecoyResidue(aa)) score++;
		}
		
		if(score > 0) ret = true;
		
		return ret;
	}


    /****************
     * Function computes the score for the HCD algorithm
     */
	void calcScore_HCD() {

        if(matchedPeaks.isEmpty()) {
            score = 0;
        }
        else {
            ModelData_HCD MD = globals.modelingMap_HCD.get(this.charge);

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

            for(PeakClass pk : matchedPeaks) {
                intensityU = MD.getLogNPdensityInt('n', pk.norm_intensity);
                distU = 0; // log of uniform distribution between -1 and 1 is zero

                char ionType = pk.matchedIonStr.charAt(0);
                intensityM = MD.getLogNPdensityInt(ionType, pk.norm_intensity);
                distM = MD.getLogNPdensityDistPos(pk.dist);

                Iscore = intensityM - intensityU;
                Dscore = distM - distU;

                if( Double.isNaN(Iscore) || Double.isInfinite(Iscore)) Iscore = 0;
                if( Double.isNaN(Dscore) || Double.isInfinite(Dscore)) Dscore = 0;

                double x = Iscore + Dscore;
                if(x < 0) x = 0; // this prevents the score from going negative

                pk.score = x;
                pk.distScore = Dscore;
                pk.intensityScore = Iscore;

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
//        if( !globals.modelingMap_CID.containsKey(this.charge) ) {
//            System.err.println("\nError! " + peptide + "/+" + this.charge +
//                    ": a CID Model does not exist for this charge state!\nExiting now.\n");
//            System.exit(0);
//        }

        double matchErr = 0;
        double a = 0, b = 0;

        matchedPeaks = new ArrayList<PeakClass>();


        int N = b_ions.size() + y_ions.size() ;

        // The same observed peak might be matchable to multiple theoretical ions for this peptide.
        // This HashMap is used to keep track of the best possible match for an observed peak
        // k = observed peak's m/z value
        // v = the best match for this observed peak based upon the abs(mz_dist)
        TDoubleObjectHashMap<PeakClass> bestMatchMap = new TDoubleObjectHashMap<PeakClass>(N);


        // Try to match y-ions first, they tend to be of higher intensity which is what we want
        for(String theo_ion : y_ions.keySet()) {
            double theo_mz  = y_ions.get(theo_ion);

            if(globals.ms2tol_units == constants.PPM_UNITS) {
                double ppmErr = globals.ms2tol / constants.PPM;
                matchErr = theo_mz * ppmErr;
            }
            else matchErr = globals.ms2tol;

            matchErr *= 0.5; // split in half
            a = globals.round_dbl((theo_mz - matchErr), 4);
            b = globals.round_dbl((theo_mz + matchErr), 4);

            // Hold candidate matching peaks here
           ArrayList<PeakClass> cand = new ArrayList<PeakClass>();
            for(int i = 0; i < obsPeakList.N; i++) {
                double obsMZ = globals.round_dbl( obsPeakList.mz[i], 4 );

                if( (obsMZ >= a) && (obsMZ <= b) ) {
                    cand.add( obsPeakList.getPeakClassInstance(i) );
                }
            }

            if(!cand.isEmpty()) {
                // Take the most intense peak as the best match
                Collections.sort(cand, PeakClass.comparator_intensity_hi2low);
                PeakClass pk = cand.get(0);

                pk.matched = true;
                pk.matchedIonStr = theo_ion;
                pk.dist = pk.mz - theo_mz;

                // Check to see if you have already assigned this peak to a theoretical ion
                if(bestMatchMap.containsKey(pk.mz)) {
                    PeakClass oldPK = bestMatchMap.get(pk.mz);
                    if(Math.abs(oldPK.dist) > (Math.abs(pk.dist))) bestMatchMap.put(pk.mz, pk);
                    oldPK = null;
                }
                else bestMatchMap.put(pk.mz, pk);
            }
            cand.clear();
        }


        // Now try to match b-ions
        for(String theo_ion : b_ions.keySet()) {
            double theo_mz  = b_ions.get(theo_ion);

            if(globals.ms2tol_units == constants.PPM_UNITS) {
                double ppmErr = globals.ms2tol / constants.PPM;
                matchErr = theo_mz * ppmErr;
            }
            else matchErr = globals.ms2tol;

            matchErr *= 0.5; // split in half
            a = globals.round_dbl((theo_mz - matchErr), 4);
            b = globals.round_dbl((theo_mz + matchErr), 4);

            // Hold candidate matching peaks here
            ArrayList<PeakClass> cand = new ArrayList<PeakClass>();
            for(int i = 0; i < obsPeakList.N; i++) {
                double obsMZ = globals.round_dbl( obsPeakList.mz[i], 4 );

                if( (obsMZ >= a) && (obsMZ <= b) ) {
                    cand.add( obsPeakList.getPeakClassInstance(i) );
                }
            }

            if(!cand.isEmpty()) {
                // Take the most intense peak as the best match
                Collections.sort(cand, PeakClass.comparator_intensity_hi2low);
                PeakClass pk = cand.get(0);

                pk.matched = true;
                pk.matchedIonStr = theo_ion;
                pk.dist = pk.mz - theo_mz;

                // Check to see if you have already assigned this peak to a theoretical ion
                if(bestMatchMap.containsKey(pk.mz)) {
                    PeakClass oldPK = bestMatchMap.get(pk.mz);
                    if(Math.abs(oldPK.dist) > (Math.abs(pk.dist))) bestMatchMap.put(pk.mz, pk);
                    oldPK = null;
                }
                else bestMatchMap.put(pk.mz, pk);
            }
            cand.clear();
        }

        // MM holds the best match for each observed peak.
        for(double mz : bestMatchMap.keys()) matchedPeaks.add( bestMatchMap.get(mz) );
        bestMatchMap.clear();
        bestMatchMap = null;
    }
}
