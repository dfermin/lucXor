/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucXor;

import java.util.Comparator;

/**
 *
 * @author dfermin
 */

// The 'Cloneable' allows you to make copies of PeakClass objects
public class PeakClass {
	double mz;
	double raw_intensity;
	double rel_intensity; // intensity normalized to be between 0 - 100
	double norm_intensity; // normalized intensity
	
	double distScore;
	double intensityScore;
	double dist;
	
	double score;
	
	boolean matched;
	String matchedIonStr;
	double matchedIonMZ; // the m/z value of the matched ion


	// Comparator to sort PeakClass objects based upon relative intensity from high to low
	public static Comparator comparator_RAW_intensity_hi2low = new Comparator<PeakClass>() {
			@Override
			public int compare(PeakClass o1, PeakClass o2) {
					return Double.compare(o2.raw_intensity, o1.raw_intensity);
			}
	};
	
	// Comparator to sort PeakClass objects based upon relative intensity from high to low
	public static Comparator comparator_intensity_hi2low = new Comparator<PeakClass>() {
			@Override
			public int compare(PeakClass o1, PeakClass o2) {
					return Double.compare(o2.rel_intensity, o1.rel_intensity);
			}
	};

	// Comparator to sort PeakClass objects based upon relative intensity from low to high
	public static Comparator comparator_intensity_low2hi = new Comparator<PeakClass>() {
			@Override
			public int compare(PeakClass o1, PeakClass o2) {
					return Double.compare(o1.rel_intensity, o2.rel_intensity);
			}
	};

	// Comparator to sort PeakClass objects based upon mz values from low to high
	public static Comparator comparator_mz = new Comparator<PeakClass>() {
			@Override
			public int compare(PeakClass o1, PeakClass o2) {
                return Double.compare(o1.mz, o2.mz);
			}
	};

	// Comparator to sort PeakClass objects based upon mz dist values from low to high
	public static Comparator comparator_mz_abs_dist = new Comparator<PeakClass>() {
			@Override
			public int compare(PeakClass o1, PeakClass o2) {
					return Double.compare(Math.abs(o1.dist), Math.abs(o2.dist));
			}
	};


	
	/**************************************************************************
	 * Function to use for comparison between PeakClass objects to remove
	 * duplicates based upon the 'mz' variable in this class.
	 */
	@Override
	public boolean equals(Object obj) {
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final PeakClass other = (PeakClass) obj;
		if (Double.doubleToLongBits(this.mz) != Double.doubleToLongBits(other.mz)) {
			return false;
		}
		return true;
	}
	/**************************************************************************/
	

	// default constructor
	PeakClass(double x, double y) { mz = x; raw_intensity = y; }



	public PeakClass() {
			mz = 0;
			rel_intensity = 0;
			distScore = 0;
			dist = 0;
			intensityScore = 0;
			matchedIonMZ = 0;
			score = 0;
			matched = false;
	}



    // Constructor used for copying the variables in one instance of a peak class to another
    public PeakClass(PeakClass other) {
        this.mz = other.mz;
        this.raw_intensity = other.raw_intensity;
        this.rel_intensity = other.rel_intensity;
        this.norm_intensity = other.norm_intensity;
        this.distScore = other.distScore;
        this.intensityScore = other.intensityScore;
        this.dist = other.dist;
        this.score = other.score;
        this.matched = other.matched;
        this.matchedIonStr = new String( other.matchedIonStr );
        this.matchedIonMZ = other.matchedIonMZ;
    }






    // clear any previously assigned score data to a peak.
	public void clear() {
		distScore = 0;
		dist = 0;
		score = 0;
		intensityScore = 0;
		matchedIonMZ = 0;
		matchedIonStr = "";
		matched = false;
	}

}
