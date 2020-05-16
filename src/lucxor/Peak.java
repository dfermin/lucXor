/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import java.util.Comparator;

/**
 * The 'Cloneable' allows you to make copies of Peak objects
 *
 * @author dfermin
 */
public class Peak {

  private double mz;
  private double rawIntensity; // original intensity
  private double relIntensity; // intensity normalized to be between 0 - 100
  private double normIntensity; // normalized intensity

  private double distScore;
  private double intensityScore;
  private double dist;

  private double score;

  private boolean matched;
  private String matchedIonStr;
  private double matchedIonMZ; // the m/z value of the matched ion


  // Comparator to sort Peak objects based upon relative intensity from high to low
  public static Comparator comparator_RAW_intensity_hi2low = (Comparator<Peak>) (o1, o2) -> Double
    .compare(o2.rawIntensity, o1.rawIntensity);

  // Comparator to sort Peak objects based upon relative intensity from high to low
  public static Comparator comparator_intensity_hi2low = (Comparator<Peak>) (o1, o2) -> Double
    .compare(o2.relIntensity, o1.relIntensity);

  // Comparator to sort Peak objects based upon relative intensity from low to high
  public static Comparator comparator_intensity_low2hi = Comparator.comparingDouble(o -> ((Peak)o).relIntensity);

  // Comparator to sort Peak objects based upon mz values from low to high
  public static Comparator comparator_mz = Comparator.comparingDouble(o -> ((Peak)o).mz);

  // Comparator to sort Peak objects based upon mz dist values from low to high
  public static Comparator comparator_mz_abs_dist = Comparator.comparingDouble(o -> Math.abs(((Peak)o).dist));

  // default constructor
  Peak(double x, double y) { mz = x; rawIntensity = y; }

  public Peak() {
    mz = 0;
    relIntensity = 0;
    distScore = 0;
    dist = 0;
    intensityScore = 0;
    matchedIonMZ = 0;
    score = 0;
    matched = false;
  }
  // Constructor used for copying the variables in one instance of a peak class to another
  public Peak(Peak other) {
    this.mz = other.mz;
    this.rawIntensity = other.rawIntensity;
    this.relIntensity = other.relIntensity;
    this.normIntensity = other.normIntensity;
    this.distScore = other.distScore;
    this.intensityScore = other.intensityScore;
    this.dist = other.dist;
    this.score = other.score;
    this.matched = other.matched;
    this.matchedIonStr = other.matchedIonStr;
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


  /**************************************************************************
   * Function to use for comparison between Peak objects to remove
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
    final Peak other = (Peak) obj;
    return Double.doubleToLongBits(this.mz) == Double.doubleToLongBits(other.mz);
  }

  public double getMz() {
    return mz;
  }

  public void setMz(double mz) {
    this.mz = mz;
  }

  public double getRawIntensity() {
    return rawIntensity;
  }

  public void setRawIntensity(double rawIntensity) {
    this.rawIntensity = rawIntensity;
  }

  public double getRelIntensity() {
    return relIntensity;
  }

  public void setRelIntensity(double relIntensity) {
    this.relIntensity = relIntensity;
  }

  public double getNormIntensity() {
    return normIntensity;
  }

  public void setNormIntensity(double normIntensity) {
    this.normIntensity = normIntensity;
  }

  public double getDistScore() {
    return distScore;
  }

  public void setDistScore(double distScore) {
    this.distScore = distScore;
  }

  public double getIntensityScore() {
    return intensityScore;
  }

  public void setIntensityScore(double intensityScore) {
    this.intensityScore = intensityScore;
  }

  public double getDist() {
    return dist;
  }

  public void setDist(double dist) {
    this.dist = dist;
  }

  public double getScore() {
    return score;
  }

  public void setScore(double score) {
    this.score = score;
  }

  public boolean isMatched() {
    return matched;
  }

  public void setMatched(boolean matched) {
    this.matched = matched;
  }

  public String getMatchedIonStr() {
    return matchedIonStr;
  }

  public void setMatchedIonStr(String matchedIonStr) {
    this.matchedIonStr = matchedIonStr;
  }

  public double getMatchedIonMZ() {
    return matchedIonMZ;
  }

  public void setMatchedIonMZ(double matchedIonMZ) {
    this.matchedIonMZ = matchedIonMZ;
  }
}
