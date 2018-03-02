package lucxor;

/**
 * Created by dfermin on 5/15/14.
 */
class ModelParameterWorkerThread implements Runnable {

  private final PSM curPSM; // the PSM modeling data will be acquired from
  private final int jobIndex;

  // Default constructor for this class
  public ModelParameterWorkerThread(PSM externalPSM, int i) {
    this.curPSM = externalPSM;
    jobIndex = i;
  }

  @Override
  public void run() {

    try {
      synchronized (this) {
        curPSM.generatePermutations(0); // generate real permutations
        curPSM.matchAllPeaks();
      }
    } catch (Exception e) {
      e.printStackTrace();
      System.exit(-1);
    }
  }

}
