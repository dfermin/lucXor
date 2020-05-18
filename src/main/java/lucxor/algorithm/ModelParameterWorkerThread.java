package lucxor.algorithm;

import lucxor.common.PSM;

/**
 * Created by dfermin on 5/15/14.
 */
public class ModelParameterWorkerThread implements Runnable {

    private final PSM curPSM; // the PSM modeling data will be acquired from

    // Default constructor for this class
    public ModelParameterWorkerThread(PSM externalPSM) {
        this.curPSM = externalPSM;
        // index of this PSM in LucXorConfiguration.psmList
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
