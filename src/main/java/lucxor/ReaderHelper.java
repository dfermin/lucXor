package lucxor;

import gnu.trove.map.hash.TIntObjectHashMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Map.Entry;
import java.util.TreeMap;
import umich.ms.datatypes.LCMSData;
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scancollection.IScanCollection;
import umich.ms.datatypes.scancollection.ScanIndex;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.mzml.MZMLFile;
import umich.ms.fileio.filetypes.mzxml.MZXMLFile;
import umich.ms.fileio.filetypes.xmlbased.AbstractXMLBasedDataSource;

/**
 * @author Dmitry Avtonomov
 */
public class ReaderHelper {

  private ReaderHelper() {}

  /**
   * Reads in spectral data from mzML files.
   */
  public static void read_mzML(TreeMap<String, TreeMap<Integer, PSM>> mapFilesToPsms, int numThreads) throws FileParsingException {

    // Iterate over the file names
    readRawXmlFiles(mapFilesToPsms, numThreads);
  }

  private static void readRawXmlFiles(TreeMap<String, TreeMap<Integer, PSM>> mapFilesToPsms, int numThreads)
      throws FileParsingException {
    for (String fn : mapFilesToPsms.keySet()) {
      String fileName = new File(fn).getName();
      System.err.printf("\n%s: ...", fileName); // beginning of info line

      final TreeMap<Integer, PSM> psmsOfInterest = mapFilesToPsms.get(fn);
      Integer threads = numThreads >= 0 ? numThreads : null;

      // read in the mzML file
      LCMSData lcmsData = null;
      try (AbstractXMLBasedDataSource<?, ?> source = getRawDataSource(fn)) {
        source.setNumThreadsForParsing(threads);
        source.setParsingTimeout(60L); // 1 minute before it times out trying to read a file
        lcmsData = new LCMSData(source);

        lcmsData.load(LCMSDataSubset.MS2_WITH_SPECTRA);
        IScanCollection scans = lcmsData.getScans();
        ScanIndex ms2ScanIndex = scans.getMapMsLevel2index().get(2);

        for (Entry<Integer, PSM> psmEntry : psmsOfInterest.entrySet()) {
          final Integer scanNum = psmEntry.getKey();
          final PSM psm = psmEntry.getValue();
          final IScan scan = ms2ScanIndex.getNum2scan().get(scanNum);
          if (scan == null) {

            continue;
          }
          SpectrumClass x = new SpectrumClass(scan.getSpectrum().getMZs(),
              scan.getSpectrum().getIntensities());
          psm.recordSpectra(x);
        }
      } finally {
        if (lcmsData != null) lcmsData.releaseMemory();
      }
      System.err.printf("\r%s: Done\n", fileName); // beginning of info line
    }
  }

  static void read_mzXML(TreeMap<String, TreeMap<Integer, PSM>> mapFilesToPsms, int numThreads) throws FileParsingException {

    // Iterate over the file names
    readRawXmlFiles(mapFilesToPsms, numThreads);
  }

  private static AbstractXMLBasedDataSource<?, ?> getRawDataSource(String fileName) {
    if (fileName.toLowerCase().endsWith(".mzml"))
      return new MZMLFile(fileName);
    if (fileName.toLowerCase().endsWith(".mzxml"))
      return new MZXMLFile(fileName);

    throw new IllegalArgumentException(String.format("Unknown input file format, file name must "
        + "end with .mzml or .mzxml"));
  }

  static TIntObjectHashMap<SpectrumClass> read_mgf(String specFile)
      throws IOException {
    TIntObjectHashMap<SpectrumClass> ret = new TIntObjectHashMap<>();

    File mgf = new File(specFile);
    BufferedReader br = new BufferedReader(new FileReader(mgf));
    String line;
    int scanNum = 0;
    SpectrumClass S;
    ArrayList<Double> mzAL = null, intensityAL = null;

    while ((line = br.readLine()) != null) {

      if (line.length() < 2) {
        continue;
      }

      if (line.startsWith("END IONS")) {
        if ((null != mzAL) && (!mzAL.isEmpty())) {

          int N = mzAL.size();
          double[] mz = new double[N];
          double[] I = new double[N];
          for (short i = 0; i < N; i++) {
            mz[i] = mzAL.get(i);
            I[i] = intensityAL.get(i);
          }
          S = new SpectrumClass(mz, I);
          ret.put(scanNum, S);

          mzAL = null;
          intensityAL = null;
        }
        scanNum = 0;
      }

      if (line.startsWith("BEGIN IONS")) {
        mzAL = new ArrayList<>();
        intensityAL = new ArrayList<>();
      }

      if (line.startsWith("TITLE=")) {
        int i = line.indexOf(".") + 1;
        int j = line.indexOf(".", i);
        String s = line.substring(i, j);
        scanNum = Integer.valueOf(s);
      }

      if (line.startsWith("CHARGE") || line.startsWith("PEPMASS")) {
        continue;
      }

      if (Character.isDigit(line.charAt(0))) {
        String[] ary = line.split("\\s+");
        double mz = MathHelper.roundDouble(Double.valueOf(ary[0]), 8);
        double I = MathHelper.roundDouble(Double.valueOf(ary[1]), 8);
        mzAL.add(mz);
        intensityAL.add(I);
      }
    }

    return ret;
  }
}
