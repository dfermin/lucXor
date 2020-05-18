package lucxor;

import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.index.Index;
import umich.ms.datatypes.index.IndexElement;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scan.StorageStrategy;
import umich.ms.datatypes.scancollection.IScanCollection;
import umich.ms.datatypes.scancollection.ScanIndex;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.LCMSDataSource;
import umich.ms.fileio.filetypes.mzml.MZMLFile;
import umich.ms.fileio.filetypes.mzml.MZMLIndex;
import umich.ms.fileio.filetypes.mzxml.MZXMLFile;

import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;

import static org.junit.Assert.*;

public class LucXorTest {

  private String file;

  @Before
  public void setUp() throws Exception {
    //file = "/Volumes/yasset_data/work/luxifor/020320_TiO2_Phospho.mzML";
    file = "D:\\ms-data\\batmass-io-errors\\HCC1806_F-2_E-1_M22_5ul_I-1_01.mzML";
  }

  @Test
  @Ignore
  public void readMZML() throws FileParsingException {

    SimpleDateFormat fmt= new SimpleDateFormat("HH:mm:ss");
    String ext = afterLastDot(file).toLowerCase();
    LCMSDataSource<?> source;
    switch (ext) {
      case "mzml":
        source = new MZMLFile(file);
        break;
      case "mzxml":
        source = new MZXMLFile(file);
        break;
      default:
        throw new UnsupportedOperationException("file not supported");
    }
    // create data structure to hold scans and load all scans
    ScanCollectionDefault scans = new ScanCollectionDefault();
    scans.setDataSource(source);

    Index<?> index = source.parseIndex();
    NavigableMap<Integer, ?> map = index.getMapByNum();
    final int batchSize = 100;

    System.out.printf("Start loading whole file @ [%s]\n", fmt.format(new Date()));
    long timeLo = System.nanoTime();

    Iterator<? extends Map.Entry<Integer, ?>> it = map.entrySet().iterator();
    while (it.hasNext()) {
      List<Integer> scanNums = new ArrayList<>();
      scanNums.add(it.next().getKey());
      while (it.hasNext() && scanNums.size() < batchSize) {
        scanNums.add(it.next().getKey());
      }
      List<IScan> parsed = source.parse(scanNums);
      System.out.printf("Read batch of %d scans, starting with: %s\n", parsed.size(), parsed.get(0).toString());
    }

    long timeHi = System.nanoTime();
    System.out.printf("Done loading whole file @ [%s]\n", fmt.format(new Date()));
    System.out.printf("Loading took %.1fs", (timeHi - timeLo)/1e9f);
    // data index, can be used to locate scans by numbers or retention times at different ms levels

    ExecutorService exec = Executors
      .newFixedThreadPool(Runtime.getRuntime().availableProcessors());
    for (final Map.Entry<Integer, IScan> kv : scans.getMapNum2scan().entrySet()) {
      exec.submit(() -> {
        IScan scan = kv.getValue();
        double spectrumIntensitySum = Arrays
          .stream(scan.getSpectrum().getIntensities()).sum();

        System.out.printf("Scan [%s], intensity sum: %.2f\n", scan.toString(),
          spectrumIntensitySum);
      });
    }
    System.out.printf("\n\nLoading took %.1fs", (timeHi - timeLo)/1e9f);
   }

  private String afterLastDot(String s) {
      int last = s.lastIndexOf('.');
      return last < 0 ? "" : s.substring(last+1);
  }
}
