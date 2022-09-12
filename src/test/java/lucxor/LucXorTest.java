package lucxor;

import lucxor.common.PSM;
import lucxor.common.PSMList;
import lucxor.common.PepXML;
import lucxor.utils.Constants;
import org.junit.Before;
import org.junit.Test;
import org.slf4j.LoggerFactory;
import umich.ms.datatypes.index.Index;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.LCMSDataSource;
import umich.ms.fileio.filetypes.mzml.MZMLFile;
import umich.ms.fileio.filetypes.mzxml.MZXMLFile;
import umich.ms.fileio.filetypes.pepxml.PepXmlParser;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.*;

import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import java.io.*;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class LucXorTest {

  private File fileMzML;
  private File filePepXML;
  private final PSMList psmList = new PSMList();

  private static final org.slf4j.Logger log = LoggerFactory.getLogger(LucXorTest.class);


  @Before
  public void setUp() throws Exception {
    URL url = LucXorTest.class.getClassLoader().getResource("tiny.pwiz.1.1.mzML");
    assert url != null;
    fileMzML = new File(url.toURI());

    url = LucXorTest.class.getClassLoader().getResource("interact-ipro.pep.xml");
    assert url != null;
    filePepXML = new File(url.toURI());
  }

  @Test
  public void readMZML() throws FileParsingException {

    SimpleDateFormat fmt= new SimpleDateFormat("HH:mm:ss");
    String ext = afterLastDot(fileMzML.getAbsolutePath()).toLowerCase();
    LCMSDataSource<?> source;
    switch (ext) {
      case "mzml":
        source = new MZMLFile(fileMzML.getAbsolutePath());
        break;
      case "mzxml":
        source = new MZXMLFile(fileMzML.getAbsolutePath());
        break;
      default:
        throw new UnsupportedOperationException("fileMzML not supported");
    }
    // create data structure to hold scans and load all scans
    ScanCollectionDefault scans = new ScanCollectionDefault();
    scans.setDataSource(source);

    Index<?> index = source.parseIndex();
    NavigableMap<Integer, ?> map = index.getMapByNum();
    final int batchSize = 100;

    System.out.printf("Start loading whole fileMzML @ [%s]\n", fmt.format(new Date()));
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
    System.out.printf("Done loading whole fileMzML @ [%s]\n", fmt.format(new Date()));
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

  @Test
  public void readPepXML() throws FileParsingException, IOException {

    Path path = Paths.get(filePepXML.getAbsolutePath());

    try (final FileInputStream fis = new FileInputStream(path.toString())) {

      //Parse metadata

      Iterator<MsmsRunSummary> it = PepXmlParser.parse(fis);
      while (it.hasNext()) {
        MsmsRunSummary runSummary = it.next();
        if(runSummary.getSearchSummary() != null && !runSummary.getSearchSummary().isEmpty()){
          List<SearchSummary> searchSummaryList = runSummary.getSearchSummary();
          PepXML.parseSummaryPTMs(searchSummaryList);
        }

        runSummary.getSpectrumQuery().forEach(query ->{
          long startScan = query.getStartScan();
          String spectrumId = query.getSpectrum();
          int charge = query.getAssumedCharge();
          query.getSearchResult().forEach(results->{
            if(results.getSearchHit().size() > 0)
              log.info("Multiple peptides reported for the same spectrum, LucXor has been tested for that");

            results.getSearchHit().forEach(hit->{
              PSM curPSM = new PSM();
              curPSM.setSpecId(spectrumId);
              curPSM.setScanNum(new Long(startScan).intValue());
              curPSM.setCharge(charge);
              curPSM.setPeptideSequence(hit.getPeptide());
              List<NameValueType> scores = hit.getSearchScore();
              PepXML.parseScores(curPSM, Constants.PEPPROPHET, scores, hit.getAnalysisResult());
              if(hit.getModificationInfo() != null)
                PepXML.addPTMs(curPSM, hit.getModificationInfo());
              if(PepXML.isValidPSM(curPSM))
                psmList.add(curPSM);
            });
          });
        });
      }
    }
  }

  private static boolean advanceReaderToNextRunSummary(XMLStreamReader xsr)
          throws XMLStreamException {
    do {
      if (xsr.next() == XMLStreamConstants.END_DOCUMENT)
        return false;
    } while (!(xsr.isStartElement() && xsr.getLocalName().equals("msms_run_summary")));

    return true;
  }
}
