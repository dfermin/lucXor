package lucxor.pride;

import lucxor.LucXorTest;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.net.URL;


public class PrideJsonRandomAccessTest {

    private File prideJsonFile;
    private PrideJsonRandomAccess randomAccess;

    @Before
    public void setUp() throws Exception {
        URL url = LucXorTest.class.getClassLoader().getResource("PRD000909_54e6c552385ce4c9b77ec3985019f6d977297bdd_ArchiveSpectrum.json");
        assert url != null;
        prideJsonFile = new File(url.toURI());
        randomAccess = new PrideJsonRandomAccess(prideJsonFile.getAbsolutePath());
    }

    @Test
    public void testParseIndex() throws IOException {
        randomAccess.parseIndex();
        Assert.assertEquals(38706, randomAccess.getNumArchiveSpectra());
    }

    @Test
    public void testReadArchiveSpectrum() throws IOException {
        randomAccess.parseIndex();
        String usi = "mzspec:PRD000909:tio2_1mg:scan:94780:SQEPVTLDFLDAELENDIK/2";
        Assert.assertEquals(usi, randomAccess.readArchiveSpectrum(usi).getUsi());

    }
}