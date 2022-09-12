package lucxor.pride;

import com.fasterxml.jackson.databind.ObjectMapper;
import uk.ac.ebi.pride.archive.dataprovider.data.spectra.ArchiveSpectrum;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * {@link PrideJsonRandomAccess} is a reader of the
 * {@link uk.ac.ebi.pride.archive.dataprovider.data.spectra.SummaryArchiveSpectrum}
 * json files.
 *
 * @author ypriverol
 */
public class PrideJsonRandomAccess {

    private final BufferedRandomAccessFile raf;

    private final Map<String, Long> index;
    private final ObjectMapper objectMapper;

    public PrideJsonRandomAccess(String fileAbsolutePath) throws IOException {
        this.raf = new BufferedRandomAccessFile(fileAbsolutePath, "r", 1024 * 100);
        this.index = new HashMap<>();
        this.objectMapper = new ObjectMapper();
    }

    /**
     * Create an index of all the spectra within an {@link ArchiveSpectrum} json file.
     * The index will be a Map with usis (identifiers) as index and values the pointer
     * where the specific spectrum starts. The index is assuming that the ArchiveSpectrum
     * file is one line json representation.
     * @throws IOException
     */
    public void parseIndex() throws IOException {

        String line;
        long pos = raf.getFilePointer();

        while( (line = raf.getNextLine()) != null){
            ArchiveSpectrum spectrum = objectMapper.readValue(line, ArchiveSpectrum.class);
            index.put(spectrum.getUsi(), pos);
            pos = raf.getFilePointer();
        }
    }

    /**
     * This class returns a specific {@link ArchiveSpectrum} for a given usi by
     * doing a lookup in the index table.
     * @param usi identifier of the spectrum
     * @return ArchiveSpectrum
     * @throws IOException
     */
    public ArchiveSpectrum readArchiveSpectrum(String usi) throws IOException {
        if(index.containsKey(usi)){
            long pos = index.get(usi);
            raf.seek(pos);
            return  objectMapper.readValue(raf.readLine(), ArchiveSpectrum.class);
        }
        return null;
    }

    /**
     * GEt the number of spectra in the {@link ArchiveSpectrum} json file.
     * @return total number fo spectra
     */
    public int getNumArchiveSpectra(){
        return index.size();
    }

    /**
     * Return all the spectra from the Json file
     * @return List of usis
     */
    public Set<String> getKeys(){
        return index.keySet();
    }
}
