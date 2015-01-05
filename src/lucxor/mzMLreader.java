package lucxor;

import gnu.trove.map.hash.TIntObjectHashMap;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import java.io.File;
import java.io.IOException;
import java.util.zip.Inflater;

import java.io.ByteArrayOutputStream;
import java.util.zip.DataFormatException;




/**
 * Created by dfermin on 4/30/14.
 */
public class mzMLreader extends DefaultHandler {

    String srcMzML;
    String mzArrayString;
    String intensityArrayString;
    String binaryDataType;
    StringBuilder tmpValue = null;
    int precision;
    int numScans;
    int scanNum;
    int msLevel;
    int mz_precision, int_precision;

    Boolean isBase64Encoded = true;
    Boolean mzStringDone, intStringDone; // true means you processed this data
    String compressionType = "none";     // default value, if nothing specified in file - assume no compression


    static TIntObjectHashMap<SpectrumStruct> spectrumMap; // k = scanNumber, v = spectrum
    private SpectrumStruct curSpectrum;


    private class SpectrumStruct {
        double[] mz_;
        double[] intensity_;
    }


    public mzMLreader(String inputF) {

        precision = 0;
        mz_precision = 0;
        int_precision = 0;

        srcMzML = inputF;
        File f = new File(srcMzML);

        if(!f.exists()) {
            System.err.print("\nmzMLreader ERROR: Unable to find '" + srcMzML + "'\n\n");
            System.exit(0);
        }

        spectrumMap = new TIntObjectHashMap<mzMLreader.SpectrumStruct>();

        tmpValue = new StringBuilder();

        SAXParserFactory factory = SAXParserFactory.newInstance();
        SAXParser parser = null;
        try {
            parser = factory.newSAXParser();
            parser.parse(srcMzML, this);
        } catch (ParserConfigurationException | SAXException | IOException e) {
            e.printStackTrace();
        }
    }



    @Override
    public void characters(char[] ac, int start, int length) throws SAXException {
        tmpValue.append(ac, start, length);
    }


    @Override
    public void startElement(String s, String s1, String element, Attributes attr) throws SAXException {

        if(element.equalsIgnoreCase("spectrumList")) {
            numScans = Integer.valueOf( attr.getValue("count") );
        }

        if(element.equalsIgnoreCase("spectrum")) {
            scanNum = Integer.valueOf( attr.getValue("index") ) + 1;
            mzStringDone = false;
            intStringDone = false;
        }

        if(element.equalsIgnoreCase("cvParam")) {
            for(int i = 0; i < attr.getLength(); i++) {
                String attr_name  = attr.getLocalName(i);
                String attr_value = attr.getValue(i);

                if(attr_name.equalsIgnoreCase("accession")) {

                    if(attr_value.equals("MS:1000523")) {
                        precision = 64;
                    }
                    if(attr_value.equals("MS:1000521")) {
                        precision = 32;
                    }

                    if(attr_value.equalsIgnoreCase("MS:1000511")) { // ms-level
                        msLevel = Integer.valueOf(attr.getValue("value"));
                        break;
                    }

                    if(attr_value.equalsIgnoreCase("MS:1000576")) { // compression enabled
                        compressionType = "zlib";
                        break;
                    }

                    if(attr_value.equalsIgnoreCase("MS:1000574")) { // compression disabled
                        compressionType = "none";
                    }

                    if(attr_value.equalsIgnoreCase("MS:1000514")) { // m/z array
                        binaryDataType = "m/z";
                        mzArrayString = "";
                        mz_precision = precision;
                        tmpValue = new StringBuilder();
                        break;
                    }

                    if(attr_value.equalsIgnoreCase("MS:1000515")) { // intensity array
                        binaryDataType = "intensity";
                        intensityArrayString = "";
                        int_precision = precision;
                        tmpValue = new StringBuilder();
                        break;
                    }
                }
            }
        }
    } // end startElement()


    @Override
    public void endElement(String s, String s1, String element) throws SAXException {

        // this conditional actually captures the compressed encoded string
        if(element.equalsIgnoreCase("binary")) {
            if( (mzStringDone == false) && (binaryDataType.equals("m/z")) ) {
                mzArrayString = tmpValue.toString();
                mzStringDone = true;
            }
            else if( (intStringDone == false) && (binaryDataType.equals("intensity")) ) {
                intensityArrayString = tmpValue.toString();
                intStringDone = true;
            }
        }


        if(element.equalsIgnoreCase("spectrum")) { // end of a record
            if(msLevel == 2) {

                double[] mz_ = null;
                double[] int_ = null;
                try {
                    mz_  = decode_string(mzArrayString, mz_precision);
                    int_ = decode_string(intensityArrayString, int_precision);

                    curSpectrum = new SpectrumStruct();

                    curSpectrum.mz_ = new double[ mz_.length ];
                    curSpectrum.intensity_ = new double[ mz_.length ];

                    curSpectrum.mz_ = mz_;
                    curSpectrum.intensity_ = int_;

                    spectrumMap.put(scanNum, curSpectrum);

                    mz_ = null;
                    int_ = null;
                    curSpectrum = null;

                } catch (IOException e) {
                    e.printStackTrace();
                } catch (DataFormatException e) {
                    e.printStackTrace();
                }
            }
        }
    }


    /*******************
     * Function returns an array of doubles that was decoded from the passed string
     * @param peaks
     * @return
     * @throws IOException
     * @throws DataFormatException
     */
    private double[] decode_string(String peaks, int precision) throws IOException, DataFormatException {

        byte[] decoded;

        if(isBase64Encoded) {
            decoded = org.apache.commons.codec.binary.Base64.decodeBase64(peaks);
        } else {
            decoded = peaks.getBytes();
        }


        if(compressionType.equals("zlib")) { // need to decompress the bytes
            decoded = zlibUncompressBuffer(decoded);
        }

        int decodedLen = decoded.length; // in bytes
        int chunkSize = precision / 8;
        int numPeaks = decodedLen / chunkSize;
        double[] retAry = new double[numPeaks];


        if (precision == 32) {

            int asInt;
            float asFloat = 0.0f;
            int offset;

            for (int i = 0; i < numPeaks; i++) {
                offset = i * chunkSize;

                asInt = (  (decoded[offset + 0] & 0xFF)) // zero shift
                        | ((decoded[offset + 1] & 0xFF) << 8)
                        | ((decoded[offset + 2] & 0xFF) << 16)
                        | ((decoded[offset + 3] & 0xFF) << 24);
                asFloat = Float.intBitsToFloat(asInt);
                retAry[i] = asFloat;
            }
        } else if (precision == 64) {
            long asLong;
            double asDouble = 0.0d;
            int offset;

            for (int i = 0; i < numPeaks; i++) {
                offset = i * chunkSize;

                asLong = ( (long) (decoded[offset + 0] & 0xFF)) // zero shift
                        | ((long) (decoded[offset + 1] & 0xFF) << 8)
                        | ((long) (decoded[offset + 2] & 0xFF) << 16)
                        | ((long) (decoded[offset + 3] & 0xFF) << 24)
                        | ((long) (decoded[offset + 4] & 0xFF) << 32)
                        | ((long) (decoded[offset + 5] & 0xFF) << 40)
                        | ((long) (decoded[offset + 6] & 0xFF) << 48)
                        | ((long) (decoded[offset + 7] & 0xFF) << 56);
                asDouble = Double.longBitsToDouble(asLong);

                retAry[i] = asDouble;
            }
        } else {
            throw new IllegalArgumentException("Precision can only be 32 or 64 bits.");
        }

        return retAry;
    }



    /**
     * Inflates zLib compressed byte[].
     * @param compressed zLib compressed bytes
     * @return inflated byte array
     * @throws IOException should never happen, ByteArrayOutputStream is in-memory
     * @throws DataFormatException in case of malformed input byte array
     */
    public static byte[] zlibUncompressBuffer(byte[] compressed) throws IOException, DataFormatException {

        Inflater decompressor = new Inflater();
        decompressor.setInput(compressed);

        ByteArrayOutputStream bos = new ByteArrayOutputStream(compressed.length);
        byte[] buf = new byte[decompressor.getRemaining() * 2];
        try {
            // Decompress the data
            while (decompressor.getRemaining() > 0) {
                int count = decompressor.inflate(buf);
                bos.write(buf, 0, count);
            }

        } finally {
            try {
                bos.close();
            } catch (IOException nope) {
                // This exception doesn't matter, but it totally should not happen
                throw nope;
            }
        }
        decompressor.end();
        byte[] result = bos.toByteArray();
        return result;
    }



    /***********
     * Function returns the number of elements in spectrum arrays
     * @param  sn
     */
    public static int getNumPeaks(int sn) {
        SpectrumStruct S = null;
        int ret = 0;

        if(spectrumMap.containsKey(sn)) {
            S = spectrumMap.get(sn);
            ret = S.mz_.length;
        }

        return ret;
    }



    /************
     * Function to get the m/z array for a given scan number
     * @param sn
     */
    public static double[] getMZ(int sn) {
        SpectrumStruct ret = null;

        if(spectrumMap.containsKey(sn)) {
            ret = spectrumMap.get(sn);
        }
        return ret.mz_;
    }



    /************
     * Function to get the intensity array for a given scan number
     * @param sn
     */
    public static double[] getIntensities(int sn) {
        SpectrumStruct ret = null;

        if(spectrumMap.containsKey(sn)) {
            ret = spectrumMap.get(sn);
        }
        return ret.intensity_;
    }
}
