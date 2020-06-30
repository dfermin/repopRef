package repopRef;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import tech.tablesaw.api.Table;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPOutputStream;

public class VCFClass {
    private VCFHeader headers;
    private VCFFileReader vcfr = null;
    private HashMap<String, VariantClass> variantMap;
    private String origVCFname = null;

    public VCFClass(File fn) {
        File tabixF = new File(fn.getAbsoluteFile() + ".tbi");
        if (!tabixF.exists()) {
            System.err.println("\nERROR: " + tabixF.getName() + " not found!\n");
            System.exit(2);
        }

        vcfr = new VCFFileReader(fn, tabixF, true);
        origVCFname = fn.getName();
    }


    public void recordObsVariants() {
        variantMap = new HashMap<>();

        // Record EACH variant found in this file
        CloseableIterator<VariantContext> it = vcfr.iterator();
        while(it.hasNext()) {
            VariantClass curVariant = new VariantClass(it.next());
            variantMap.put( curVariant.getID(), curVariant );
        }
        it.close();
    }


     public void backfillRefCalls(Table obsVariants, String outputDir) throws IOException {
       // Record the header lines
       headers = new VCFHeader();
       headers = vcfr.getFileHeader();

       // check to make sure the output folder exists
       File outFolder = new File(outputDir);
       if( !outFolder.isDirectory() ) {
           System.err.println("\nERROR: Output folder " + outputDir + " does not exists\n");
           System.exit(5);
       }

       // Extract the sample identifier from the file name
       String patientID = origVCFname.replace(".vcf.gz", "");

       // create the output file
       // From: https://stackoverflow.com/questions/10887828/string-to-gzipoutputstream
       String new_outF = outputDir + File.separator + origVCFname.replace(".vcf.gz", ".backfilledREF.vcf.gz");
       GZIPOutputStream gzOut = new GZIPOutputStream(new FileOutputStream(new_outF));
       BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(gzOut, "UTF-8"));

       ArrayList<String> headerLines = new ArrayList<String>();
       headerLines.add("##fileformat=VCF4.2");
       headerLines.add( headers.getFormatHeaderLine("GT").toString() );
       headerLines.add( headers.getFormatHeaderLine("AD").toString() );
       headerLines.add( headers.getFormatHeaderLine("DP").toString() );
       headerLines.add( headers.getInfoHeaderLine("AC").toString() );
       headerLines.add( headers.getInfoHeaderLine("AF").toString() );
       headerLines.add( headers.getInfoHeaderLine("AN").toString() );
       for(VCFFilterHeaderLine h : headers.getFilterLines()) { headerLines.add( h.toString() ); }
       for(VCFContigHeaderLine h : headers.getContigLines()) { headerLines.add( h.toString() ); }

       headerLines.add( headers.getMetaDataLine("DRAGENCommandLine").toString() );
       headerLines.add( headers.getMetaDataLine("reference").toString() );
       headerLines.add("##repopRefCmd=java -jar repopRef.jar -s 2 -r " + Main.variantListFileName + " -o " + outputDir );
       headerLines.add("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + patientID);

       for(String line : headerLines) {
            bw.append(line);
            bw.newLine();
       }

       for(int i = 0; i < obsVariants.rowCount(); i++) {
           String chrom = obsVariants.get(i, 0);
           int pos = Integer.valueOf( obsVariants.get(i, 1) );
           String ref = obsVariants.get(i, 2);
           String alts = obsVariants.get(i, 3);

           VariantClass curVariant = null;
           CloseableIterator<VariantContext> it = vcfr.query(chrom, pos, pos);
           if(it.hasNext()) { // this happens if the current variant is present in this VCF file.
               VariantContext vc = it.next();
               curVariant = new VariantClass(vc);
               bw.append( curVariant.getLine() );
               bw.newLine();
           } else {
               // If you get here, that means this variant is not found in this VCF file
               curVariant = new VariantClass();
               bw.append( curVariant.getLine_REF(chrom, pos, ref, alts) );
               bw.newLine();
           }
       }

       bw.close();

//       for(String k : obsVariants) {
//           String chrom = k.split(":")[0];
//           int pos = Integer.valueOf( k.split(":")[1]);
//           CloseableIterator<VariantContext> it = vcfr.query(chrom, pos, pos);
//           // If this happens the current VCF file has this variant.
//           // So just report it out as is.
//           if(it.hasNext()) {
//               VariantContext vc = it.next();
//               VariantClass curVariant = new VariantClass(vc);
//               bw.append( curVariant.getLine() );
//               bw.newLine();
//           } else {
//               // If you got here, that means this is a variant that the current sample is a homogyzous reference for
//               //VariantClass curVariant = new VariantClass(chrom, pos);
//           }
//       }
    }

    // Returns a HashSet of unique variant calls from the given VCF file.
    public HashSet<String> getAllObsVariants() {
        HashSet<String> ret = new HashSet<>();

        for(String k : this.variantMap.keySet()) {
            ret.add( this.variantMap.get(k).getStep1_line() );
        }
        return(ret);
    }

    public HashMap<String, VariantClass> getVariantMap() { return(this.variantMap); }

    public VCFHeader getHeaders() { return(this.headers); }


}
