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


     public void backfillRefCalls(Table obsVariants) throws IOException {
         // Record the header lines
         HashSet<String> seenVariants = new HashSet<>(); // keep track of variants that you've done in this VCF file
         headers = new VCFHeader();
         headers = vcfr.getFileHeader();

         // Extract the sample identifier from the file name
         String patientID = origVCFname.replace(".vcf.gz", "");

         ArrayList<String> headerLines = new ArrayList<String>();
         headerLines.add("##fileformat=VCF4.2");
         headerLines.add(headers.getFormatHeaderLine("GT").toString());
         headerLines.add(headers.getFormatHeaderLine("AD").toString());
         headerLines.add(headers.getFormatHeaderLine("DP").toString());
         headerLines.add(headers.getInfoHeaderLine("AC").toString());
         headerLines.add(headers.getInfoHeaderLine("AF").toString());
         headerLines.add(headers.getInfoHeaderLine("AN").toString());
         for (VCFFilterHeaderLine h : headers.getFilterLines()) {
             headerLines.add(h.toString());
         }
         for (VCFContigHeaderLine h : headers.getContigLines()) {
             headerLines.add(h.toString());
         }

         headerLines.add(headers.getMetaDataLine("DRAGENCommandLine").toString());
         headerLines.add(headers.getMetaDataLine("reference").toString());
         headerLines.add("##repopRefCmd=java -jar repopRef.jar -s 2 -r " + Main.variantListFileName + " -i " + Main.inputPath);
         headerLines.add("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + patientID);

         for (String line : headerLines) {
             System.out.println(line);
         }

         for (int i = 0; i < obsVariants.rowCount(); i++) {
             String chrom = obsVariants.get(i, 0);
             int pos = Integer.valueOf(obsVariants.get(i, 1));
             String ref = obsVariants.get(i, 2);
             String alts = obsVariants.get(i, 3);

             String k = chrom + ":" + String.valueOf(pos);
             if( !seenVariants.contains(k) ) {
                 VariantClass curVariant = null;
                 CloseableIterator<VariantContext> it = vcfr.query(chrom, pos, pos);
                 if (it.hasNext()) { // this happens if the current variant is present in this VCF file.
                     VariantContext vc = it.next();
                     curVariant = new VariantClass(vc);
                     System.out.println(curVariant.getLine());
                 } else {
                     // If you get here, that means this variant is not found in this VCF file
                     curVariant = new VariantClass();
                     System.out.println(curVariant.getLine_REF(chrom, pos, ref, alts));
                 }
                 // We are done with this variant so move on
                 seenVariants.add(k);
             }
         }
     }


    // Returns a HashSet of unique variant calls from the given VCF file.
    public HashSet<String> getAllObsVariants() {
        HashSet<String> ret = new HashSet<>();

        HashMap<String, HashSet<String> > ALTmap = new HashMap<>(); // k = variant_loci, v = set of all it's alternative alleles

        // The same variant may have different alternative alleles.
        // This code merges all of the alternative allele strings into a single line.
        for(String k : this.variantMap.keySet()) {
            VariantClass vc = this.variantMap.get(k);
            String alt = vc.getALT();
            if( !ALTmap.containsKey(k) ) {
                ALTmap.put(k, new HashSet<String>());
            }
            ALTmap.get(k).add(alt);
        }

        for(String k : ALTmap.keySet()) {
            HashSet<String> all_alts = new HashSet<>();
            for(String s : ALTmap.get(k)) {
                for(String part : s.split(",")) { all_alts.add(part); }
            }
            String ref = this.variantMap.get(k).getREF();
            String merged_alts = k.replace(":", "\t") + "\t" + ref + "\t" + String.join(",", all_alts);
            ret.add(merged_alts);
        }
        return(ret);
    }


    public HashMap<String, VariantClass> getVariantMap() { return(this.variantMap); }

    public VCFHeader getHeaders() { return(this.headers); }


}
