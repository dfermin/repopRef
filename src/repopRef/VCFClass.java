package repopRef;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.util.HashMap;

public class VCFClass {
    private VCFHeader headers;
    private HashMap<String, VariantClass> variantMap;

    public VCFClass(File fn, boolean firstIter) {
        File tabixF = new File(fn.getAbsoluteFile() + ".tbi");
        if(!tabixF.exists()) {
            System.err.println("\nERROR: " + tabixF.getName() + " not found!\n");
            System.exit(2);
        }

        System.err.println("Reading " + fn.getName());
        VCFFileReader vcfr = new VCFFileReader(fn, tabixF, true);

        // Record the header lines
        if(firstIter) {
            headers = new VCFHeader();
            headers = vcfr.getFileHeader();
        }

        // Record EACH variant found in this file
        variantMap = new HashMap<>();
        CloseableIterator<VariantContext> it = vcfr.iterator();
        while(it.hasNext()) {
            VariantContext vc = it.next();
            VariantClass v = new VariantClass(vc);
            variantMap.put(v.getID(), v);
        }
        it.close();
    }


    public HashMap<String, VariantClass> getVariantMap() { return(this.variantMap); }

    public VCFHeader getHeaders() { return(this.headers); }

}
