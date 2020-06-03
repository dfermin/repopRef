package repopRef;

import htsjdk.variant.vcf.*;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

public class Main {

    static public String inputPath;
    static public Map<String, VariantClass> allVariants = null;
    static public ArrayList<String> headerLines;

    public static void main(String[] args) {

        if(args.length == 0) {
            System.err.println("\nUSAGE: java -jar repopRef.jar -i <folder of input VCFs> | bgzip -c > STDOUT\n\n");
            System.exit(1);
        }

        for(int i = 0; i < args.length - 1; i++) {
           int j = i + 1;
           if(args[i].equals("-i")) inputPath = args[j];
        }

        System.err.println("Source folder: " + inputPath);
        File dir = new File(inputPath);
        boolean firstIter = true;

        int numFiles = 0;
        for(File fn : dir.listFiles()) {
            if(fn.getName().endsWith(".vcf.gz")) numFiles++;
        }
        int ctr = 1;

        for(File fn : dir.listFiles()) {
            if(fn.getName().endsWith(".vcf.gz")) {
                System.err.print(ctr + " of " + numFiles + ": ");
                VCFClass curVCF = new VCFClass(fn, firstIter);
                ctr++;
                recordVariants(curVCF, firstIter);

                if(firstIter) {
                    recordHeaderLines(curVCF.getHeaders());
                }
                firstIter = false;
            }
        }

        // Write the new variants to disk
        writeVariants();
        System.err.println("\nDone!\n");
    }

    private static void recordHeaderLines(VCFHeader headers) {
        headerLines = new ArrayList<>();
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
        headerLines.add("#repopRefCmd=java -jar repopRef.jar -i " + inputPath );
    }


    private static void writeVariants() {
        System.err.println("\nWriting to STDOUT\n");

        System.out.println("##fileformat=VCFv4.2");
        for(String line : headerLines) { System.out.println("##" + line); }

        // first get all the patient IDs
        TreeSet<String> samples = new TreeSet<>();
        for(VariantClass vc : allVariants.values()) {
            samples.addAll(vc.getAllSampleIDs());
        }

        ArrayList<String> outputHeader = new ArrayList<>();
        outputHeader.add("#CHROM");
        outputHeader.add("POS");
        outputHeader.add("ID");
        outputHeader.add("REF");
        outputHeader.add("ALT");
        outputHeader.add("QUAL");
        outputHeader.add("FILTER");
        outputHeader.add("INFO");
        outputHeader.add("FORMAT");

        for(String s : samples) { outputHeader.add(s); }

        System.out.println(String.join("\t", outputHeader));
        // Iterate over each variant
        int ctr = 0;
        int N = allVariants.size();
        for(VariantClass vc : allVariants.values()) {
            vc.printCalls(samples);
            ctr++;

            if( (ctr % 10000) == 0 ) {
                System.err.println("\rWriting out calls " + ctr + " of " + N);
            }
        }
    }

    // This function combines the variant calls from mulitple files and stores them
    // Into the static variable allVariants
    private static void recordVariants(VCFClass curVCF, boolean firstIter) {
        if(firstIter) allVariants = new HashMap<>();
        HashMap<String, VariantClass> curVCFgenotypes = curVCF.getVariantMap();
        for(String k : curVCFgenotypes.keySet()) {
            VariantClass vc = curVCFgenotypes.get(k);

            // Check to see if this variant is already in the 'allVariants' map
            // If it is, just add the current sample's genotype to the genotypeMap
            if( allVariants.containsKey(k) ) {
                allVariants.get(k).addGenotype( vc.getPatientGenotype() );
            }
            else {
                // Otherwise this is a new variant entry for 'allVariants'
                allVariants.put(k, vc);
            }
        }
    }

}
