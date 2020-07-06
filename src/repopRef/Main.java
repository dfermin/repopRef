package repopRef;

import htsjdk.variant.vcf.*;
import tech.tablesaw.api.ColumnType;
import tech.tablesaw.api.Table;
import tech.tablesaw.io.csv.CsvReadOptions;

import java.io.*;
import java.util.*;

public class Main {

    static public String inputPath = null;
    static public String variantListFileName = null;
    static public String outputDir = null;
    static public int stepNumber;
    static public int totNumberOfFiles;
    static public HashSet<String> obsVariants;
    static public Table allObservedVariants;
    static public ArrayList<String> headerLines;

    public static void main(String[] args) {

        parseCommandLineArgs(args);

        System.err.println("Source folder: " + inputPath);
        File dir = new File(inputPath);

        totNumberOfFiles = 0;
        for(File fn : dir.listFiles()) {
            if(fn.getName().endsWith(".vcf.gz")) totNumberOfFiles++;
        }

        switch (stepNumber) {
            case 1:
                try {
                    runStep1(dir);
                } catch (IOException e) {
                    e.printStackTrace();
                }
                break;

            case 2:
                try {
                    runStep2(dir);
                } catch (IOException e) {
                    e.printStackTrace();
                }
                break;

            default:
                break;
                // nothing
        }

        System.exit(1);

        // Write the new variants to disk
        //writeVariants();
        System.err.println("\nDone!\n");
    }


    // This function creates the new VCF files with the reference calls put back in
    private static void runStep2(File dir) throws IOException {

        readInVariantList(); // read in all the observed variants from the file created in step 1

        int ctr = 1;
        for(File fn : dir.listFiles()) {
            if(fn.getName().endsWith(".vcf.gz")) {
                System.err.println("Step 2: " + ctr + " of " + totNumberOfFiles  + ": " + fn.getName());
                VCFClass curVCF = new VCFClass(fn);
                curVCF.backfillRefCalls(allObservedVariants, outputDir);
                ctr++;
            }
        }
    }

    private static void readInVariantList() throws IOException {
        File inF = new File(variantListFileName);
        if( !inF.exists() ) {
            System.err.println("\nERROR: Unable to find " + variantListFileName + "\n");
            System.exit(4);
        }

        /*==========
        Using tablesaw (https://jtablesaw.github.io/tablesaw/gettingstarted)
        This reads in the file produced by step 1 and sorts it by chromosome and position
         ===========*/
        ColumnType[] types = {ColumnType.CATEGORY, ColumnType.INTEGER, ColumnType.CATEGORY, ColumnType.CATEGORY};
        CsvReadOptions.CsvReadOptionsBuilder builder = CsvReadOptions.builder(variantListFileName)
                .separator('\t')
                .header(true)
                .columnTypes(types);
        CsvReadOptions options = builder.build();
        allObservedVariants = Table.read().csv(options).sortOn("chrom", "pos");

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

/*
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
*/

    // This function records all the variants observed across all the VCF files in the input folder
    public static void runStep1(File dir) throws IOException {
        obsVariants = new HashSet<>();

        // Iterate over each VCF file in 'dir'
        int ctr = 1;
        for(File fn : dir.listFiles()) {
            // Record every distinct variant observed across all the files
            if(fn.getName().endsWith(".vcf.gz")) {
                System.err.println("Step 1: " + ctr + " of " + totNumberOfFiles + ": " + fn.getName());
                VCFClass curVCF = new VCFClass(fn);
                curVCF.recordObsVariants();
                obsVariants.addAll( curVCF.getAllObsVariants() );
                ctr++;
            }
        }

        // Create the name for the temporary output file
        String outFname = dir.getName() + ".allObsVariants.step1";
        File outF = new File(outFname);
        FileWriter fw = new FileWriter(outF);
        BufferedWriter bw = new BufferedWriter(fw);
        bw.write("chrom\tpos\tref\talt\n");
        for(String k : obsVariants) {
            bw.write(k + "\n");
        }
        bw.close();
        System.err.println("\n" + obsVariants.size() + " variants written to " + outFname + "\n");
        obsVariants.clear();
    }


    // This function combines the variant calls from mulitple files and stores them
    // Into the static variable allVariants
//    private static void recordVariants(VCFClass curVCF, boolean firstIter) {
//        if(firstIter) allVariants = new HashMap<>();
//        HashMap<String, VariantClass> curVCFgenotypes = curVCF.getVariantMap();
//        for(String k : curVCFgenotypes.keySet()) {
//            VariantClass vc = curVCFgenotypes.get(k);
//
//            // Check to see if this variant is already in the 'allVariants' map
//            // If it is, just add the current sample's genotype to the genotypeMap
//            if( allVariants.containsKey(k) ) {
//                allVariants.get(k).addGenotype( vc.getPatientGenotype() );
//            }
//            else {
//                // Otherwise this is a new variant entry for 'allVariants'
//                allVariants.put(k, vc);
//            }
//        }
//    }

    public static void parseCommandLineArgs(String[] args) {

        if(args.length == 0) {
            System.err.println("\nUSAGE: java -jar repopRef.jar -s <step number> -i <folder of input VCFs> -r <output from step1>");
            System.err.println("\tStep 1: Create text file of all variant coordinates found among all VCF files");
            System.err.println("\tStep 2: Repopulate the reference calls in all VCF files.");
            System.err.println("\t        You need the output file from Step 1 for Step 2\n");
            System.exit(1);
        }

        stepNumber = 0;
        for(int i = 0; i < args.length - 1; i++) {
            int j = i + 1;
            if(args[i].equals("-i")) inputPath = args[j];
            if(args[i].equals("-s")) stepNumber = Integer.valueOf(args[j]);
            if(args[i].equals("-r")) variantListFileName = args[j];
            if(args[i].equals("-o")) outputDir = args[j];
        }

        if(stepNumber == 0) {
            System.err.println("\nERROR: You must specify a step number (1,2 or 3)\n");
            System.exit(1);
        }

        if( (stepNumber == 2) && (variantListFileName == null) ) {
            System.err.println("\nERROR: Step 2 requires the output file produced by step 1.");
            System.err.println("Please either run step 1 or provide the file with the '-r' argument\n");
            System.exit(2);
        }

        if( (stepNumber == 2) && (outputDir == null) ) {
            System.err.println("\nERROR: Step 2 requires an output directory.");
            System.err.println("Please create the output folder and provide it as the '-o' argument\n");
            System.exit(3);
        }
    }
}
