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
    static public int stepNumber;
    static public HashSet<String> obsVariants;
    static public Table allObservedVariants;

    public static void main(String[] args) {

        parseCommandLineArgs(args);

        switch (stepNumber) {
            case 1:
                try {
                    runStep1();
                } catch (IOException e) {
                    e.printStackTrace();
                }
                break;

            case 2:
                try {
                    runStep2();
                } catch (IOException e) {
                    e.printStackTrace();
                }
                break;

            default:
                break;
                // nothing
        }

        System.err.println("\nDone!\n");
    }


    // This function creates the new VCF files with the reference calls put back in
    private static void runStep2() throws IOException {

        readInVariantList(); // read in all the observed variants from the file created in step 1
        File fn = new File(inputPath);
        System.err.println("Step 2: " + fn.getName());
        VCFClass curVCF = new VCFClass(fn);
        curVCF.backfillRefCalls(allObservedVariants);
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
    public static void runStep1() throws IOException {
        obsVariants = new HashSet<>();

        File dir = new File(inputPath);
        int totNumberOfFiles = 0;
        for(File fn : dir.listFiles()) {
            if(fn.getName().endsWith(".vcf.gz")) totNumberOfFiles++;
        }

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

        collapseMultiAllelicSites();

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


    // Function to collapse multi-allelic sites recorded into obsVariants
    public static void collapseMultiAllelicSites() {
        HashMap<String, HashSet<String> > ALTmap = new HashMap<>(); // k = variant_loci, v = set of all it's alternative alleles

        // The same variant may have different alternative alleles.
        // This code merges all of the alternative allele strings into a single line.
        for(String line : obsVariants) {
            String[] parts = line.split("\t");
            String k = parts[0] + "\t" + parts[1] + "\t" + parts[2];
            String alt = parts[3];

            if( !ALTmap.containsKey(k) ) {
                ALTmap.put(k, new HashSet<String>());
            }
            ALTmap.get(k).add(alt);
        }

        obsVariants.clear();
        for(String k : ALTmap.keySet()) {
            HashSet<String> all_alts = new HashSet<>();
            for(String s: ALTmap.get(k)) {
                for(String part : s.split(",")) { all_alts.add(part); }
            }
            String merged = k + "\t" + String.join(",", all_alts);
            obsVariants.add(merged);
        }
    }


    public static void parseCommandLineArgs(String[] args) {

        if(args.length == 0) {
            System.err.println("\nUSAGE: java -jar repopRef.jar -s <step number> -i <input> -r <output from step1>");
            System.err.println("\t-i\tInput either a folder of VCF files for step 1 or an individual VCF file for step 2");
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

        if( (stepNumber == 2) && (new File(inputPath).isDirectory()) ) {
            System.err.println("\nERROR: For Step 2 you must provide a single VCF file as input.");
            System.exit(3);
        }
    }
}
