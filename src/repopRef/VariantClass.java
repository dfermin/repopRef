package repopRef;

import htsjdk.variant.variantcontext.*;
import java.text.DecimalFormat;
import java.util.*;

public class VariantClass {
    private String chrom;
    private long pos;
    private Allele REF;
    private List<Allele> ALT;
    private String uniqID;
    private CommonInfo CI;
    private HashMap<String, Genotype> nonRefCalls = null; // k = sampleID, v = 0/1, 1/0, 1/1
    private HashMap<String, Integer> readDepth = null; // k = sampleID, v = DP for this variant in this sample

    public VariantClass(VariantContext vc) {

       this.chrom = vc.getContig();
       this.pos = vc.getStart();
       this.REF = vc.getReference();
       this.ALT = vc.getAlternateAlleles();
       this.CI = vc.getCommonInfo();
       // construct unique ID for this variant
       List<String> alts = new ArrayList<>();
       for(Allele a : this.ALT) { alts.add(a.getDisplayString()); }
       uniqID = chrom + ":" + this.pos + "_" + this.REF.getDisplayString() + ">" + String.join(",", alts);

       nonRefCalls = new HashMap<>();
       readDepth = new HashMap<>();
       for(Genotype g : vc.getGenotypes()) {
           String patient = g.getSampleName();
           nonRefCalls.put(patient, g);
           readDepth.put(patient, g.getDP());
       }
    }

    public String getID() {
        return(this.uniqID);
    }

    public void addGenotype(Genotype g) {
        nonRefCalls.put(g.getSampleName(), g);
        readDepth.put(g.getSampleName(), g.getDP());
    }

    public ArrayList<String> getAllSampleIDs() {
        ArrayList<String> ret = new ArrayList<>();
        for(Genotype g : this.nonRefCalls.values()) {
            ret.add(g.getSampleName());
        }
        return(ret);
    }


    public Genotype getPatientGenotype() {
        Genotype ret = null;
        // When you call this function, it should only contain one genotype.
        if(nonRefCalls.size() > 1) {
            System.err.println("\nERROR!: multiple genotypes found for" + this.getID() + "\n");
            System.exit(3);
        }

        // In theory you don't need a for loop since it's only one element but this code is easier to read
        for(Genotype g : nonRefCalls.values()) {
            ret = g;
        }
        return(ret);
    }

    public Genotype getPatientGenotype(String k) {
        Genotype ret = null;
        return(nonRefCalls.get(k));
    }

    // Function will create a genotype object populated with only reference calls
    public Genotype createReferenceGenotype(String targetSample) {

        // get 1 patient ID that has this variant
        Map.Entry<String, Genotype> entry = nonRefCalls.entrySet().iterator().next();
        String k = entry.getKey();
        Genotype g = entry.getValue();

        List<Allele> ref_alleles = new ArrayList<>();
        ref_alleles.add(REF);
        ref_alleles.add(REF);

        Genotype ret = GenotypeBuilder.create(targetSample, ref_alleles);
        return(ret);
    }


    // This function makes a new AC, AN, and AF value to go into the INFO column of the VCF
    public String buildNewInfoString(TreeSet<String> samples) {
        String ret = "";
        double AC = 0;
        double AN = 0;
        double AF = 0;

        for(String curSample : samples) {
            Genotype g = null;
            if(this.nonRefCalls.containsKey(curSample)) {
                g = this.nonRefCalls.get(curSample);
                if(g.isHet()) { AC += 1; }
                else if(g.isHomVar()) { AC += 2; }
            }
        }

        AN = samples.size() * 2.0;
        AF = AC / AN;

        String AC_str = "AC=" + String.valueOf(AC);
        String AN_str = "AN=" + String.valueOf(AN);
        String AF_str = "AF=" + String.valueOf(AF);

        ret = AC_str + ";" + AN_str + ";" + AF_str;
        return(ret);
    }


    public void printCalls(TreeSet<String> samples) {
        // construct the ALT allele string. Comma separate multiple alleles
        ArrayList<String> alt_alleles = new ArrayList<>();
        for(Allele a : this.ALT) {
            alt_alleles.add(a.getDisplayString());
        }

        // construct the INFO field string from CI object
        ArrayList<String> info_list = new ArrayList<>();
        for(Map.Entry<String, Object> e : CI.getAttributes().entrySet()) {
            String v = "";
            Object o = e.getValue();
            if(o instanceof ArrayList) {
                ArrayList<Object> OL = (ArrayList<Object>) o;
                ArrayList<String> tmp = new ArrayList<>();
                for(int i = 0; i < OL.size(); i++) {
                    tmp.add((String) OL.get(i));
                }
                v = String.join(",", tmp);
            }
            else {
                v = (String) e.getValue();
            }

            String pair = e.getKey() + "=" + v;
            info_list.add(pair);
        }

        // Get FILTER string
        String filterValue = "PASS";
        for(String s : CI.getFilters()) {
            filterValue = s;
        }

        DecimalFormat df = new DecimalFormat("#.##");
        //String QUALscore = df.format(CI.getPhredScaledQual());

        ArrayList<String> output = new ArrayList<>();
        output.add(this.chrom); //#CHROM
        output.add(String.valueOf(this.pos)); //POS
        output.add("."); // ID field
        output.add(this.REF.getDisplayString()); // REF sequence
        output.add(String.join(",", alt_alleles)); // ALT sequence
        output.add("."); // QUAL
        output.add(filterValue); // FILTER
        output.add(buildNewInfoString(samples)); // INFO

        output.add("GT:AD:DP"); // FORMAT string

        for(String curSample : samples) {
            Genotype g = null;
            String AD = "-1,-1";
            if(this.nonRefCalls.containsKey(curSample)) {
                g = this.nonRefCalls.get(curSample);
            } else {
                g = this.createReferenceGenotype(curSample);
            }
            // construct the column string for this sample
            if(null != g.getAD()) {
                ArrayList<String> adList = new ArrayList<>();
                for (int i = 0; i < g.getAD().length; i++) {
                    adList.add(String.valueOf(g.getAD()[i]));
                }
                AD = String.join(",", adList);
            }

            int DP = g.getDP();
            if(DP == -1) { // compute an estimate for the read depth for the reference calls
                int readSum = 0;
                int N = 0;
                for(int v : this.readDepth.values()) {
                    readSum += v;
                    N++;
                }
                DP = readSum/N;
            }

            String genotype = "0/0";
            if( g.isNoCall() ) genotype = "./.";
            if(g.isHet()) genotype = "0/1";
            if(g.isHomVar()) genotype = "1/1";

            String genoStr = genotype + ":" + AD + ":" + String.valueOf( DP );
            output.add(genoStr);
        }
        String line = (String.join("\t", output));
        System.out.println(line);
        output = null;
    }
}
