package repopRef;

import htsjdk.variant.variantcontext.*;
import java.text.DecimalFormat;
import java.util.*;

public class VariantClass {
    private String chrom;
    private long pos;
    private String REF;
    private List<Allele> ALT;
    private String alt_str;
    private String uniqID;
    private String ID;
    private String qual;
    private String filter;
    private CommonInfo CI;
    private Genotype geno;
    //private HashMap<String, Genotype> nonRefCalls = null; // k = sampleID, v = 0/1, 1/0, 1/1
    //private HashMap<String, Integer> readDepth = null; // k = sampleID, v = DP for this variant in this sample

    public VariantClass(VariantContext vc) {

       this.chrom = vc.getContig();
       this.pos = vc.getStart();
       this.ID = vc.getID();
       this.REF = vc.getReference().getDisplayString();
       this.ALT = vc.getAlternateAlleles();
       this.qual = String.valueOf( vc.getPhredScaledQual() );
       this.CI = vc.getCommonInfo();

       if(vc.getFilters().isEmpty()) this.filter = "PASS";
       else {
           for(String s : vc.getFilters()) {
               this.filter = s;
               break;
           }
       }

       // construct unique ID for this variant
       List<String> alts = new ArrayList<>();
       for(Allele a : this.ALT) { alts.add(a.getDisplayString()); }
       this.alt_str = String.join(",", alts);
       uniqID = this.chrom + ":" + this.pos;

       for(Genotype g : vc.getGenotypes()) {
           geno = g;
       }
    }

    public VariantClass() {
        this.chrom = null;
        this.pos = -1;
        this.ID = null;
        this.REF = null;
        this.ALT = null;
        this.alt_str = null;
    }

    public String getStep1_line() {
        String ret = this.chrom + "\t" + this.pos + "\t" + this.REF + "\t" + this.alt_str;
        return(ret);
    }

    private String formatCommonInfo() {

        ArrayList<String> ary = new ArrayList<>();
        for(String s : this.CI.getAttributes().keySet()) {
            String pair = s + "=" + this.CI.getAttribute(s);
            ary.add(pair);
        }
        String ret = String.join(";", ary);
        return(ret);
    }

    public String getLine() {
        List<String> ret = new ArrayList<>();

        ret.add(this.chrom);
        ret.add(String.valueOf(this.pos));
        ret.add(this.ID);
        ret.add(this.REF);
        ret.add(alt_str);
        ret.add(this.qual);
        ret.add(this.filter);
        ret.add(formatCommonInfo());
        ret.add("GT:AD:DP");

        String DP_str = String.valueOf( geno.getDP() );
        ArrayList<String> ary = new ArrayList<>();
        for(int i : geno.getAD()) { ary.add(String.valueOf(i));  }
        String AD_str = String.join(",", ary);

        String geno_str = getGenotypeCode() + ":" + AD_str + ":" + DP_str;
        ret.add(geno_str);

        String ret_str = String.join("\t", ret);
        return(ret_str);
    }


    private String getGenotypeCode() {
        String ret = "";
        if(this.geno.isNoCall()) ret = "./.";
        if(this.geno.isHet()) ret = "0/1";
        if(this.geno.isHomVar()) ret = "1/1";
        return(ret);
    }


    // This version of the function is specifically for REFERENCE calls that don't exist in the VCF file
    public String getLine_REF(String chrom, int pos, String ref, String alt) {
       List<String> ret = new ArrayList<>();

       ret.add(chrom);
       ret.add( String.valueOf(pos) );
       ret.add("."); // ID
       ret.add(ref);
       ret.add(alt);
       ret.add("."); // qual
       ret.add("REF_CALL"); // filter
       ret.add("."); // INFO
       ret.add("GT:AD:DP");
       ret.add("0/0:-1,-1:-1");

       String ret_str = String.join("\t", ret);
       return(ret_str);
    }


    public String getID() {
        return(this.uniqID);
    }
/*
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

 */
}
