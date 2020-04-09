package AoUtils.Triads;

import AoUtils.AoFile;
import pgl.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/**
 * Utilities to process triads gene of wheat reference genome v1.0 (IWGSC v1.0)
 * Other functions may be provided later
 * @author AoyueBi
 */

public class Triadsgenes {
    //TriadID	A	B	D	synteny	Expressed
    //T000001	TraesCS7A02G243100	TraesCS7B02G148400	TraesCS7D02G241900	syntenic	TRUE
    //T000002	TraesCS7A02G360600	TraesCS7B02G267100	TraesCS7D02G362400	syntenic	TRUE

    private static HashMap<String, String> triadsGeneAMap = null;
    private static HashMap<String, String> triadsGeneBMap = null;
    private static HashMap<String, String> triadsGeneDMap = null;
    private static HashMap<String, String> geneTriadsMap = null;
    private static HashMap<String, String> triadsSyntenicMap = null;
    private static HashMap<String, String> triadsExpressedMap = null;
    static List<String> triads = new ArrayList<>();
    static List<String> geneA = new ArrayList<>();
    static List<String> geneB = new ArrayList<>();
    static List<String> geneD = new ArrayList<>();

    private static void readDBfile(){
        triadsGeneAMap = new HashMap<>();
        triadsGeneBMap = new HashMap<>();
        triadsGeneDMap = new HashMap<>();
        geneTriadsMap = new HashMap<>();
        triadsSyntenicMap = new HashMap<>();
        triadsExpressedMap = new HashMap<>();
        String infileS = "/Users/Aoyue/Documents/Data/wheat/gene/homoeologsGene111/triadGenes1.1.txt";
        try {
            BufferedReader br = AoFile.readFile(infileS);
            String header = br.readLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String triadID = l.get(0);
                String genea = l.get(1);
                String geneb = l.get(2);
                String gened = l.get(3);
                String synteny = l.get(4);
                String expressed = l.get(5);
                triads.add(triadID);
                geneA.add(genea);
                geneB.add(geneb);
                geneD.add(gened);
                triadsGeneAMap.put(triadID,genea);
                triadsGeneBMap.put(triadID,geneb);
                triadsGeneDMap.put(triadID,gened);
                geneTriadsMap.put(genea,triadID);
                geneTriadsMap.put(geneb,triadID);
                geneTriadsMap.put(gened,triadID);
                triadsSyntenicMap.put(triadID,synteny);
                triadsExpressedMap.put(triadID,expressed);
            }
            br.close();
            System.out.println("Finished reading triads files");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public static boolean ifTriads (String gene){
        boolean out = false;
        readDBfile();
        List<String> l = new ArrayList<>();
        l.addAll(geneA);
        l.addAll(geneB);
        l.addAll(geneD);
        Collections.sort(l);
        int index = Collections.binarySearch(l,gene);
        if (index>-1) out = true;
        if (out) System.out.println("This gene is a triad gene");
        else System.out.println("This gene is not a triad gene");
        return out;
    }

    /**
     *
     * @param gene
     * @return
     */
    public static boolean ifExpressedBasedGene(String gene){
        readDBfile();
        boolean out = false;
        String triadID = getTriadID(gene);
        String result = triadsExpressedMap.get(triadID);
        if (result.equals("TRUE"))out=true; //FALSE
        return out;
    }

    /**
     *
     * @param triadID
     * @return
     */
    public static boolean ifExpressedBasedTriadID(String triadID){
        readDBfile();
        boolean out = false;
        String result = triadsExpressedMap.get(triadID);
        if (result.equals("TRUE"))out=true;
        return out;
    }

    /**
     *
     * @param gene
     * @return
     */
    public static boolean ifSyntenicBasedGene(String gene){
        readDBfile();
        boolean out = false;
        String triadID = getTriadID(gene);
        String result = triadsSyntenicMap.get(triadID);
        if (result.equals("syntenic"))out=true; //non-syntenic
        return out;
    }

    /**
     *
     * @param triadID
     * @return
     */
    public static boolean ifSyntenicBasedTriadID(String triadID){
        readDBfile();
        boolean out = false;
        String result = triadsSyntenicMap.get(triadID);
        if (result.equals("syntenic"))out=true;
        return out;
    }


    /**
     *
     * @param gene
     */
    public static String getTriadID(String gene){
        readDBfile();
        String out;
        out = geneTriadsMap.get(gene);
        return out;
    }

    /**
     *
     * @param triadID
     * @return
     */
    public static String getGeneinAsub(String triadID){
        readDBfile();
        String out;
        out = triadsGeneAMap.get(triadID);
        return out;
    }

    /**
     *
     * @param triadID
     * @return
     */
    public static String getGeneinBsub(String triadID){
        readDBfile();
        String out;
        out = triadsGeneBMap.get(triadID);
        return out;
    }

    /**
     *
     * @param triadID
     * @return
     */
    public static String getGeneinDsub(String triadID){
        readDBfile();
        String out;
        out = triadsGeneDMap.get(triadID);
        return out;
    }

    public static List<String> getGenesonTriadsID(String triadID){
        List<String> out = new ArrayList<>();
        out.add(getGeneinAsub(triadID));
        out.add(getGeneinBsub(triadID));
        out.add(getGeneinDsub(triadID));

        return out;
    }


}
