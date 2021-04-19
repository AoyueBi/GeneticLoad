package AoUtils.Triads;

import AoUtils.AoFile;

import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import pgl.infra.anno.gene.GeneFeature;

/**
 * Utilities to process triadsList gene of wheat reference genome v1.0 (IWGSC v1.0)
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
    public List<String> triadsList = new ArrayList<>();
    static List<String> genesList = new ArrayList<>();
    static List<String> geneAList = new ArrayList<>();
    static List<String> geneBList = new ArrayList<>();
    static List<String> geneDList = new ArrayList<>();


    //        GeneDB genedb = new GeneDB(); //需要修改
    //        Triadsgenes tg = new Triadsgenes();  //需要修改

    public Triadsgenes(){
        this.readDBfile();
    }

    private void readDBfile(){
        triadsGeneAMap = new HashMap<>();
        triadsGeneBMap = new HashMap<>();
        triadsGeneDMap = new HashMap<>();
        geneTriadsMap = new HashMap<>();
        triadsSyntenicMap = new HashMap<>();
        triadsExpressedMap = new HashMap<>();
//        String infileS = "/Users/Aoyue/Documents/Data/wheat/gene/homoeologsGene111/triadGenes1.1.txt";
//        String infileS = "/Users/Aoyue/Documents/Data/wheat/gene/homoeologsGene111/triadGene1.1_V2.txt";
        String infileS = "/Users/Aoyue/Documents/Data/wheat/gene/homoeologsGene111/triadGenes1.1_V3.txt";

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
                triadsList.add(triadID);
                geneAList.add(genea);
                geneBList.add(geneb);
                geneDList.add(gened);
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
            genesList.addAll(geneAList);genesList.addAll(geneBList);genesList.addAll(geneDList);
            Collections.sort(genesList);
            System.out.println("There is totally " + genesList.size() + " genes in " + triadsList.size() + " triads");
            System.out.println("Finished reading triadsList files");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     *
     * @return
     */
    public int getTriadNum(){
        return triadsList.size();
    }

    /**
     *
     * @return
     */
    public int getGeneNum(){
        return triadsList.size()*3;
    }

    /**
     * there are
     */
    public void checkGenesNotInPGF(){
        String infileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        GeneFeature gf = new GeneFeature(infileS);
        int a = gf.getGeneNumber();
        System.out.println(a);
//pgf genes -> triads
    }

    /**
     *
     * @param gene
     * @return
     */
    public boolean ifTriads (String gene){
        boolean out = false;
        int index = Collections.binarySearch(genesList,gene);
        if (index>-1) out = true;
//        if (out) System.out.println("This gene is a triad gene");
//        else System.out.println("This gene is not a triad gene");
        return out;
    }

    /**
     *
     * @param gene
     * @return
     */
    public boolean ifExpressedBasedGene(String gene){
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
    public boolean ifExpressedBasedTriadID(String triadID){
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
    public boolean ifSyntenicBasedGene(String gene){
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
    public boolean ifSyntenicBasedTriadID(String triadID){
        boolean out = false;
        String result = triadsSyntenicMap.get(triadID);
        if (result.equals("syntenic"))out=true;
        return out;
    }


    /**
     *
     * @param gene
     */
    public String getTriadID(String gene){
        String out;
        out = geneTriadsMap.get(gene);
        return out;
    }

    /**
     *
     * @param triadID
     * @return
     */
    public String getGeneinAsub(String triadID){
        String out;
        out = triadsGeneAMap.get(triadID);
        return out;
    }

    /**
     *
     * @param triadID
     * @return
     */
    public String getGeneinBsub(String triadID){
        String out;
        out = triadsGeneBMap.get(triadID);
        return out;
    }

    /**
     *
     * @param triadID
     * @return
     */
    public String getGeneinDsub(String triadID){
        String out;
        out = triadsGeneDMap.get(triadID);
        return out;
    }

    /**
     *
     * @param triadID
     * @return
     */
    public List<String> getGenesonTriadsID(String triadID){
        List<String> out = new ArrayList<>();
        out.add(getGeneinAsub(triadID));
        out.add(getGeneinBsub(triadID));
        out.add(getGeneinDsub(triadID));

        return out;
    }


}
