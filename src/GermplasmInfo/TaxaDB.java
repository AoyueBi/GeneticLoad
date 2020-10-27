package GermplasmInfo;

/**
 * @author AoyueBi
 * @data 2020-06-22 10:37
 */

import AoUtils.AoFile;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/**
 * Utilities to record the taxa info and invoke any feature quickly
 */
public class TaxaDB {
    private static String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
    //************** 建立对应的 HashMap *********************//
    private static HashMap<String, String> hmTaxaID = null;
    private static HashMap<String, String> hmTaxaMeanDepth = null;
    private static HashMap<String, String> hmTaxaGenomeType = null;
    private static HashMap<String, String> hmTaxaHeterozygousProportion = null;
    private static HashMap<String, String> hmTaxaMissingRate = null;
    private static HashMap<String, String> hmDxy_geneticDivergency = null;
    private static HashMap<String, String> hmTaxaCountry = null;
    private static HashMap<String, String> hmTaxaPart_Continent = null;
    private static HashMap<String, String> hmTaxaContinent = null;
    private static HashMap<String, String> hmTaxaContinent_forTree = null;
    private static HashMap<String, String> hmTaxaIndex_forTree = null;
    private static HashMap<String, String> hmTaxaIndexforMutationBurden = null;
    private static HashMap<String, String> hmTaxaPCA_group = null;
    private static HashMap<String, String> hmTaxaTreeValidatedGroupbyPloidy = null;
    private static HashMap<String, String> hmTaxaTreeValidatedGroupbySubspecies = null;
    private static HashMap<String, String> hmTaxaContinent_by7 = null;

    public List<String> taxaList = new ArrayList<>();
    List<String> columnList = new ArrayList<>(); //查看某一列含有多少个set


    public TaxaDB(){
//        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt");
        this.readDBfile();
    }

    public void readDBfile(){
        //************** 初始化 HashMap *********************//
        hmTaxaID = new HashMap<>();
        hmTaxaMeanDepth = new HashMap<>();
        hmTaxaGenomeType = new HashMap<>();
        hmTaxaHeterozygousProportion = new HashMap<>();
        hmTaxaMissingRate = new HashMap<>();
        hmDxy_geneticDivergency = new HashMap<>();
        hmTaxaCountry = new HashMap<>();
        hmTaxaPart_Continent = new HashMap<>();
        hmTaxaContinent = new HashMap<>();
        hmTaxaContinent_forTree = new HashMap<>();
        hmTaxaIndex_forTree = new HashMap<>();
        hmTaxaIndexforMutationBurden = new HashMap<>();
        hmTaxaPCA_group = new HashMap<>();
        hmTaxaTreeValidatedGroupbyPloidy = new HashMap<>();
        hmTaxaTreeValidatedGroupbySubspecies = new HashMap<>();
        hmTaxaContinent_by7 = new HashMap<>();

        try {
            BufferedReader br = AoFile.readFile(infileS);
            String header = br.readLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String taxa = l.get(0);
                taxaList.add(taxa);
                //**************** add element *******************//
                hmTaxaID.put(taxa,l.get(1));
                hmTaxaMeanDepth.put(taxa,l.get(2));
                hmTaxaGenomeType.put(taxa,l.get(3));
                hmTaxaHeterozygousProportion.put(taxa,l.get(4));
                hmTaxaMissingRate.put(taxa,l.get(5));
               hmDxy_geneticDivergency.put(taxa,l.get(6));
                hmTaxaCountry.put(taxa,l.get(7));
                hmTaxaPart_Continent.put(taxa,l.get(8));
                hmTaxaContinent.put(taxa,l.get(9));
                hmTaxaContinent_forTree.put(taxa,l.get(10));
                hmTaxaIndex_forTree.put(taxa,l.get(11));
                hmTaxaIndexforMutationBurden.put(taxa,l.get(12));
                hmTaxaPCA_group.put(taxa,l.get(13));
                hmTaxaTreeValidatedGroupbyPloidy.put(taxa,l.get(14));
                hmTaxaTreeValidatedGroupbySubspecies.put(taxa,l.get(15));
                hmTaxaContinent_by7.put(taxa,l.get(16));
            }
            br.close();
            System.out.println("There is totally " + taxaList.size() + " taxa in " + infileS);
            System.out.println("********************* Finished reading taxaList files ***************************");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

//**************************** method for taxa feature ***********************************************//
    public String getID(String taxa){return hmTaxaID.get(taxa);}
    public String getMeanDepth(String taxa){return hmTaxaMeanDepth.get(taxa);}
    public String getGenomeType(String taxa){return hmTaxaGenomeType.get(taxa);}
    public String getHeterozygousProportion(String taxa){return hmTaxaHeterozygousProportion.get(taxa);}
    public String getMissingRate(String taxa){return hmTaxaMissingRate.get(taxa);}
    public String getDxy_geneticDivergency(String taxa){return hmDxy_geneticDivergency.get(taxa);}
    public String getCountry(String taxa){return hmTaxaCountry.get(taxa);}
    public String getPart_Continent(String taxa){return hmTaxaPart_Continent.get(taxa);}
    public String getContinent(String taxa){return hmTaxaContinent.get(taxa);}
    public String getContinent_forTree(String taxa){return hmTaxaContinent_forTree.get(taxa);}
    public String getIndex_forTree(String taxa){return hmTaxaIndex_forTree.get(taxa);}
    public String getIndexforMutationBurden(String taxa){return hmTaxaIndexforMutationBurden.get(taxa);}
    public String getPCA_group(String taxa){return hmTaxaPCA_group.get(taxa);}
    public String getTreeValidatedGroupbyPloidy(String taxa){return hmTaxaTreeValidatedGroupbyPloidy.get(taxa);}
    public String getTreeValidatedGroupbySubspecies(String taxa){return hmTaxaTreeValidatedGroupbySubspecies.get(taxa);}
    public String getContinent_by7(String taxa){return hmTaxaContinent_by7.get(taxa);}

    /**
     * get header list from the db file
     * @return
     */
    public static List<String> checkHeaderList(){ return AoFile.getheader(infileS);}

    /**
     *
     * get String list from column index
     * @param columnIndex
     * @return
     */
    public static List<String> getStringList( int columnIndex){
        List<String> out = new ArrayList<>();
        try {
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            String temp = br.readLine(); //read header
            List<String> l = new ArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String goal = l.get(columnIndex);
                out.add(goal);
                cnt++;
            }
            br.close();
            System.out.println("Total num in the list is    " + cnt + "\t" + out.size());
            Collections.sort(out);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return out;
    }

    public void mkcycle(){
        List<String> l = AoFile.getheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt");
        for (int i = 1; i < l.size(); i++) {
            String columnName = l.get(i);
            System.out.println("public String get" + columnName + "(String taxa){return hmTaxa" + columnName + ".get(taxa);}");
        }

    }

}
