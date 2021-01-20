package ExonAnnotation;

import AoUtils.AoFile;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.*;

public class ExonAnnotation {

    String exonAnnotationFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/001_exonSNP_anno_addGroup.txt.gz";
    String ID;
    int chr;
    int pos;
    String sub;
    String ancestral;
    String derivedSift;
    String gerp;
    String variantsGroup;
    String ifRefisAnc;
    String loadGroup;

    //************** 建立对应的 HashMap *********************//
    private static HashMap<String, String> hmIDVariantsGroup = null;
    private static HashMap<String, String> hmIDIfRefisAnc = null;
    private static HashMap<String, String> hmIDLoadGroup = null;
    private static HashMap<String, String> hmIDSub = null;
    Set<String> variantsGroupSet = new HashSet<>();
    Set<String> loadGroupSet = new HashSet<>();
    List<String> idList = new ArrayList<>();



    public ExonAnnotation(){
        this.readExonAnnotation();
    }

    public void readExonAnnotation(){

        hmIDVariantsGroup = new HashMap<>();
        hmIDIfRefisAnc = new HashMap<>();
        hmIDLoadGroup = new HashMap<>();
        hmIDSub = new HashMap<>();

        try {
            BufferedReader br = AoFile.readFile(exonAnnotationFileS);
            String header = br.readLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                ID = l.get(0);
                chr = Integer.parseInt(l.get(1));
                pos = Integer.parseInt(l.get(2));
                ancestral = l.get(15);
                derivedSift = l.get(16);
                gerp = l.get(20);
                sub = l.get(22);
                variantsGroup = l.get(23);
                ifRefisAnc = l.get(24);
                loadGroup = l.get(25);

                idList.add(ID);
                variantsGroupSet.add(variantsGroup);
                loadGroupSet.add(loadGroup);
                hmIDVariantsGroup.put(ID,variantsGroup);
                hmIDIfRefisAnc.put(ID,ifRefisAnc);
                hmIDLoadGroup.put(ID,loadGroup);
                hmIDSub.put(ID,sub);
            }
            Collections.sort(idList);


            br.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public String getVariantsGroup(String ID){return hmIDVariantsGroup.get(ID);}
    public String getIfRefisAnc(String ID){return hmIDIfRefisAnc.get(ID);}
    public String getLoadGroup(String ID){return hmIDLoadGroup.get(ID);}
    public String getSub(String ID){return hmIDSub.get(ID);}

    public String[] variantsGroupArray(){return variantsGroupSet.toArray(new String[variantsGroupSet.size()]);}
    public String[] loadGroupArray(){return loadGroupSet.toArray(new String[loadGroupSet.size()]);}


    public boolean ifinExonAnnotation(String ID){
        int index = Collections.binarySearch(idList,ID);
        if (index>-1) {
            return true;
        }else {
            return false;
        }
    }
}
