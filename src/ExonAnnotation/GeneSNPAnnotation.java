package ExonAnnotation;

import AoUtils.AoFile;
import gnu.trove.list.array.TCharArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.lang.StringUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.util.*;

public class GeneSNPAnnotation { //## 注意：请手动选择如何定义有害突变。

    String exonAnnotationFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/006_geneSNPAnnotation_merge/011_geneSNPAnno.txt.gz";
    String ID;
    String Chr;
    String Pos;
    String Ref;
    String Alt;
    String Major;
    String Minor;
    String Maf;
    String Transcript;
    String Ancestral;
    String DAF;
    String gerp;
    String Region;
    String variantType;
    String Alt_SIFT;
    String Ref_SIFT;
    String sift;
    String Effect_VEP;
    String Impact_VEP;
    String Effect_snpEff;
    String Impact_snpEff;
    String Rank_Total;
    String gerp16way;
    String AMINO_POS;
    String UniProtKB_AC;
    String LIST_S2;
    String phyloP;
    String Amino_Size;
    String phyloP_RefMask;
    String AlleleAge_J;
    String AlleleAge_M;
    String AlleleAge_R;
    String Alt_PolyPhen2_prediction;
    String Alt_PolyPhen2_class;
    String altPolyPhen;
    String Ref_PolyPhen2_prediction;
    String Ref_PolyPhen2_class;
    String Ref_PolyPhen2_prob;
    String Derived_PolyPhen2_prediction;
    String Derived_PolyPhen2_class;
    String derivedPolyPhen;
    ////// 以下需要在程序中计算
    String variantsGroup;
    String loadGroup;
    String ifRefisAnc;
    String sub;
    int syn = 0;

    //************** 建立对应的 HashMap *********************//
    private static HashMap<String, String> hmIDVariantsGroup = null;
    private static HashMap<String, String> hmIDIfRefisAnc = null;
    private static HashMap<String, String> hmIDLoadGroup = null;
    private static HashMap<String, String> hmIDSub = null;
    Set<String> variantsGroupSet = new HashSet<>();
    Set<String> loadGroupSet = new HashSet<>();
    List<String> idList = new ArrayList<>();

    String variantType1 = "001_synonymous";
    String variantType2 = "002_nonsynonymous";
    String variantType3 = "003_nonsynGERPandDerivedSIFT";
    String variantType4 = "004_nonsynDerivedSIFT";
    String variantType5 = "005_GERP";
    String variantType6 = "006_SIFT_StopGain"; //只是SIFT中为
    String variantType7 = "007_VEP";
    String variantType8 = "008_snpEff";
    String variantType9 = "009_VEP_stopGained";
    String variantType10 = "010_GERP16way_1"; //1
    String variantType11 = "011_GERP16way_1.2_max";
    String variantType12 = "012_GERP16way_1.2_2.5";
    String variantType13 = "013_GERP16way_2.5_max";
    String variantType14 = "014_GERP16way_2.14_max";
    String variantType15 = "015_LIST_S2";
    String variantType16 = "016_PhyloP1.5_RefNotMask";
    String variantType17 = "017_PhyloP1.5_RefMask";
    String variantType18 = "018_GERP16way_2.14andSIFT"; //2.14
    String variantType19 = "019_Alt_PolyPhen2";
    String variantType20 = "020_Derived_PolyPhen2";
    String variantType21 = "021_GERP16way_2.14andPolyPhen2"; //2.14
    String variantType22 = "022_GERP16way_1.78";
    String variantType23 = "023_GERP16way_1.78andPolyPhen2";
    String variantType24 = "024_GERP16way_1andPolyPhen2";
    String variantType25 = "025_GERP16way_3.1";
    String variantType26 = "026_Derived_PolyPhen2_probably";
    String variantType27 = "027_GERP16way_1andDerived_PolyPhen2_probably";
    String variantType28 = "028_GERP16way_1.78andDerived_PolyPhen2_probably";
    String variantType29 = "029_GERP16way_2.14andDerived_PolyPhen2_probably";
    String variantType30 = "030_GERP16way_3.1andDerived_PolyPhen2_probably";
    String variantType31 = "031_GERP16way_1.5andDerived_PolyPhen2_0.5";


    public GeneSNPAnnotation(String exonAnnotationFileS) {
//        String type = variantType14;
        String type = variantType31;
        this.exonAnnotationFileS=exonAnnotationFileS;
        this.readExonAnnotation(type);
    }

    public GeneSNPAnnotation() {
//        String type = variantType14;
        String type = variantType31;
        this.readExonAnnotation(type);
    }


    public void readExonAnnotation(String type) { //## 注意：请手动选择如何定义有害突变。
        hmIDVariantsGroup = new HashMap<>();
        hmIDIfRefisAnc = new HashMap<>();
        hmIDLoadGroup = new HashMap<>();
        hmIDSub = new HashMap<>();


        int chrNum = 42;
        String[] subs = {"A", "A", "B", "B", "D", "D", "A", "A", "B", "B", "D", "D", "A", "A", "B", "B", "D", "D", "A", "A", "B", "B", "D", "D", "A", "A", "B", "B", "D", "D", "A", "A", "B", "B", "D", "D", "A", "A", "B", "B", "D", "D"};
        HashMap<Integer, String> hmchr2sub = new HashMap<>();
        for (int i = 0; i < 42; i++) {
            hmchr2sub.put(i + 1, subs[i]);
        }

        try {
            BufferedReader br = AoFile.readFile(exonAnnotationFileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++;
                l = PStringUtils.fastSplit(temp);
                ID = l.get(0);
                Chr = l.get(1);
                Pos = l.get(2);
                Ref = l.get(3);
                Alt = l.get(4);
                Major = l.get(5);
                Minor = l.get(6);
                Maf = l.get(7);
                Transcript = l.get(8);
                Ancestral = l.get(9);
                DAF = l.get(10);
                gerp = l.get(11);
                Region = l.get(12);
                variantType = l.get(13);
                Alt_SIFT = l.get(14);
                Ref_SIFT = l.get(15);
                sift = l.get(16);
                Effect_VEP = l.get(17);
                Impact_VEP = l.get(18);
                Effect_snpEff = l.get(19);
                Impact_snpEff = l.get(20);
                Rank_Total = l.get(21);
                gerp16way = l.get(22);
                AMINO_POS = l.get(23);
                UniProtKB_AC = l.get(24);
                LIST_S2 = l.get(25);
                phyloP = l.get(26);
                Amino_Size = l.get(27);
                phyloP_RefMask = l.get(28);
                AlleleAge_J = l.get(29);
                AlleleAge_M = l.get(30);
                AlleleAge_R = l.get(31);
                Alt_PolyPhen2_prediction = l.get(32);
                Alt_PolyPhen2_class = l.get(33);
                altPolyPhen = l.get(34);
                Ref_PolyPhen2_prediction = l.get(35);
                Ref_PolyPhen2_class = l.get(36);
                Ref_PolyPhen2_prob = l.get(37);
                Derived_PolyPhen2_prediction = l.get(38);
                Derived_PolyPhen2_class = l.get(39);
                derivedPolyPhen = l.get(40);
//************************************************* 第一步： 全基因组水平上的分组判断 *****************************************************************************//
                sub = hmchr2sub.get(Integer.parseInt(Chr));
                hmIDSub.put(ID, sub);
                idList.add(ID);

                if (Ref.equals(Ancestral)) {
                    ifRefisAnc = "Anc";
                    hmIDIfRefisAnc.put(ID, ifRefisAnc);
                }
                if (Alt.equals(Ancestral)) {
                    ifRefisAnc = "Der";
                    hmIDIfRefisAnc.put(ID, ifRefisAnc);
                }
                if (!Ancestral.equals(Ref) && !Ancestral.equals(Alt)) {
                    ifRefisAnc = "NA";
                    hmIDIfRefisAnc.put(ID, ifRefisAnc);
                }

                if (!Region.equals("CDS")) {
                    if (Region.equals("Intron")) {
                        variantsGroup = "Intron";
                        hmIDVariantsGroup.put(ID, variantsGroup);
                    }
                    if (Region.equals("UTR_3")) {
                        variantsGroup = "UTR_3";
                        hmIDVariantsGroup.put(ID, variantsGroup);
                    }
                    if (Region.equals("UTR_5")) {
                        variantsGroup = "UTR_5";
                        hmIDVariantsGroup.put(ID, variantsGroup);
                    }
                    variantsGroupSet.add(variantsGroup);

                } else if (Region.equals("CDS")) {
                    if (variantType.equals("NONSYNONYMOUS")) {
                        variantsGroup = "Nonsyn";
                        hmIDVariantsGroup.put(ID, variantsGroup);
                    }
                    if (variantType.equals("SYNONYMOUS")) {
                        variantsGroup = "Syn";
                        hmIDVariantsGroup.put(ID, variantsGroup);
                        if (sub.equals("D")) {
                            syn++;
                        }
                    }
                    if (variantType.equals("STOP-GAIN")) {
                        variantsGroup = "Nonsense";
                        hmIDVariantsGroup.put(ID, variantsGroup);
                    }
                    if (Impact_VEP.equals("HIGH") && Effect_VEP.contains("splice_")) {
                        variantsGroup = "Essential splice";
                        hmIDVariantsGroup.put(ID, variantsGroup);
                    }
                    variantsGroupSet.add(variantsGroup);
                }

//************************************************* 第二步： 全基因组水平上含有祖先等位位点的分组判断 *****************************************************************************//
                if (!ifRefisAnc.equals("NA")) {
                    if (variantType.equals("SYNONYMOUS")) {
                        loadGroup = "Synonymous";
                        hmIDLoadGroup.put(ID, loadGroup);
                        loadGroupSet.add(loadGroup);
                    }
                    if (variantType.equals("NONSYNONYMOUS")) {
                        loadGroup = "Nonsynonymous";
                        hmIDLoadGroup.put(ID, loadGroup);
                        loadGroupSet.add(loadGroup);
                    }
                    if (type.equals("014_GERP16way_2.14_max")) {
                        if (variantType.equals("NONSYNONYMOUS")) {
                            if (!gerp16way.startsWith("N")) {
                                double gerp16wayd = Double.parseDouble(gerp16way);
                                if (gerp16wayd >= 2.14) {
                                    //说明必须满足gerp大于等于 2.14
                                    loadGroup = "Deleterious";
                                    hmIDLoadGroup.put(ID, loadGroup);
                                    loadGroupSet.add(loadGroup);
                                }
                            }
                        }
                    }
                    if (type.equals("031_GERP16way_1.5andDerived_PolyPhen2_0.5")) {
                        if (variantType.equals("NONSYNONYMOUS")) {
                            if (!gerp16way.startsWith("N") && (!derivedPolyPhen.startsWith("N"))) {
                                double gerp16wayd = Double.parseDouble(gerp16way);
                                double derivedPolyPhend = Double.parseDouble(derivedPolyPhen);
                                if (gerp16wayd >= 1.5 && derivedPolyPhend >= 0.5) { //这里注意判断条件一定要认真核对！！！
                                    loadGroup = "Deleterious";
                                    hmIDLoadGroup.put(ID, loadGroup);
                                    loadGroupSet.add(loadGroup);
                                }
                            }
                        }
                    }
                }
            }
            br.close();
            Collections.sort(idList); //一定要排序，不然不能进行搜索功能
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    public String getVariantsGroup(String ID) {
        return hmIDVariantsGroup.get(ID);
    }

    public String getIfRefisAnc(String ID) {
        return hmIDIfRefisAnc.get(ID);
    }

    public String getLoadGroup(String ID) {
        return hmIDLoadGroup.get(ID);
    }

    public Set<String> getLoadGroupSet() {
        return loadGroupSet;
    }

    public HashMap<String, String> gethmIDLoadgroup() {
        return hmIDLoadGroup;
    }

    public String getSub(String ID) {
        return hmIDSub.get(ID);
    }


    public String[] variantsGroupArray() {
        return variantsGroupSet.toArray(new String[variantsGroupSet.size()]);
    }

    public String[] loadGroupArray() {
        return loadGroupSet.toArray(new String[loadGroupSet.size()]);
    }


    public boolean ifinExonAnnotation(String ID) {
        int index = Collections.binarySearch(idList, ID);
        if (index > -1) {
            return true;
        } else {
            return false;
        }
    }

    public Integer getSyn() {
        return syn;
    }

    public Integer getIDlength() {
        return idList.size();
    }
}
