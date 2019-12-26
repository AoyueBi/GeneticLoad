/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package GermplasmInfo;

import AoUtils.AoFile;
import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/**
 *
 * @author Aoyue
 */
public class GermplasmInfo {

    public GermplasmInfo() {
//        this.addIftaxaonVMap2Anno();
//        this.addTreeValidatedGroupbyPloid();
        this.addTreeValidatedGroupbySuspecies();

    }

    //
    public void addTreeValidatedGroupbySuspecies(){
        String groupDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/008_treeValidatedGroup_bySubspecies";
        HashMap<String,String> hm = new AoFile().getHashMapwithFilename(groupDirS);

        String outfileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/004_wheatVMapII_germplasmInfoAddSubspecies_20191225.txt";
        String infileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/003_wheatVMapII_germplasmInfoAddPloidy_20191225.txt";

        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            bw.write(temp + "\t" + "TreeValidatedGroupbySubspecies");
            bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String taxa = l.get(4);
                cnt++;
                String value = hm.get(taxa);
                if (value == null || value == "") { //indicated there is no value
                    bw.write(temp + "\t" + "NA");
                    bw.newLine();
                } else {//indicated that id for VMap2
                    bw.write(temp + "\t" + value);
                    bw.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
            new AoFile().readheader(outfileS);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }


    }

    /**
     * 四倍体除去波斯(9个)不分析，剩余四倍体合计 187-9=178个
     * 六倍体除去新疆（5个）和spelt（14个）不分析，剩余六倍体合计419-19=400个
     *
     * 按照倍性区分：
     * Hexaploid
     * ExclusionHexaploid
     * Tetraploid
     * ExclusionTetraploid
     */
    public void addTreeValidatedGroupbyPloid() {
        String groupDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/008_treeValidatedGroup_byPloidy";
        HashMap<String,String> hm = new AoFile().getHashMapwithFilename(groupDirS);

        String infileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/002_wheatVMapII_germplasmInfo_20191225.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/003_wheatVMapII_germplasmInfoAddPloidy_20191225.txt";
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            bw.write(temp + "\t" + "TreeValidatedGroupbyPloidy");
            bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String taxa = l.get(4);
                cnt++;
                String value = hm.get(taxa);
                if (value == null || value == "") { //indicated there is no value
                    bw.write(temp + "\t" + "NA");
                    bw.newLine();
                } else {//indicated that id for VMap2
                    bw.write(temp + "\t" + value);
                    bw.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }


    }

    public void addIftaxaonVMap2Anno() {
        //在ABD六倍体中，去除 IG140057 PI583718 PI534284 Beaqle Caruton
        //在AB四倍体中，去除 PI466930 PI466959 PI272522
        //在D二倍体中，去除 Ku-2071 TA2462 AE430
        String[] removeTaxa = {"IG140057", "PI583718", "PI534284", "Beaqle", "Caruton", "PI466930", "PI466959_2", "PI272522", "KU-2071", "TA2462", "AE430"};
        Arrays.sort(removeTaxa);
        String infileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/001_wheatVMapII_germplasmInfo_20191206.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/002_wheatVMapII_germplasmInfo_20191225.txt";
        new AoFile().readheader(infileS);
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            bw.write(temp + "\t" + "IfOnVmap2");
            bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String taxa = l.get(4);
                cnt++;
                int index = Arrays.binarySearch(removeTaxa, taxa);
                if (index > -1) { //搜索到了，说明是被删除的样本
                    bw.write(temp + "\t" + "0");
                    bw.newLine();
                } else {//index <0, indicated that id for VMap2
                    bw.write(temp + "\t" + "1");
                    bw.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

}
