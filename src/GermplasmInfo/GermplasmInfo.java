/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package GermplasmInfo;

import AoUtils.AoFile;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.*;

/**
 *
 * @author Aoyue
 */
public class GermplasmInfo {

    public GermplasmInfo() {
//        this.addIftaxaonVMap2Anno();
//        this.addTreeValidatedGroupbyPloid();
//        this.addTreeValidatedGroupbySuspecies();
//        this.addInfo();

//        this.getWild_emmer_South2();
//        this.getEuropeanLandrace();

//        this.mergeTxt();

//        this.addColumntoTaxaDB();
//        this.addGroup();
//        this.summaryGroupbyContinent();

    }

    /**
     *
     */
    public void summaryGroupbyContinent(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
        HashMap<String,String> hmPartConti2newConti = new HashMap<>();
        String value = null;
        AoFile.readheader(infileS);
        String[] partContinentArray = AoFile.getStringArraybySet(infileS,8);
        for (int i = 0; i < partContinentArray.length; i++) {
            System.out.println(partContinentArray[i]);
            String key = partContinentArray[i];
            if (key.equals("Africa")){
                value = "Africa";
                hmPartConti2newConti.put(key,value);
            }
            if (key.equals("Central Asia")){
                value = "Central and South Asia";
                hmPartConti2newConti.put(key,value);
            }
            if (key.equals("Central Europe")){
                value = "Europe";
                hmPartConti2newConti.put(key,value);
            }
            if (key.equals("East Asia")){
                value = "East Asia";
                hmPartConti2newConti.put(key,value);
            }
            if (key.equals("Eastern Europe")){
                value = "Europe";
                hmPartConti2newConti.put(key,value);
            }
            if (key.equals("NA")){
                value = "NA";
                hmPartConti2newConti.put(key,value);
            }
            if (key.equals("North America")){
                value = "America";
                hmPartConti2newConti.put(key,value);
            }
            if (key.equals("Northern Europe")){
                value = "Europe";
                hmPartConti2newConti.put(key,value);
            }
            if (key.equals("Oceania")){
                value = "Oceania";
                hmPartConti2newConti.put(key,value);
            }
            if (key.equals("South America")){
                value = "America";
                hmPartConti2newConti.put(key,value);
            }
            if (key.equals("South Asia")){
                value = "Central and South Asia";
                hmPartConti2newConti.put(key,value);
            }
            if (key.equals("Southeast Europe")){
                value = "Europe";
                hmPartConti2newConti.put(key,value);
            }
            if (key.equals("Southern Europe")){
                value = "Europe";
                hmPartConti2newConti.put(key,value);
            }
            if (key.equals("Western Asia")){
                value = "Western Asia";
                hmPartConti2newConti.put(key,value);
            }
            if (key.equals("Western Europe")){
                value = "Europe";
                hmPartConti2newConti.put(key,value);
            }
        }
        System.out.println(partContinentArray.length);

        AoFile.addColumbyString(infileS,8,hmPartConti2newConti,"Continent_by7");

    }

    public void addGroup(){
//        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
//        AoFile.readheader(dbfileS);
//        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/032_pca/001_input/001_matrix.txt";
//        HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,0,3);
//        AoFile.addColumbyString(taxaFileS,0,hm,"GenomeType");

        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
        AoFile.readheader(dbfileS);
//        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/032_pca/001_input/002_matrix_Asub.txt";
//        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/032_pca/001_input/002_matrix_Bsub.txt";
//        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/032_pca/001_input/002_matrix_Dsub.txt";
//        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/032_pca/001_input/003_matrix_hexa.txt";
//        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/032_pca/001_input/004_matrix_tetra.txt";
        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/032_pca/001_input/005_matrix_DD.txt";
        int[] columnIndexes = {3,8,10,12,15,16};
        HashMap<String,String>[] hm = new AoFile().getHashMapsStringKey(dbfileS,0,columnIndexes);
        AoFile.addColumsbyString(taxaFileS,0,hm,"\tGenomeType\tPart_Continent\tContinent_forTree\tIndexforMutationBurden\tSubspecies\tContinent_by7");

    }

    /**
     * 向新建立的 vmap2.0 taxa info DB 添加各种信息库
     */
    public void addColumntoTaxaDB(){
//        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/002_QC/002_merge/001_taxa_QC.txt";
//        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
////        HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,0,1);
////        AoFile.addColumbyString(taxaFileS,0,hm,"HeterozygousProportion");
//        HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,0,2);
//        AoFile.addColumbyString(taxaFileS,0,hm,"MissingRate");

//        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/007_qualityCheck/004_addInfo/001_taxa_QC.txt.gz";
//        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
//        HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,0,1);
//        AoFile.addColumbyString(taxaFileS,0,hm,"HeterozygousProportion_onHapscanner");

//        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt";
//        AoFile.readheader(dbfileS);
//        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/temp/taxa_InfoDB.txt";
//        HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,0,12);
//        AoFile.addColumbyString(taxaFileS,0,hm,"MeanDepth_onVmap2.1");

//        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/026_depth/001_file/depth_PopDepth_VCF.txt";
//        AoFile.readheader(dbfileS);
//        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/temp/taxa_InfoDB.txt";
//        HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,0,1);
//        AoFile.addColumbyString(taxaFileS,0,hm,"MeanDepth_fromPopdepth");

//        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/031_VMap2.0_QC/005_dxy/dxy.txt";
//        AoFile.readheader(dbfileS);
//        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
//        HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,0,2);
//        AoFile.addColumbyString(taxaFileS,0,hm,"Dxy_geneticDivergency");

        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList_addCS2017.txt";
        AoFile.readheader(dbfileS);
        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
        HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,0,11);
        AoFile.addColumbyString(taxaFileS,0,hm,"TreeValidatedGroupbySubspecies");

    }

    public void mergeTxt(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/013_indiHeter";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/013_indiHeter/merge/heter_indivi.txt";
        AoFile.mergeTxtbysuffix(infileDirS,outfileS,".txt");
    }

    /**
     *  there 12 taxa (European materials) clustered on Asia branch, we need to remove that.
     */
    public void getEuropeanLandrace(){

        String landrace_exclusionEuropean = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/009_treeValidatedFroup_byRegion/002_Landrace_European/Landrace_exclusionEuropean.txt";
        List<String> lonAsia = new AoFile().getStringListwithoutHeader(landrace_exclusionEuropean,0);
        try {
            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/009_treeValidatedFroup_byRegion/002_Landrace_European/Landrace_European.txt";
            new AoFile().readheader(infileS);
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String taxa = l.get(0);
                String continent = l.get(4);
                String subspecies = l.get(11);
                if(subspecies.equals("Landrace")){
                    if(continent.equals("Europe")){
                        int index = Collections.binarySearch(lonAsia,taxa);
                        if (index < 0){
                            bw.write(taxa);
                            bw.newLine();
                        }
                    }
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

    /**
     * based the tree, I divide the wild emmer group into Wild_emmer_North  Wild_emmer_South1  Wild_emmer_South2
     */
    public void getWild_emmer_South2(){

        HashMap<String,String> hm = new AoFile().getHashMapwithFileDirs("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/009_Wild_emmer");
        List<String> we = new ArrayList<>(hm.keySet());
        Collections.sort(we);
        try {
            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/009_Wild_emmer/Wild_emmer_South2.txt";
            new AoFile().readheader(infileS);
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String taxa = l.get(0);
                String subspecies = l.get(11);
                if(subspecies.equals("Wild_emmer")){
                    int index = Collections.binarySearch(we,taxa);
                    if (index < 0){
                        bw.write(taxa);
                        bw.newLine();
                    }
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

    /**
     * add info to taxaList file, the info db file from /Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/wheatVMapII_germplasmInfo_20191225.txt
     */
    public void addInfo(){
//        String dbfileS ="/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/wheatVMapII_germplasmInfo_20191225.txt";
//        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";

        //        HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,4,5);
//        new AoFile().addColumbyString(taxaFileS,0,hm,"Ploidy");

//        HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,4,24);
//        new AoFile().addColumbyString(taxaFileS,0,hm,"PCA_group");

//        HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,4,24);
//        new AoFile().addColumbyString(taxaFileS,0,hm,"PCA_group");

//                HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,4,26);
//        new AoFile().addColumbyString(taxaFileS,0,hm,"TreeValidatedGroupbyPloidy");

//        HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,4,27);
//        new AoFile().addColumbyString(taxaFileS,0,hm,"TreeValidatedGroupbySubspecies");

//        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_summary.txt";
//        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt";

        //添加个体杂合度信息
//        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/013_indiHeter/merge/heter_indivi.txt";
//        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt";
//        HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,0,3);
//        AoFile.addColumbyString(taxaFileS,0,hm,"MeanDepth");
//        AoFile.addColumbyString(taxaFileS,0,hm,"Heterozygosity_Indivi");

        //添加到CS的 genetic divergence
        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/016_dxy/dxy.txt";
        String taxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt";
        HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,0,1);
        AoFile.addColumbyString(taxaFileS,0,hm,"Dxy_geneticDivergency");

    }

    //
    public void addTreeValidatedGroupbySuspecies(){
        String groupDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/008_treeValidatedGroup_bySubspecies";
        HashMap<String,String> hm = new AoFile().getHashMapwithFileDirs(groupDirS);

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
        HashMap<String,String> hm = new AoFile().getHashMapwithFileDirs(groupDirS);

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
